#!/usr/bin/env python3
"""Build a MinKNOW playback bulk FAST5 from read-oriented POD5 files.

This is intentionally a two-stage tool. ``make-background`` extracts short,
verified open-pore intervals from a real bulk FAST5 once. ``convert`` can then
reuse that small cache with any compatible POD5 run.

Only POD5 metadata are scanned during indexing. Signal is loaded one channel at
a time while the output is written, so the complete input signal is never held
in memory.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import glob
import json
import math
import os
import shlex
import sqlite3
import sys
import uuid
import warnings
from contextlib import ExitStack
from pathlib import Path
from typing import Iterable, Sequence

import h5py
import numpy as np
import pod5
import vbz_h5py_plugin

vbz_h5py_plugin.register_plugin()
warnings.filterwarnings(
    "ignore",
    message=r"Call to deprecated function .*Scaling fields were unused.*",
    category=DeprecationWarning,
)

STATE_PORE = 5
STATE_STRAND = 7
STATE_UNCLASSIFIED = 200
READ_CLASS_STRAND = ord("S")
END_REASON = {
    "unknown": 0,
    "partial": 1,
    "mux_change": 2,
    "unblock_mux_change": 3,
    "data_service_unblock_mux_change": 4,
    "signal_positive": 5,
    "signal_negative": 6,
}


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def expand_inputs(values: Sequence[str]) -> list[Path]:
    found: list[Path] = []
    for value in values:
        path = Path(value)
        if path.is_dir():
            found.extend(sorted(path.glob("*.pod5")))
        elif any(c in value for c in "*?["):
            found.extend(Path(x) for x in sorted(glob.glob(value)))
        elif path.is_file():
            found.append(path)
        else:
            raise FileNotFoundError(value)
    unique = list(dict.fromkeys(x.resolve() for x in found))
    if not unique:
        raise ValueError("No POD5 files found")
    return unique


def set_string_attrs(group: h5py.Group, values: dict) -> None:
    for key, value in values.items():
        if value is None:
            continue
        if isinstance(value, (dt.datetime, dt.date)):
            value = value.isoformat()
        group.attrs[str(key)] = str(value)


def compression_args(name: str) -> dict:
    if name == "vbz":
        return {"compression": 32020, "compression_opts": (0, 2, 1, 1)}
    if name == "gzip":
        return {"compression": "gzip", "compression_opts": 1, "shuffle": True}
    return {}


def make_background(args: argparse.Namespace) -> None:
    source = Path(args.template).resolve()
    output = Path(args.output).resolve()
    if output.exists() and not args.force:
        raise FileExistsError(f"Refusing to overwrite {output}; use --force")
    output.parent.mkdir(parents=True, exist_ok=True)
    wanted = int(round(args.seconds_per_channel * args.sample_rate))
    minimum = int(round(args.minimum_segment_seconds * args.sample_rate))
    margin = int(round(args.edge_margin_seconds * args.sample_rate))

    with h5py.File(source, "r") as src, h5py.File(output, "w") as dst:
        sample_rate = int(src["/Meta"].attrs["sample_rate"])
        if args.sample_rate != sample_rate:
            raise ValueError(
                f"Template is {sample_rate} Hz, not requested {args.sample_rate} Hz"
            )
        dst.attrs.update(
            format="pod5-bulk-background-v1",
            source_file=source.name,
            sample_rate=sample_rate,
            state_selection="StateData summary_state == pore (5)",
            created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
        )
        bg = dst.create_group("background")
        schema = dst.create_group("schema")
        # Preserve exact compound dtypes without copying the large datasets.
        for source_path in (
            "/IntermediateData/Channel_1/Reads",
            "/MultiplexData/Channel_1/Multiplex",
            "/StateData/Channel_1/States",
            "/Device/AsicCommands",
            "/Device/MetaData",
        ):
            if source_path in src:
                parent, name = source_path.rsplit("/", 1)
                out_parent = schema.require_group(parent.lstrip("/"))
                d = src[source_path]
                out_parent.create_dataset(name, shape=(0,), maxshape=(None,), dtype=d.dtype)

        if "/Device/MetaData" in src:
            count = min(wanted, len(src["/Device/MetaData"]))
            schema.create_dataset(
                "device_metadata_seed",
                data=src["/Device/MetaData"][:count],
                compression="gzip",
                compression_opts=1,
            )

        raw = src["/Raw"]
        channels = sorted(int(k.split("_")[-1]) for k in raw if k.startswith("Channel_"))
        failures = []
        successful = []
        for channel in channels:
            state_path = f"/StateData/Channel_{channel}/States"
            signal_path = f"/Raw/Channel_{channel}/Signal"
            if state_path not in src or signal_path not in src:
                failures.append(channel)
                continue
            states = src[state_path][:]
            signal = src[signal_path]
            candidates = []
            for i in range(len(states) - 1):
                if int(states[i]["summary_state"]) != STATE_PORE:
                    continue
                start = int(states[i]["acquisition_raw_index"]) + margin
                end = int(states[i + 1]["acquisition_raw_index"]) - margin
                if end - start >= minimum:
                    candidates.append((end - start, start, end))
            candidates.sort(reverse=True)
            if not candidates:
                failures.append(channel)
                continue
            group = bg.create_group(f"Channel_{channel}")
            meta = raw[f"Channel_{channel}/Meta"].attrs
            for key in ("digitisation", "offset", "range", "sample_rate", "description"):
                if key in meta:
                    group.attrs[key] = meta[key]
            remaining = wanted
            used = 0
            for _, start, end in candidates:
                if remaining <= 0 or used >= args.max_segments:
                    break
                take = min(remaining, end - start)
                group.create_dataset(
                    f"segment_{used}",
                    data=signal[start : start + take],
                    dtype=np.int16,
                    compression="gzip",
                    compression_opts=1,
                    shuffle=True,
                )
                group[f"segment_{used}"].attrs["source_start_sample"] = start
                remaining -= take
                used += 1
            # A shorter genuine segment is preferable to fabricating samples; it will
            # simply cycle more often in the output.
            if used == 0:
                failures.append(channel)
                del bg[f"Channel_{channel}"]
            else:
                group.attrs["cached_samples"] = wanted - remaining
                successful.append(channel)
                if channel % 32 == 0:
                    log(f"Cached open-pore background through channel {channel}")
        if failures:
            if args.strict_channels:
                raise RuntimeError(
                    "No qualifying pore-state background for channels: "
                    + ",".join(map(str, failures))
                )
            if not successful:
                raise RuntimeError("No qualifying pore-state background was found")
            fallbacks = {}
            for channel in failures:
                donor = min(successful, key=lambda good: (abs(good - channel), good))
                bg[f"Channel_{channel}"] = bg[f"Channel_{donor}"]
                fallbacks[str(channel)] = donor
            dst.attrs["fallback_channels_json"] = json.dumps(fallbacks, sort_keys=True)
            log(
                "Channels without qualifying pore intervals use the nearest donor "
                f"channel: {fallbacks}"
            )
    log(f"Wrote reusable background cache: {output}")


INDEX_SCHEMA = """
CREATE TABLE files (id INTEGER PRIMARY KEY, path TEXT UNIQUE NOT NULL);
CREATE TABLE reads (
 file_id INTEGER NOT NULL, read_id TEXT NOT NULL, channel INTEGER NOT NULL,
 well INTEGER NOT NULL, start INTEGER NOT NULL, samples INTEGER NOT NULL,
 read_number INTEGER NOT NULL, end_reason TEXT NOT NULL, forced INTEGER NOT NULL,
 median_before REAL, offset REAL NOT NULL, scale REAL NOT NULL,
 num_events INTEGER, tracked_scale REAL, tracked_shift REAL,
 predicted_scale REAL, predicted_shift REAL, reads_since_mux INTEGER,
 time_since_mux REAL, PRIMARY KEY(file_id, read_id));
CREATE INDEX reads_channel_start ON reads(channel, start);
CREATE TABLE run_info (acquisition_id TEXT PRIMARY KEY, json TEXT NOT NULL);
"""


def run_info_dict(ri) -> dict:
    return {
        "acquisition_id": ri.acquisition_id,
        "acquisition_start_time": ri.acquisition_start_time.isoformat(),
        "adc_max": ri.adc_max,
        "adc_min": ri.adc_min,
        "context_tags": dict(ri.context_tags),
        "experiment_name": ri.experiment_name,
        "flow_cell_id": ri.flow_cell_id,
        "flow_cell_product_code": ri.flow_cell_product_code,
        "protocol_name": ri.protocol_name,
        "protocol_run_id": ri.protocol_run_id,
        "protocol_start_time": ri.protocol_start_time.isoformat(),
        "sample_id": ri.sample_id,
        "sample_rate": ri.sample_rate,
        "sequencing_kit": ri.sequencing_kit,
        "sequencer_position": ri.sequencer_position,
        "sequencer_position_type": ri.sequencer_position_type,
        "software": ri.software,
        "system_name": ri.system_name,
        "system_type": ri.system_type,
        "tracking_id": dict(ri.tracking_id),
    }


def pair_value(pair, name: str) -> float:
    value = getattr(pair, name, float("nan"))
    return float(value) if value is not None else float("nan")


def build_index(inputs: list[Path], database: Path, exclude_forced: bool) -> dict:
    if database.exists():
        database.unlink()
    con = sqlite3.connect(database)
    con.executescript(INDEX_SCHEMA)
    total = kept = removed = 0
    end_reasons = collections.defaultdict(lambda: {"total": 0, "retained": 0, "excluded": 0})
    source_start = source_end = None
    try:
        for file_id, path in enumerate(inputs):
            con.execute("INSERT INTO files(id,path) VALUES (?,?)", (file_id, str(path)))
            rows = []
            with pod5.Reader(path) as reader:
                for read in reader.reads():
                    total += 1
                    forced = int(read.end_reason.forced)
                    start = int(read.start_sample)
                    end = start + int(read.sample_count)
                    source_start = start if source_start is None else min(source_start, start)
                    source_end = end if source_end is None else max(source_end, end)
                    counts = end_reasons[read.end_reason.name]
                    counts["total"] += 1
                    # Validate run identity even for reads excluded from the signal.
                    ri = read.run_info
                    info = run_info_dict(ri)
                    con.execute(
                        "INSERT OR IGNORE INTO run_info VALUES (?,?)",
                        (ri.acquisition_id, json.dumps(info, sort_keys=True)),
                    )
                    if exclude_forced and forced:
                        removed += 1
                        counts["excluded"] += 1
                        continue
                    counts["retained"] += 1
                    rows.append(
                        (
                            file_id, str(read.read_id), read.pore.channel, read.pore.well,
                            read.start_sample, read.sample_count, read.read_number,
                            read.end_reason.name, forced, float(read.median_before),
                            float(read.calibration.offset), float(read.calibration.scale),
                            int(read.num_minknow_events),
                            pair_value(read.tracked_scaling, "scale"),
                            pair_value(read.tracked_scaling, "shift"),
                            pair_value(read.predicted_scaling, "scale"),
                            pair_value(read.predicted_scaling, "shift"),
                            int(read.num_reads_since_mux_change),
                            float(read.time_since_mux_change),
                        )
                    )
                    kept += 1
                    if len(rows) >= 10000:
                        con.executemany("INSERT INTO reads VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)
                        rows.clear()
                if rows:
                    con.executemany("INSERT INTO reads VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)
            con.commit()
            log(f"Indexed {path.name}")
        runs = con.execute("SELECT acquisition_id,json FROM run_info").fetchall()
        if len(runs) != 1:
            raise ValueError(
                f"Inputs contain {len(runs)} acquisition IDs; one output must represent one run"
            )
        min_start, max_end, max_channel = con.execute(
            "SELECT min(start), max(start+samples), max(channel) FROM reads"
        ).fetchone()
        if min_start is None:
            raise ValueError("No reads remain after filtering")
        bad = []
        for channel, in con.execute("SELECT DISTINCT channel FROM reads"):
            vals = con.execute(
                "SELECT DISTINCT offset,scale FROM reads WHERE channel=?", (channel,)
            ).fetchall()
            if len(vals) != 1:
                bad.append(channel)
        if bad:
            raise ValueError(f"Calibration changes within channels {bad}; unsupported without rescaling")
        info = json.loads(runs[0][1])
        return {
            "total": total, "kept": kept, "removed": removed,
            "end_reasons": dict(end_reasons),
            "source_min_start": source_start, "source_max_end": source_end,
            "min_start": int(min_start), "max_end": int(max_end),
            "max_channel": int(max_channel), "run_info": info,
        }
    finally:
        con.close()


def load_background(cache: h5py.File, channel: int, target_offset: float, target_scale: float):
    key = f"/background/Channel_{channel}"
    if key not in cache:
        raise KeyError(f"Background cache has no channel {channel}")
    group = cache[key]
    source_offset = float(group.attrs["offset"])
    source_scale = float(group.attrs["range"]) / float(group.attrs["digitisation"])
    segments = []
    for name in sorted(group, key=lambda x: int(x.rsplit("_", 1)[1])):
        raw = group[name][:].astype(np.float32)
        # ONT calibration: pA = (ADC + offset) * scale.
        pa = (raw + source_offset) * source_scale
        converted = np.rint(pa / target_scale - target_offset)
        segments.append(np.clip(converted, -32768, 32767).astype(np.int16))
    return segments


def background_chunk(segments: list[np.ndarray], start: int, size: int) -> np.ndarray:
    total = sum(len(x) for x in segments)
    if total == 0:
        raise ValueError("Empty background cache")
    out = np.empty(size, dtype=np.int16)
    position = start % total
    written = 0
    offsets = np.cumsum([0] + [len(x) for x in segments])
    while written < size:
        index = int(np.searchsorted(offsets, position, side="right") - 1)
        within = position - int(offsets[index])
        take = min(size - written, len(segments[index]) - within)
        out[written : written + take] = segments[index][within : within + take]
        written += take
        position = (position + take) % total
    return out


def create_structure(out: h5py.File, cache: h5py.File, info: dict, duration: int) -> None:
    out.attrs["file_version"] = "bulk-playback-reconstruction-v1"
    meta = out.create_group("Meta")
    meta.attrs["duration_samples"] = np.uint64(duration)
    meta.attrs["sample_rate"] = np.uint32(info["sample_rate"])
    user = meta.create_group("User")
    provenance = {
        "generator": "pod5_to_bulk_fast5.py",
        "warning": "Auxiliary tables are reconstructed; Device/AsicCommands is empty.",
    }
    text = json.dumps(provenance, sort_keys=True)
    user.create_dataset("conversion_provenance", data=np.array([text], dtype=h5py.string_dtype()))
    ugk = out.create_group("UniqueGlobalKey")
    set_string_attrs(ugk.create_group("context_tags"), info["context_tags"])
    tracking = dict(info["tracking_id"])
    tracking["is_simulated"] = "1"
    set_string_attrs(ugk.create_group("tracking_id"), tracking)
    out.create_group("Raw")
    out.create_group("IntermediateData")
    out.create_group("MultiplexData")
    out.create_group("StateData")
    device = out.create_group("Device")
    for name in ("AsicCommands", "MetaData"):
        source = f"/schema/Device/{name}"
        if source in cache:
            device.create_dataset(name, shape=(0,), maxshape=(None,), dtype=cache[source].dtype)


def write_device_metadata(out: h5py.File, cache: h5py.File, duration: int, chunk: int) -> None:
    if "/schema/Device/MetaData" not in cache or "/schema/device_metadata_seed" not in cache:
        return
    if "/Device/MetaData" in out:
        del out["/Device/MetaData"]
    seed = cache["/schema/device_metadata_seed"][:]
    dset = out["/Device"].create_dataset(
        "MetaData", shape=(duration,), maxshape=(None,), dtype=seed.dtype,
        chunks=(chunk,), compression="gzip", compression_opts=1, shuffle=True,
    )
    for start in range(0, duration, chunk):
        n = min(chunk, duration - start)
        indexes = (np.arange(n) + start) % len(seed)
        dset[start : start + n] = seed[indexes]


def scalar_or_nan(value):
    return value if value is not None else np.nan


def write_auxiliary_channel(out: h5py.File, cache: h5py.File, channel: int, rows, origin: int, duration: int):
    # Reads: one reconstructed record per POD5 read, rather than MinKNOW's internal segments.
    read_dtype = cache["/schema/IntermediateData/Channel_1/Reads"].dtype
    arr = np.zeros(len(rows), dtype=read_dtype)
    names = set(read_dtype.names or ())
    for i, row in enumerate(rows):
        start = int(row[4]) - origin
        length = max(0, min(duration, start + int(row[5])) - max(0, start))
        values = {
            "read_id": row[1].encode(), "end_reason": END_REASON.get(row[7], 0),
            "read_number": row[6], "read_start": max(0, start), "read_length": length,
            "classification": READ_CLASS_STRAND, "pen_classification": READ_CLASS_STRAND,
            "modal_classification": READ_CLASS_STRAND, "channel_config": row[3],
            "median_before": row[9], "median": row[9],
            "num_reads_since_mux_change": row[17], "time_since_mux_change": row[18],
        }
        for key, value in values.items():
            if key in names:
                arr[key][i] = value
        for key, scale_i, shift_i in (
            ("tracked_scaling", 13, 14), ("predicted_scaling", 15, 16)
        ):
            if key in names and arr.dtype[key].names:
                arr[key][i]["scale"] = row[scale_i]
                arr[key][i]["shift"] = row[shift_i]
    group = out["/IntermediateData"].create_group(f"Channel_{channel}")
    group.create_group("Meta")
    group.create_dataset("Reads", data=arr, maxshape=(None,), chunks=True, compression="gzip", compression_opts=1)

    mux_dtype = cache["/schema/MultiplexData/Channel_1/Multiplex"].dtype
    mux_events = [(0, 0)]
    previous = 0
    for row in rows:
        start, well = max(0, int(row[4]) - origin), int(row[3])
        if 0 <= start < duration and well != previous:
            mux_events.append((start, well)); previous = well
    mux = np.zeros(len(mux_events), dtype=mux_dtype)
    for i, (start, well) in enumerate(mux_events):
        mux["approx_raw_index"][i] = start
        if "well" in mux.dtype.names:
            mux["well"][i] = well
        if "well_id" in mux.dtype.names:
            mux["well_id"][i] = well
    group = out["/MultiplexData"].create_group(f"Channel_{channel}")
    group.create_group("Meta")
    group.create_dataset("Multiplex", data=mux, maxshape=(None,), chunks=True, compression="gzip", compression_opts=1)

    state_dtype = cache["/schema/StateData/Channel_1/States"].dtype
    events = [(0, STATE_PORE)]
    for row in rows:
        start = max(0, int(row[4]) - origin)
        end = min(duration, int(row[4]) - origin + int(row[5]))
        if start < duration:
            events.append((start, STATE_STRAND))
        if end < duration:
            events.append((end, STATE_PORE))
    # At equal indices, strand takes precedence over pore.
    merged = {}
    for index, state in events:
        merged[index] = max(state, merged.get(index, -1))
    events = sorted(merged.items())
    states = np.zeros(len(events), dtype=state_dtype)
    for i, (index, state) in enumerate(events):
        for field in ("acquisition_raw_index", "analysis_raw_index"):
            if field in states.dtype.names:
                states[field][i] = index
        if "trigger_time" in states.dtype.names:
            states["trigger_time"][i] = 0
        states["summary_state"][i] = state
    group = out["/StateData"].create_group(f"Channel_{channel}")
    group.create_group("Meta")
    group.create_dataset("States", data=states, maxshape=(None,), chunks=True, compression="gzip", compression_opts=1)


def positive_seconds(value: str) -> float:
    seconds = float(value)
    if not math.isfinite(seconds) or seconds <= 0:
        raise argparse.ArgumentTypeError("duration must be finite and greater than zero")
    return seconds


def nonnegative_seconds(value: str) -> float:
    seconds = float(value)
    if not math.isfinite(seconds) or seconds < 0:
        raise argparse.ArgumentTypeError("window start must be finite and nonnegative")
    return seconds


def window_predicate(origin: int, duration: int, complete_only: bool) -> tuple[str, tuple]:
    """Half-open intersection, optionally restricted to fully contained reads."""
    end = origin + duration
    condition = "start < ? AND start+samples > ?"
    params = (end, origin)
    if complete_only:
        condition += " AND start >= ? AND start+samples <= ?"
        params += (origin, end)
    return condition, params


def resolve_timeline(args, summary: dict) -> tuple[int, int]:
    """Choose bounds independently of which read signals are retained."""
    rate = int(summary["run_info"]["sample_rate"])
    if rate <= 0:
        raise ValueError("Sample rate must be positive")
    if args.timeline == "source":
        first, last = summary["source_min_start"], summary["source_max_end"]
    else:
        first, last = summary["min_start"], summary["max_end"]
    origin = int(first) if args.time_origin == "rebase" else 0
    origin += int(round(args.window_start_seconds * rate))
    if origin >= int(last):
        raise ValueError("Window start must precede the end of the selected timeline")
    duration = int(last) - origin
    if args.duration_seconds is not None:
        requested = int(round(args.duration_seconds * rate))
        if requested < duration:
            raise ValueError(
                "--duration-seconds cannot shorten the selected timeline; "
                "use --max-duration-seconds for a test window"
            )
        duration = requested
    if args.max_duration_seconds is not None:
        duration = min(duration, int(round(args.max_duration_seconds * rate)))
    if duration <= 0:
        raise ValueError("Output duration must contain at least one sample")
    return origin, duration


def conversion_report(args, summary, inputs, output, cache_path, database, origin, duration, channels):
    """Describe the completed artifact, without changing conversion policy."""
    info = summary["run_info"]
    rate = int(info["sample_rate"])
    hours = max(1, math.ceil(duration / (rate * 3600)))
    product = str(info.get("flow_cell_product_code") or "<PRODUCT_CODE>").upper()
    kit = str(info.get("sequencing_kit") or "<KIT>").upper()
    flags = [
        "--position", "MS00000", "--product-code", product, "--kit", kit,
        "--experiment-duration", str(hours),
        "--fastq", "--bam", "--pod5", "--verbose", "--basecalling",
        "--simulation", f"/tmp/.dorado/{output.name}",
    ]
    with sqlite3.connect(f"{database.as_uri()}?mode=ro", uri=True) as con:
        intersecting, crossing = con.execute(
            "SELECT count(*), coalesce(sum(start < ? OR start+samples > ?), 0) "
            "FROM reads WHERE channel BETWEEN 1 AND ? AND start < ? AND start+samples > ?",
            (origin, origin + duration, channels, origin + duration, origin),
        ).fetchone()
    return {
        "report_version": 1,
        "status": "conversion_completed",
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "output": str(output), "output_bytes": output.stat().st_size,
        "inputs": [{"path": str(p), "bytes": p.stat().st_size} for p in inputs],
        "input_file_count": len(inputs),
        "background_cache": str(cache_path), "index": str(database),
        "run_info": info,
        "reads": {
            "total": summary["total"], "retained_in_index": summary["kept"],
            "excluded_forced": summary["removed"], "end_reasons": summary["end_reasons"],
            "intersecting_output": intersecting,
            "boundary_crossing_reads": crossing,
            "selected_for_overlay": intersecting - crossing if args.complete_reads_only else intersecting,
            "excluded_at_output_boundary": crossing if args.complete_reads_only else 0,
            "clipped_at_output_boundary": 0 if args.complete_reads_only else crossing,
            "retained_outside_output": summary["kept"] - intersecting,
        },
        "timeline": {
            "sample_rate": rate, "channels": channels,
            "source_first_read_start_sample": summary["source_min_start"],
            "source_last_read_end_sample": summary["source_max_end"],
            "retained_first_read_start_sample": summary["min_start"],
            "retained_last_read_end_sample": summary["max_end"],
            "output_origin_sample": origin, "duration_samples": duration,
            "duration_seconds": duration / rate, "duration_hours": duration / rate / 3600,
        },
        "settings": {key: getattr(args, key) for key in (
            "timeline", "time_origin", "duration_seconds", "exclude_forced", "auxiliary", "device_metadata",
            "compression", "chunk_samples", "channels", "max_duration_seconds",
            "window_start_seconds", "complete_reads_only",
        )},
        "recommended_simulation": {
            "flags": flags, "shell_flags": shlex.join(flags),
            "notes": [
                "Append these flags to your installed MinKNOW start_protocol.py command.",
                "Position MS00000 and the container simulation path are editable assumptions.",
                "Experiment duration is in hours, rounded up from the generated signal duration (minimum 1).",
                "This is a protocol time limit, not a guarantee of playback end-of-file behavior.",
                "Kit and product code come from POD5 metadata; verify installed protocol support.",
                "Replace <KIT> or <PRODUCT_CODE> if source metadata is missing.",
                "FASTQ output may be compressed; compare read/base counts rather than file sizes.",
            ],
        },
        "limitations": [
            "Conversion completion is not a signal-integrity or playback validation.",
            "POD5 cannot establish experiment stop time after the last supplied read.",
            "Intersecting read counts do not establish unique recovery of overlapping reads.",
        ],
    }


def convert(args: argparse.Namespace) -> None:
    inputs = expand_inputs(args.inputs)
    output = Path(args.output).resolve()
    cache_path = Path(args.background_cache).resolve()
    tmp_dir = Path(args.tmp_dir).resolve()
    tmp_dir.mkdir(parents=True, exist_ok=True)
    database = tmp_dir / "pod5_bulk_index.sqlite"
    report_path = output.with_suffix(output.suffix + ".report.json")
    for path in (output, report_path):
        if path.exists() and not args.force:
            raise FileExistsError(f"Refusing to overwrite {path}; use --force")
    # Do not leave an old success report beside an interrupted forced rebuild.
    report_path.unlink(missing_ok=True)

    summary = build_index(inputs, database, args.exclude_forced)
    info = summary["run_info"]
    origin, duration = resolve_timeline(args, summary)
    channels = args.channels or max(512, summary["max_channel"])
    log(json.dumps({**summary, "run_info": info["acquisition_id"], "origin": origin, "duration": duration}, indent=2))

    con = sqlite3.connect(database)
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        with h5py.File(cache_path, "r") as cache, h5py.File(output, "w") as out, ExitStack() as stack:
            if int(cache.attrs["sample_rate"]) != int(info["sample_rate"]):
                raise ValueError("POD5 and background cache sample rates differ")
            create_structure(out, cache, info, duration)
            readers = {
                file_id: stack.enter_context(pod5.Reader(path))
                for file_id, path in con.execute("SELECT id,path FROM files")
            }
            comp = compression_args(args.compression)
            chunk_size = min(args.chunk_samples, max(1, duration))
            condition, window_params = window_predicate(origin, duration, args.complete_reads_only)
            for channel in range(1, channels + 1):
                rows = con.execute(
                    "SELECT file_id,read_id,channel,well,start,samples,read_number,end_reason,forced,"
                    "median_before,offset,scale,num_events,tracked_scale,tracked_shift,predicted_scale,"
                    "predicted_shift,reads_since_mux,time_since_mux FROM reads "
                    f"WHERE channel=? AND {condition} ORDER BY start",
                    (channel, *window_params),
                ).fetchall()
                if rows:
                    target_offset, target_scale = float(rows[0][10]), float(rows[0][11])
                else:
                    bg_group = cache[f"/background/Channel_{channel}"]
                    target_offset = float(bg_group.attrs["offset"])
                    target_scale = float(bg_group.attrs["range"]) / float(bg_group.attrs["digitisation"])
                segments = load_background(cache, channel, target_offset, target_scale)

                signals = {}
                by_file = collections.defaultdict(list)
                for row in rows:
                    by_file[row[0]].append(row[1])
                for file_id, ids in by_file.items():
                    for read in readers[file_id].reads(selection=ids):
                        signals[str(read.read_id)] = np.asarray(read.signal, dtype=np.int16)
                if len(signals) != len(rows):
                    raise RuntimeError(f"Failed to retrieve all signals for channel {channel}")

                raw_group = out["/Raw"].create_group(f"Channel_{channel}")
                meta = raw_group.create_group("Meta")
                digitisation = int(info["adc_max"]) - int(info["adc_min"]) + 1
                meta.attrs.update(
                    description="Grouper", digitisation=float(digitisation),
                    offset=float(target_offset), range=float(digitisation * target_scale),
                    sample_rate=float(info["sample_rate"]),
                )
                dset = raw_group.create_dataset(
                    "Signal", shape=(duration,), maxshape=(None,), dtype=np.int16,
                    chunks=(chunk_size,), **comp,
                )
                row_index = 0
                for start in range(0, duration, chunk_size):
                    end = min(duration, start + chunk_size)
                    data = background_chunk(segments, start, end - start)
                    while row_index < len(rows) and rows[row_index][4] - origin + rows[row_index][5] <= start:
                        row_index += 1
                    j = row_index
                    while j < len(rows) and rows[j][4] - origin < end:
                        row = rows[j]
                        read_start = int(row[4]) - origin
                        signal = signals[row[1]]
                        left = max(start, read_start)
                        right = min(end, read_start + len(signal))
                        if right > left:
                            data[left - start : right - start] = signal[left - read_start : right - read_start]
                        j += 1
                    dset[start:end] = data
                del signals
                if args.auxiliary == "reconstructed":
                    write_auxiliary_channel(out, cache, channel, rows, origin, duration)
                if channel % 16 == 0 or channel == channels:
                    log(f"Wrote channel {channel}/{channels}")
            if args.auxiliary == "reconstructed" and args.device_metadata:
                write_device_metadata(out, cache, duration, chunk_size)
            out.flush()
    finally:
        con.close()
    report = conversion_report(
        args, summary, inputs, output, cache_path, database, origin, duration, channels,
    )
    pending = report_path.with_suffix(report_path.suffix + ".tmp")
    pending.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    pending.replace(report_path)
    log(f"Conversion report: {report_path}")
    if report["reads"]["selected_for_overlay"] == 0:
        log("WARNING: No retained reads selected for this window/channel range; output is background-only.")
    log("Recommended simulation flags (edit position/container path if needed):")
    log(report["recommended_simulation"]["shell_flags"])
    log(f"Wrote bulk playback FAST5: {output}")
    log(f"Index retained for audit: {database}")


def validate(args: argparse.Namespace) -> None:
    path = Path(args.fast5)
    with h5py.File(path, "r") as f:
        required = ["/Meta", "/Raw", "/UniqueGlobalKey/context_tags", "/UniqueGlobalKey/tracking_id"]
        missing = [x for x in required if x not in f]
        if missing:
            raise ValueError(f"Missing: {missing}")
        duration = int(f["/Meta"].attrs["duration_samples"])
        rate = int(f["/Meta"].attrs["sample_rate"])
        channels = sorted(k for k in f["/Raw"] if k.startswith("Channel_"))
        lengths = {len(f[f"/Raw/{ch}/Signal"]) for ch in channels}
        if lengths != {duration}:
            raise ValueError(f"Channel lengths do not all match duration: {lengths}")
        print(json.dumps({
            "file": str(path), "channels": len(channels), "duration_samples": duration,
            "sample_rate": rate, "duration_seconds": duration / rate,
            "auxiliary_groups": [x for x in ("IntermediateData", "MultiplexData", "StateData", "Device") if x in f],
        }, indent=2))


def positive_integer(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("value must be greater than zero")
    return number


def audit(args: argparse.Namespace) -> None:
    from signal_integrity_audit import audit as run_audit
    run_audit(args, sys.modules[__name__])


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    b = sub.add_parser("make-background", help="Extract reusable open-pore signal from a bulk FAST5")
    b.add_argument("template")
    b.add_argument("output")
    b.add_argument("--seconds-per-channel", type=float, default=10.0)
    b.add_argument("--minimum-segment-seconds", type=float, default=1.0)
    b.add_argument("--edge-margin-seconds", type=float, default=0.1)
    b.add_argument("--max-segments", type=int, default=8)
    b.add_argument("--sample-rate", type=int, default=5000)
    b.add_argument(
        "--strict-channels", action="store_true",
        help="fail instead of using a nearest-channel pore segment when a channel has none",
    )
    b.add_argument("--force", action="store_true")
    b.set_defaults(func=make_background)

    c = sub.add_parser("convert", help="Convert one run's POD5 files to bulk playback FAST5")
    c.add_argument("inputs", nargs="+")
    c.add_argument("--output", required=True)
    c.add_argument("--background-cache", required=True)
    c.add_argument("--tmp-dir", default="tmp")
    c.add_argument(
        "--timeline", choices=("source", "retained"), default="source",
        help="Timing bounds from all supplied reads (default) or only retained reads",
    )
    c.add_argument(
        "--time-origin", choices=("rebase", "absolute"), default="absolute",
        help="Start at acquisition zero (default), or the first read in the selected timeline",
    )
    c.add_argument(
        "--window-start-seconds", type=nonnegative_seconds, default=0.0,
        help="Test-window offset from the selected time origin (default: 0 seconds)",
    )
    c.add_argument(
        "--complete-reads-only", action="store_true",
        help="Skip reads crossing either output boundary instead of clipping them; independent of forced-read filtering",
    )
    c.add_argument(
        "--exclude-forced", action=argparse.BooleanOptionalAction, default=False,
        help="Opt in to excluding forced-ended reads (e.g. unblocks/mux changes); included by default",
    )
    c.add_argument("--auxiliary", choices=("none", "reconstructed"), default="reconstructed")
    c.add_argument("--device-metadata", action=argparse.BooleanOptionalAction, default=True)
    c.add_argument("--compression", choices=("vbz", "gzip", "none"), default="vbz")
    c.add_argument("--chunk-samples", type=int, default=180480)
    c.add_argument("--channels", type=int, help="Testing only; defaults to at least 512")
    duration_options = c.add_mutually_exclusive_group()
    duration_options.add_argument(
        "--duration-seconds", type=positive_seconds,
        help="Explicit output length from the window start; may extend but not shorten the remaining timeline",
    )
    duration_options.add_argument(
        "--max-duration-seconds", type=positive_seconds,
        help="Testing only: cap output length from the window start; clips reads unless --complete-reads-only is set",
    )
    c.add_argument("--force", action="store_true")
    c.set_defaults(func=convert)

    v = sub.add_parser("validate", help="Perform structural checks without scanning signal")
    v.add_argument("fast5")
    v.set_defaults(func=validate)

    a = sub.add_parser("audit", help="Compare retained POD5 ADC samples and calibration against a completed bulk FAST5")
    a.add_argument("fast5")
    a.add_argument("--report", help="Conversion report (default: <fast5>.report.json)")
    a.add_argument("--mode", choices=("sampled", "full"), default="sampled")
    a.add_argument("--reads", type=positive_integer, default=1000, help="Sampled-mode read budget (default: 1000); ignored in full mode")
    a.add_argument("--seed", type=int, default=42, help="Reproducible sampled selection seed (default: 42)")
    a.add_argument("--chunk-samples", type=positive_integer, default=180480, help="FAST5 comparison slice size; source memory is bounded by one read")
    a.add_argument("--tmp-dir", help="Existing parent directory for a private temporary metadata index (default: system temp)")
    a.add_argument("--force", action="store_true", help="Replace an existing <fast5>.audit.json, never the FAST5 or conversion report")
    a.set_defaults(func=audit)
    return p


def main() -> None:
    args = parser().parse_args()
    try:
        args.func(args)
    except Exception as error:
        log(f"ERROR: {error}")
        if os.environ.get("POD5_BULK_DEBUG"):
            raise
        raise SystemExit(1)


if __name__ == "__main__":
    main()
