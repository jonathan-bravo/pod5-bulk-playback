"""Read-only post-conversion audit of retained POD5 signal and calibration.

The source metadata index is rebuilt in a private temporary directory: the
converter's reusable tmp index is deliberately not trusted. Synthetic background
and reconstructed auxiliary state history are outside this audit's scope.
"""
from __future__ import annotations

import collections
import datetime as dt
import hashlib
import heapq
import json
import math
import os
from pathlib import Path
import sqlite3
import tempfile

import h5py
import numpy as np
import pod5


DETAIL_LIMIT = 100


def prepare_index(con, converter, settings, origin, duration, channels):
    condition, params = converter.window_predicate(origin, duration, settings.complete_reads_only)
    con.execute(
        "CREATE TABLE audit_reads AS SELECT rowid AS audit_id, * FROM reads "
        f"WHERE channel BETWEEN 1 AND ? AND {condition}",
        (channels, *params),
    )
    con.execute("CREATE INDEX audit_order ON audit_reads(channel,start,audit_id)")
    con.execute("CREATE INDEX audit_file ON audit_reads(file_id,audit_id)")
    con.execute("CREATE UNIQUE INDEX audit_identity ON audit_reads(audit_id)")
    con.execute("CREATE TABLE overlap_reads(audit_id INTEGER PRIMARY KEY)")
    examples = []
    events = 0
    channel = None
    furthest_end = -1
    previous = None
    # A running maximum is enough to identify every read involved in an overlap,
    # without enumerating a potentially quadratic number of overlapping pairs.
    for row in con.execute(
        "SELECT audit_id,file_id,read_id,channel,start,samples FROM audit_reads "
        "ORDER BY channel,start,audit_id"
    ):
        left, right = max(origin, row[4]), min(origin + duration, row[4] + row[5])
        if row[3] != channel:
            channel, furthest_end, previous = row[3], -1, None
        if left < furthest_end:
            events += 1
            con.executemany("INSERT OR IGNORE INTO overlap_reads VALUES (?)", [(row[0],), (previous[0],)])
            if len(examples) < DETAIL_LIMIT:
                examples.append({
                    "channel": channel,
                    "read_ids": [previous[2], row[2]],
                    "file_ids": [previous[1], row[1]],
                    "source_start_sample": left, "source_end_sample": min(right, furthest_end),
                })
        if right > furthest_end:
            furthest_end, previous = right, row
    con.commit()
    return {
        "scope": "all selected read intervals, including in sampled mode",
        "reads_involved": con.execute("SELECT count(*) FROM overlap_reads").fetchone()[0],
        "overlap_events": events,
        "examples": examples,
        "note": "Events are not a count of every overlapping pair. Overlaps are conflicts even if shared samples agree.",
    }


def choose_reads(con, mode, budget, seed, origin, duration):
    total = con.execute("SELECT count(*) FROM audit_reads").fetchone()[0]
    con.execute("CREATE TABLE audit_selection(audit_id INTEGER PRIMARY KEY)")
    if mode == "full" or budget >= total:
        con.execute("INSERT INTO audit_selection SELECT audit_id FROM audit_reads")
    else:
        representatives = {}
        reservoir = []
        for ident, file_id, read_id, channel, start, samples in con.execute(
            "SELECT audit_id,file_id,read_id,channel,start,samples FROM audit_reads"
        ):
            score = int.from_bytes(hashlib.blake2b(
                f"{seed}:{file_id}:{read_id}".encode(), digest_size=8,
            ).digest(), "big")
            candidate = (score, ident)
            bucket = min(23, max(0, (start - origin) * 24 // duration))
            groups = [("file", file_id), ("channel", channel), ("time", bucket)]
            if start < origin or start + samples > origin + duration:
                groups.append(("boundary", 0))
            for group in groups:
                if group not in representatives or candidate < representatives[group]:
                    representatives[group] = candidate
            entry = (-score, -ident)
            if len(reservoir) < budget:
                heapq.heappush(reservoir, entry)
            elif entry > reservoir[0]:
                heapq.heapreplace(reservoir, entry)
        # Prefer representatives; a small budget may not cover every group.
        chosen = sorted(set(representatives.values()))[:budget]
        ids = {ident for _, ident in chosen}
        for score, ident in sorted((-score, -ident) for score, ident in reservoir):
            if len(ids) >= budget:
                break
            ids.add(ident)
        con.executemany("INSERT INTO audit_selection VALUES (?)", [(ident,) for ident in sorted(ids)])
    con.commit()
    selected, files, channels = con.execute(
        "SELECT count(*),count(DISTINCT file_id),count(DISTINCT channel) "
        "FROM audit_reads JOIN audit_selection USING(audit_id)"
    ).fetchone()
    digest = hashlib.sha256()
    for file_id, read_id in con.execute(
        "SELECT file_id,read_id FROM audit_reads JOIN audit_selection USING(audit_id) ORDER BY file_id,read_id"
    ):
        digest.update(f"{file_id}:{read_id}\n".encode())
    return {"eligible_reads": total, "selected_reads": selected,
            "selected_files": files, "selected_channels": channels,
            "selection_sha256": digest.hexdigest(),
            "sampling_note": "Representatives across input files, channels and 24 time bins are prioritized; small budgets may not cover every group.",
            "all_eligible_reads_selected": selected == total}


def compare_signals(con, bulk, origin, duration, info, chunk_samples, log=None):
    counts = collections.Counter()
    details = []
    checked_channels = set()
    source_paths = dict(con.execute("SELECT id,path FROM files"))

    def problem(kind, row, **extra):
        counts[kind] += 1
        if len(details) < DETAIL_LIMIT:
            details.append({"kind": kind, "file_id": row[1], "source_file": source_paths[row[1]],
                            "read_id": row[2], "channel": row[3], **extra})

    for file_id, path in con.execute(
        "SELECT DISTINCT files.id,files.path FROM files "
        "JOIN audit_reads ON files.id=audit_reads.file_id JOIN audit_selection USING(audit_id) "
        "ORDER BY files.id"
    ):
        if log:
            log(f"Audit: comparing signal from {path}")
        with pod5.Reader(path) as reader:
            cursor = con.execute(
                "SELECT audit_id,file_id,read_id,channel,start,samples,offset,scale "
                "FROM audit_reads JOIN audit_selection USING(audit_id) "
                "WHERE file_id=? ORDER BY audit_id", (file_id,),
            )
            while rows := cursor.fetchmany(128):
                pending = {row[2]: row for row in rows}
                # Only this small ID batch and one decompressed source read are held.
                for read in reader.reads(selection=list(pending)):
                    read_id = str(read.read_id)
                    if read_id not in pending:
                        raise ValueError(f"Unexpected or duplicate read returned from {path}: {read_id}")
                    row = pending.pop(read_id)
                    _, _, _, channel, start, samples, offset, scale = row
                    if (read.start_sample, read.sample_count, read.pore.channel) != (start, samples, channel):
                        raise ValueError(f"Source metadata changed during audit: {path}: {read_id}")
                    raw = bulk[f"Raw/Channel_{channel}"]
                    if channel not in checked_channels:
                        checked_channels.add(channel)
                        attrs = raw["Meta"].attrs
                        expected_digitisation = int(info["adc_max"]) - int(info["adc_min"]) + 1
                        expected = {"offset": offset, "digitisation": expected_digitisation,
                                    "range": expected_digitisation * scale,
                                    "sample_rate": int(info["sample_rate"])}
                        for key, value in expected.items():
                            actual = float(attrs.get(key, float("nan")))
                            if not math.isclose(actual, value, rel_tol=1e-9, abs_tol=1e-9):
                                problem("calibration_mismatches", row, field=key, expected=value,
                                        actual=actual if math.isfinite(actual) else None)
                    signal = np.asarray(read.signal)
                    if len(signal) != samples:
                        raise ValueError(f"Source signal length mismatch: {path}: {read_id}")
                    left, right = max(start, origin), min(start + samples, origin + duration)
                    counts["checked_reads"] += 1
                    counts["checked_samples"] += right - left
                    mismatches = 0
                    first = None
                    for begin in range(left, right, chunk_samples):
                        stop = min(right, begin + chunk_samples)
                        expected = signal[begin-start:stop-start]
                        actual = raw["Signal"][begin-origin:stop-origin]
                        different = np.flatnonzero(expected != actual)
                        mismatches += len(different)
                        if first is None and len(different):
                            index = int(different[0])
                            first = {"source_sample": begin + index, "bulk_sample": begin + index - origin,
                                     "expected_adc": int(expected[index]), "actual_adc": int(actual[index])}
                    if mismatches:
                        counts["mismatched_samples"] += mismatches
                        problem("mismatched_reads", row, mismatched_samples=mismatches, first_difference=first)
                    else:
                        counts["matching_reads"] += 1
                    del signal
                for row in pending.values():
                    problem("missing_reads", row)
    keys = ("checked_reads", "checked_samples", "matching_reads", "mismatched_reads",
            "mismatched_samples", "missing_reads", "calibration_mismatches")
    return {**{key: counts[key] for key in keys}, "calibration_channels_checked": len(checked_channels),
            "details": details, "detail_limit": DETAIL_LIMIT,
            "sample_count_note": "Samples are counted per read; overlapping samples can be counted more than once."}


def audit(args, converter):
    bulk_path = Path(args.fast5).resolve()
    conversion_path = Path(args.report).resolve() if args.report else Path(str(bulk_path) + ".report.json")
    output = Path(str(bulk_path) + ".audit.json")
    if output.is_symlink():
        raise ValueError("Refusing to write an audit report through a symlink")
    if output in (bulk_path, conversion_path.resolve()):
        raise ValueError("Audit output collides with the FAST5 or conversion report")
    if output.exists() and not args.force:
        raise FileExistsError(f"Refusing to overwrite {output}; use --force")
    write_report = True
    result = {
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "audit_version": 1, "status": "error", "fast5": str(bulk_path),
        "conversion_report": str(conversion_path), "mode": args.mode,
        "requested_reads": args.reads, "seed": args.seed,
        "scope": "Retained raw ADC signal and selected-channel calibration; not background or auxiliary history.",
        "limitations": [
            "Sampled success does not prove integrity of unchecked reads.",
            "Source identity checks use paths, sizes and metadata, not original conversion-time cryptographic hashes.",
            "POD5 sources must remain unchanged while auditing.",
        ],
    }
    try:
        report = json.loads(conversion_path.read_text())
        if report.get("report_version") != 1 or report.get("status") != "conversion_completed":
            raise ValueError("Expected a completed version-1 conversion report")
        inputs = [Path(entry["path"]).resolve() for entry in report["inputs"]]
        if output in inputs:
            write_report = False
            raise ValueError("Audit output collides with an input")
        for entry, path in zip(report["inputs"], inputs):
            if path.stat().st_size != entry["bytes"]:
                raise ValueError(f"Source size differs from conversion report: {path}")
        # Recover exact saved settings, not today's CLI defaults. Early reports
        # predate timeline selection and used retained bounds.
        settings = converter.parser().parse_args([
            "convert", "unused", "--output", "unused", "--background-cache", "unused",
        ])
        saved = report["settings"]
        for key, value in saved.items():
            if not hasattr(settings, key):
                raise ValueError(f"Unknown conversion setting: {key}")
            setattr(settings, key, value)
        settings.timeline = saved.get("timeline", "retained")
        settings.complete_reads_only = saved.get("complete_reads_only", False)
        settings.window_start_seconds = saved.get("window_start_seconds", 0.0)
        if "exclude_forced" not in saved or "time_origin" not in saved:
            raise ValueError("Conversion report lacks required filtering/origin settings")
        with tempfile.TemporaryDirectory(prefix="pod5-signal-audit-", dir=args.tmp_dir) as tmp:
            database = Path(tmp) / "index.sqlite"
            converter.log("Audit: rebuilding source metadata index in a private temporary directory")
            summary = converter.build_index(inputs, database, settings.exclude_forced)
            if summary["run_info"] != report["run_info"]:
                raise ValueError("Source run metadata differs from conversion report")
            for key, report_key in (("total", "total"), ("kept", "retained_in_index"), ("removed", "excluded_forced")):
                if summary[key] != report["reads"][report_key]:
                    raise ValueError(f"Source read count differs from conversion report: {key}")
            if summary["end_reasons"] != report["reads"]["end_reasons"]:
                raise ValueError("Source end-reason counts differ from conversion report")
            origin, duration = converter.resolve_timeline(settings, summary)
            timeline = report["timeline"]
            channels = settings.channels or max(512, summary["max_channel"])
            expected_timeline = {
                "output_origin_sample": origin, "duration_samples": duration, "channels": channels,
                "sample_rate": summary["run_info"]["sample_rate"],
                "source_first_read_start_sample": summary["source_min_start"],
                "source_last_read_end_sample": summary["source_max_end"],
                "retained_first_read_start_sample": summary["min_start"],
                "retained_last_read_end_sample": summary["max_end"],
            }
            for key, value in expected_timeline.items():
                if timeline[key] != value:
                    raise ValueError(f"Reconstructed timeline differs from conversion report: {key}")
            result["timeline"] = expected_timeline
            result["source_files"] = [str(path) for path in inputs]
            with sqlite3.connect(database) as con, h5py.File(bulk_path, "r") as bulk:
                if int(bulk["Meta"].attrs["duration_samples"]) != duration or int(bulk["Meta"].attrs["sample_rate"]) != timeline["sample_rate"]:
                    raise ValueError("FAST5 duration/sample rate differs from conversion report")
                expected_channels = {f"Channel_{channel}" for channel in range(1, channels + 1)}
                if set(bulk["Raw"]) != expected_channels:
                    raise ValueError("FAST5 channels differ from conversion report")
                for name in expected_channels:
                    signal = bulk[f"Raw/{name}/Signal"]
                    if signal.shape != (duration,) or signal.dtype.kind != "i" or signal.dtype.itemsize != 2:
                        raise ValueError(f"Incorrect signal shape or dtype: {name}")
                result["overlaps"] = prepare_index(con, converter, settings, origin, duration, channels)
                result["selection"] = choose_reads(con, args.mode, args.reads, args.seed, origin, duration)
                eligible = result["selection"]["eligible_reads"]
                result["overlaps"]["reads_involved_fraction"] = (
                    result["overlaps"]["reads_involved"] / eligible if eligible else 0.0
                )
                converter.log(f"Audit: checking {result['selection']['selected_reads']} read signals")
                result["comparison"] = compare_signals(
                    con, bulk, origin, duration, summary["run_info"], args.chunk_samples, converter.log,
                )
                comparisons = result["comparison"]
                if any(comparisons[key] for key in ("mismatched_reads", "missing_reads", "calibration_mismatches")):
                    result["status"] = "failed"
                elif result["overlaps"]["reads_involved"]:
                    result["status"] = "overlap_conflict"
                elif comparisons["checked_reads"] == 0:
                    result["status"] = "inconclusive"
                else:
                    result["status"] = "passed" if result["selection"]["all_eligible_reads_selected"] else "passed_sampled"
    except Exception as error:
        result["status"] = "error"
        result["error"] = str(error)
    if not write_report:
        raise ValueError(result["error"])
    # This sidecar is the only persistent file written by the audit.
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=output.parent,
                                     prefix=output.name + ".", suffix=".tmp", delete=False) as handle:
        pending = Path(handle.name)
        try:
            handle.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
            handle.flush()
            handle.close()
            os.replace(pending, output)
        finally:
            pending.unlink(missing_ok=True)
    if "comparison" in result:
        counts = result["comparison"]
        converter.log(
            f"Audit checked {counts['checked_reads']} reads / {counts['checked_samples']} samples; "
            f"mismatched reads: {counts['mismatched_reads']}; "
            f"reads involved in overlaps: {result['overlaps']['reads_involved']}"
        )
    if "error" in result:
        converter.log(f"Audit error: {result['error']}")
    converter.log(f"Signal audit: {result['status']}; report: {output}")
    if result["status"] in ("passed", "passed_sampled"):
        return
    raise SystemExit(2 if result["status"] in ("inconclusive", "error") else 1)
