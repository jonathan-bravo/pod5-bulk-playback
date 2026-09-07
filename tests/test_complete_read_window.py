"""Deterministic boundary fixtures; synthetic ADC values are not biological models."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import uuid

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("window_converter", ROOT / "pod5_to_bulk_fast5.py")
converter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(converter)
CACHE = ROOT / "background_r10_4_1_5khz.h5"


def args_for(*options):
    return converter.parser().parse_args([
        "convert", "input.pod5", "--output", "out.fast5", "--background-cache", str(CACHE),
        *options,
    ])


@pytest.mark.parametrize("options,expected", [
    ([], False),
    (["--exclude-forced"], True),
    (["--no-exclude-forced"], False),
    (["--max-duration-seconds", "30", "--complete-reads-only"], False),
    (["--max-duration-seconds", "30", "--complete-reads-only", "--exclude-forced"], True),
])
def test_forced_reads_included_by_default_independently_of_window(options, expected):
    assert args_for(*options).exclude_forced is expected


def test_window_offset_uses_selected_origin():
    summary = {
        "run_info": {"sample_rate": 10}, "source_min_start": 100, "source_max_end": 1000,
        "min_start": 200, "max_end": 800,
    }
    for options, expected in [
        ([], 50),
        (["--time-origin", "rebase"], 150),
        (["--timeline", "retained", "--time-origin", "rebase"], 250),
    ]:
        args = args_for(*options, "--window-start-seconds", "5", "--max-duration-seconds", "30")
        assert converter.resolve_timeline(args, summary) == (expected, 300)
    # A cap never pads beyond the remaining source timeline.
    assert converter.resolve_timeline(
        args_for("--window-start-seconds", "95", "--max-duration-seconds", "30"), summary,
    ) == (950, 50)
    for offset in (100, 101):
        with pytest.raises(ValueError, match="Window start must precede"):
            converter.resolve_timeline(args_for("--window-start-seconds", str(offset)), summary)


@pytest.mark.parametrize("value", ["-1", "nan", "inf", "-inf"])
def test_invalid_window_start(value):
    with pytest.raises(SystemExit):
        args_for(f"--window-start-seconds={value}")


@pytest.mark.parametrize("complete,exclude_forced,window_length", [
    (False, True, 300), (True, True, 300), (True, False, 300),
    (True, True, 30),  # all-background window: no read fits entirely
])
def test_window_signal_auxiliary_and_report(tmp_path, monkeypatch, capsys, complete, exclude_forced, window_length):
    source = tmp_path / "input.pod5"
    source.touch()
    output = tmp_path / "window.fast5"
    origin = 200
    end = origin + window_length
    info = {
        "acquisition_id": "test", "sample_rate": 5000, "adc_min": -32768, "adc_max": 32767,
        "context_tags": {}, "tracking_id": {},
        "sequencing_kit": "sqk-rbk114-24", "flow_cell_product_code": "FLO-MIN114",
    }
    with converter.h5py.File(CACHE, "r") as cache:
        meta = cache["background/Channel_1"].attrs
        offset = float(meta["offset"])
        scale = float(meta["range"]) / float(meta["digitisation"])

    # Includes exact boundary contacts, both crossing directions, a spanning
    # read, and a fully contained forced-ended read on otherwise separate intervals.
    layout = [
        (1, 100, 100, False), (1, 150, 100, False), (1, 260, 50, False),
        (1, 330, 20, True), (1, 450, 100, False), (1, 500, 100, False),
        (2, 200, 50, False), (2, 450, 50, False),
        (3, 100, 500, False), (4, 180, 50, False),
    ]
    reads = []
    for number, (channel, start, length, forced) in enumerate(layout, 1):
        reads.append(SimpleNamespace(
            start_sample=start, sample_count=length, read_id=uuid.UUID(int=number),
            signal=np.arange(length, dtype=np.int16) + np.int16(number * 1000),
            end_reason=SimpleNamespace(forced=forced, name="unblock_mux_change" if forced else "signal_positive"),
            run_info=SimpleNamespace(acquisition_id="test"),
            pore=SimpleNamespace(channel=channel, well=1), read_number=number,
            median_before=0, calibration=SimpleNamespace(offset=offset, scale=scale),
            num_minknow_events=0, tracked_scaling=None, predicted_scaling=None,
            num_reads_since_mux_change=0, time_since_mux_change=0,
        ))
    loaded = []

    class Reader:
        def __init__(self, path):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def reads(self, selection=None):
            if selection is None:
                return iter(reads)
            loaded.extend(selection)
            return (read for read in reads if str(read.read_id) in selection)

    monkeypatch.setattr(converter.pod5, "Reader", Reader)
    monkeypatch.setattr(converter, "run_info_dict", lambda ri: info)
    args = args_for(
        "--output", str(output), "--tmp-dir", str(tmp_path / "index"),
        "--channels", "4", "--compression", "gzip", "--chunk-samples", "31",
        "--no-device-metadata", "--window-start-seconds", str(origin / 5000),
        "--max-duration-seconds", str(window_length / 5000),
        *(["--complete-reads-only"] if complete else []),
        *(["--exclude-forced"] if exclude_forced else []),
    )
    args.inputs = [str(source)]
    converter.convert(args)

    retained = [r for r in reads if not (exclude_forced and r.end_reason.forced)]
    intersects = [r for r in retained if r.start_sample < end and r.start_sample + r.sample_count > origin]
    crossing = [r for r in intersects if r.start_sample < origin or r.start_sample + r.sample_count > end]
    crossing_ids = {r.read_id for r in crossing}
    selected = [r for r in intersects if not complete or r.read_id not in crossing_ids]
    assert ("output is background-only" in capsys.readouterr().err) == (not selected)
    assert set(loaded) == {str(r.read_id) for r in selected}
    report = json.loads(Path(str(output) + ".report.json").read_text())
    assert report["reads"]["selected_for_overlay"] == len(selected)
    assert report["reads"]["excluded_at_output_boundary"] == (len(crossing) if complete else 0)
    assert report["reads"]["clipped_at_output_boundary"] == (0 if complete else len(crossing))
    assert report["reads"]["excluded_forced"] == int(exclude_forced)
    assert report["timeline"]["output_origin_sample"] == origin
    assert report["timeline"]["duration_samples"] == window_length

    with converter.h5py.File(output, "r") as out, converter.h5py.File(CACHE, "r") as cache:
        for channel in range(1, 5):
            raw = out[f"Raw/Channel_{channel}"]
            meta = raw["Meta"].attrs
            segments = converter.load_background(
                cache, channel, float(meta["offset"]), float(meta["range"]) / float(meta["digitisation"]),
            )
            expected = converter.background_chunk(segments, 0, window_length)
            channel_reads = sorted([r for r in selected if r.pore.channel == channel], key=lambda r: r.start_sample)
            for read in channel_reads:
                left, right = max(origin, read.start_sample), min(end, read.start_sample + read.sample_count)
                expected[left-origin:right-origin] = read.signal[left-read.start_sample:right-read.start_sample]
            assert np.array_equal(raw["Signal"][:], expected)
            aux = out[f"IntermediateData/Channel_{channel}/Reads"][:]
            assert len(aux) == len(channel_reads)
            for record, read in zip(aux, channel_reads):
                assert record["read_id"].decode() == str(read.read_id)
                assert record["read_start"] == max(0, read.start_sample - origin)
                assert record["read_length"] == min(end, read.start_sample + read.sample_count) - max(origin, read.start_sample)
            # A read crossing the left boundary on channel 4 ends at sample 30,
            # not its original full length of 50 samples.
            if channel == 4 and not complete:
                states = out["StateData/Channel_4/States"][:]
                assert states[-1]["acquisition_raw_index"] == 30
