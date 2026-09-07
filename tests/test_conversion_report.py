import importlib.util
import json
import shlex
import sqlite3
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("converter", ROOT / "pod5_to_bulk_fast5.py")
converter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(converter)


@pytest.mark.parametrize("seconds,hours", [(30, 1), (3600, 1), (3601, 2), (75707.253, 22)])
def test_report_flags_and_output_counts(tmp_path, seconds, hours):
    source = tmp_path / "source.pod5"
    source.write_bytes(b"input")
    output = tmp_path / "playback with spaces.fast5"
    output.write_bytes(b"output")
    database = tmp_path / "index.sqlite"
    origin = 100
    duration = round(seconds * 5000)
    with sqlite3.connect(database) as con:
        con.execute("CREATE TABLE reads(channel INTEGER, start INTEGER, samples INTEGER)")
        con.executemany("INSERT INTO reads VALUES (?,?,?)", [
            (1, origin, 100),                   # included, complete
            (1, origin + duration - 10, 20),    # clipped
            (2, origin, 100),                   # omitted channel
            (1, origin + duration, 100),        # outside window
        ])
    args = converter.parser().parse_args([
        "convert", str(source), "--output", str(output),
        "--background-cache", "cache.h5", "--channels", "1",
        "--max-duration-seconds", str(seconds), "--exclude-forced",
    ])
    summary = {
        "run_info": {"sample_rate": 5000, "sequencing_kit": "sqk-rbk114-24",
                     "flow_cell_product_code": "FLO-MIN114"},
        "total": 6, "kept": 4, "removed": 2,
        "end_reasons": {"unblock": {"total": 2, "retained": 0, "excluded": 2}},
        "source_min_start": 0, "source_max_end": origin + duration + 100,
        "min_start": origin, "max_end": origin + duration + 100,
    }
    report = converter.conversion_report(
        args, summary, [source], output, tmp_path / "cache.h5", database,
        origin, duration, 1,
    )
    report = json.loads(json.dumps(report))
    flags = report["recommended_simulation"]["flags"]
    assert shlex.split(report["recommended_simulation"]["shell_flags"]) == flags
    assert flags[flags.index("--kit") + 1] == "SQK-RBK114-24"
    assert flags[flags.index("--experiment-duration") + 1] == str(hours)
    assert flags[-1] == "/tmp/.dorado/playback with spaces.fast5"
    assert report["reads"]["intersecting_output"] == 2
    assert report["reads"]["clipped_at_output_boundary"] == 1
    assert report["reads"]["retained_outside_output"] == 2
    assert report["reads"]["excluded_forced"] == 2
    assert report["timeline"]["duration_seconds"] == seconds
    assert report["output_bytes"] == 6


@pytest.mark.parametrize("exclude,kept,removed", [(True, 1, 1), (False, 2, 0)])
def test_index_counts_excluded_reasons_and_source_timing(tmp_path, monkeypatch, exclude, kept, removed):
    info = {"acquisition_id": "run", "sample_rate": 5000}
    monkeypatch.setattr(converter, "run_info_dict", lambda ri: info)

    def read(start, forced):
        return SimpleNamespace(
            start_sample=start, sample_count=100, read_id=str(start),
            end_reason=SimpleNamespace(forced=forced, name="unblock" if forced else "signal_positive"),
            run_info=SimpleNamespace(acquisition_id="run"),
            pore=SimpleNamespace(channel=1, well=1), read_number=start,
            median_before=0, calibration=SimpleNamespace(offset=0, scale=1),
            num_minknow_events=0, tracked_scaling=None, predicted_scaling=None,
            num_reads_since_mux_change=0, time_since_mux_change=0,
        )

    class Reader:
        def __init__(self, path):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def reads(self):
            return iter([read(0, True), read(200, False)])

    monkeypatch.setattr(converter.pod5, "Reader", Reader)
    summary = converter.build_index([tmp_path / "source.pod5"], tmp_path / "index.sqlite", exclude)
    assert (summary["total"], summary["kept"], summary["removed"]) == (2, kept, removed)
    assert summary["source_min_start"] == 0
    assert summary["source_max_end"] == 300
    assert summary["min_start"] == (200 if exclude else 0)
    assert summary["end_reasons"]["unblock"]["excluded"] == removed


@pytest.mark.parametrize("timeline,origin_mode,origin,duration", [
    ("source", "absolute", 0, 1000),
    ("source", "rebase", 100, 900),
    ("retained", "absolute", 0, 800),
    ("retained", "rebase", 200, 600),
])
def test_independent_timeline_and_origin(timeline, origin_mode, origin, duration):
    args = converter.parser().parse_args([
        "convert", "input.pod5", "--output", "out.fast5", "--background-cache", "cache.h5",
        "--timeline", timeline, "--time-origin", origin_mode,
    ])
    summary = {
        "run_info": {"sample_rate": 10}, "source_min_start": 100, "source_max_end": 1000,
        "min_start": 200, "max_end": 800,
    }
    assert converter.resolve_timeline(args, summary) == (origin, duration)
    if timeline == "source":
        summary.update(min_start=400, max_end=600)
        assert converter.resolve_timeline(args, summary) == (origin, duration)


def test_explicit_duration_and_test_cap():
    base = ["convert", "input.pod5", "--output", "out.fast5", "--background-cache", "cache.h5"]
    summary = {
        "run_info": {"sample_rate": 10}, "source_min_start": 100, "source_max_end": 1000,
        "min_start": 200, "max_end": 800,
    }
    parser = converter.parser()
    assert converter.resolve_timeline(parser.parse_args(base), summary) == (0, 1000)
    for seconds in (100, 120):
        args = parser.parse_args(base + ["--duration-seconds", str(seconds)])
        assert converter.resolve_timeline(args, summary) == (0, seconds * 10)
    args = parser.parse_args(base + ["--time-origin", "rebase", "--duration-seconds", "120"])
    assert converter.resolve_timeline(args, summary) == (100, 1200)
    args = parser.parse_args(base + ["--duration-seconds", "99"])
    with pytest.raises(ValueError, match="cannot shorten"):
        converter.resolve_timeline(args, summary)
    args = parser.parse_args(base + ["--max-duration-seconds", "30"])
    assert converter.resolve_timeline(args, summary) == (0, 300)
    args = parser.parse_args(base + ["--max-duration-seconds", "0.001"])
    with pytest.raises(ValueError, match="at least one sample"):
        converter.resolve_timeline(args, summary)
    with pytest.raises(SystemExit):
        parser.parse_args(base + ["--duration-seconds", "120", "--max-duration-seconds", "30"])


@pytest.mark.parametrize("value", ["0", "-1", "nan", "inf", "-inf"])
@pytest.mark.parametrize("flag", ["--duration-seconds", "--max-duration-seconds"])
def test_invalid_duration_arguments(value, flag):
    with pytest.raises(SystemExit):
        converter.parser().parse_args([
            "convert", "input.pod5", "--output", "out.fast5", "--background-cache", "cache.h5",
            f"{flag}={value}",
        ])


def test_excluded_reads_cannot_hide_another_acquisition(tmp_path, monkeypatch):
    monkeypatch.setattr(converter, "run_info_dict", lambda ri: {"acquisition_id": ri.acquisition_id})

    class Reader:
        def __init__(self, path):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def reads(self):
            for run in ("one", "two"):
                yield SimpleNamespace(
                    start_sample=0, sample_count=100,
                    end_reason=SimpleNamespace(forced=True, name="unblock"),
                    run_info=SimpleNamespace(acquisition_id=run),
                )

    monkeypatch.setattr(converter.pod5, "Reader", Reader)
    with pytest.raises(ValueError, match="2 acquisition IDs"):
        converter.build_index([tmp_path / "source.pod5"], tmp_path / "index.sqlite", True)
