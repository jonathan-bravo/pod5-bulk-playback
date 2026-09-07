"""Tiny real POD5/FAST5 fixtures with deterministic, non-biological ADC ramps."""
import datetime as dt
import hashlib
import json
from pathlib import Path
import sqlite3
import uuid

import h5py
import numpy as np
import pod5
import pytest

import pod5_to_bulk_fast5 as converter
import signal_integrity_audit as audit_module

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "background_r10_4_1_5khz.h5"
pytestmark = pytest.mark.filterwarnings("ignore:.*Scaling fields were unused.*:DeprecationWarning")


def make_conversion(tmp_path, flags=(), overlap=False, different_overlap=False):
    now = dt.datetime(2026, 1, 1, tzinfo=dt.timezone.utc)
    info = pod5.RunInfo(
        acquisition_id="test-run", acquisition_start_time=now, adc_max=32767, adc_min=-32768,
        context_tags={}, experiment_name="audit-test", flow_cell_id="test", flow_cell_product_code="FLO-MIN114",
        protocol_name="test", protocol_run_id="test-protocol", protocol_start_time=now, sample_id="test",
        sample_rate=5000, sequencing_kit="sqk-rbk114-24", sequencer_position="test",
        sequencer_position_type="MinION", software="test", system_name="test", system_type="test", tracking_id={},
    )
    layout = [(1, 10, 100, False), (1, 200, 100, True), (1, 450, 100, False),
              (2, 20, 150, False), (2, 250, 100, False), (2, 600, 100, False)]
    if overlap:
        layout.append((1, 50, 100, False))
    reads = []
    for number, (channel, start, length, forced) in enumerate(layout, 1):
        signal = np.arange(start, start + length, dtype=np.int16)
        if different_overlap and number == len(layout):
            signal += 1000
        reads.append(pod5.Read(
            read_id=uuid.UUID(int=number), pore=pod5.Pore(channel, 1, "not_set"),
            calibration=pod5.Calibration(0, 1), read_number=number, start_sample=start, median_before=0,
            end_reason=pod5.EndReason(pod5.EndReasonEnum.UNBLOCK_MUX_CHANGE if forced else pod5.EndReasonEnum.SIGNAL_POSITIVE, forced),
            run_info=info, signal=signal,
        ))
    sources = [tmp_path / "a.pod5", tmp_path / "b.pod5"]
    for source, records in zip(sources, (reads[:3], reads[3:])):
        with pod5.Writer(source) as writer:
            writer.add_reads(records)
    output = tmp_path / "bulk.fast5"
    args = converter.parser().parse_args([
        "convert", *map(str, sources), "--output", str(output), "--background-cache", str(CACHE),
        "--tmp-dir", str(tmp_path / "conversion-index"), "--channels", "2", "--compression", "gzip",
        "--auxiliary", "none", "--chunk-samples", "31", *flags,
    ])
    converter.convert(args)
    return output, sources


def run_audit(output, *options, expected="passed"):
    args = converter.parser().parse_args(["audit", str(output), "--mode", "full", "--chunk-samples", "17", *options])
    if expected in ("passed", "passed_sampled"):
        args.func(args)
    else:
        with pytest.raises(SystemExit) as caught:
            args.func(args)
        assert caught.value.code == (2 if expected in ("inconclusive", "error") else 1)
    result = json.loads(Path(str(output) + ".audit.json").read_text())
    assert result["status"] == expected, result
    return result


@pytest.mark.parametrize("flags,checked", [
    ([], 6),
    (["--exclude-forced"], 5),
    (["--window-start-seconds", "0.01", "--max-duration-seconds", "0.05"], 4),
    (["--window-start-seconds", "0.01", "--max-duration-seconds", "0.05", "--complete-reads-only"], 1),
    (["--timeline", "retained", "--time-origin", "rebase"], 6),
    (["--duration-seconds", "1"], 6),
])
def test_real_signal_audit_and_read_only_inputs(tmp_path, flags, checked):
    output, sources = make_conversion(tmp_path, flags)
    index = tmp_path / "conversion-index" / "pod5_bulk_index.sqlite"
    # Deliberately replace the reusable conversion index. Audit must not use it.
    index.write_bytes(b"unrelated later conversion")
    paths = [output, *sources, Path(str(output) + ".report.json"), index]
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result = run_audit(output)
    assert result["comparison"]["checked_reads"] == checked
    assert result["comparison"]["matching_reads"] == checked
    assert result["overlaps"]["reads_involved"] == 0
    assert before == {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def test_sampled_reproducibility(tmp_path):
    output, _ = make_conversion(tmp_path)
    first = run_audit(output, "--mode", "sampled", "--reads", "2", expected="passed_sampled")
    second = run_audit(output, "--mode", "sampled", "--reads", "2", "--force", expected="passed_sampled")
    assert first["selection"] == second["selection"]
    assert first["comparison"]["checked_reads"] == 2
    args = converter.parser().parse_args(["audit", str(output)])
    with pytest.raises(FileExistsError):
        args.func(args)


@pytest.mark.parametrize("kind", ["signal", "calibration", "shape"])
def test_corruption(tmp_path, kind):
    output, _ = make_conversion(tmp_path)
    with h5py.File(output, "r+") as bulk:
        if kind == "signal":
            bulk["Raw/Channel_1/Signal"][10] = -123
        elif kind == "calibration":
            bulk["Raw/Channel_1/Meta"].attrs["offset"] = 7
        else:
            bulk["Raw/Channel_1/Signal"].resize((699,))
    result = run_audit(output, expected="error" if kind == "shape" else "failed")
    if kind == "signal":
        assert result["comparison"]["mismatched_samples"] == 1
        first = result["comparison"]["details"][0]["first_difference"]
        assert first == {"source_sample": 10, "bulk_sample": 10, "expected_adc": 10, "actual_adc": -123}
    elif kind == "calibration":
        assert result["comparison"]["calibration_mismatches"] == 1


@pytest.mark.parametrize("different", [False, True])
def test_overlaps_not_silently_passed(tmp_path, different):
    output, _ = make_conversion(tmp_path, overlap=True, different_overlap=different)
    result = run_audit(output, expected="failed" if different else "overlap_conflict")
    assert result["overlaps"]["reads_involved"] == 2
    if not different:
        sampled = run_audit(output, "--force", "--mode", "sampled", "--reads", "1", expected="overlap_conflict")
        assert sampled["overlaps"]["reads_involved"] == 2


def test_missing_source_and_stale_report(tmp_path):
    output, sources = make_conversion(tmp_path)
    source = sources[0]
    source.rename(source.with_suffix(".moved"))
    result = run_audit(output, expected="error")
    assert "No such file" in result["error"]
    source.with_suffix(".moved").rename(source)
    report_path = Path(str(output) + ".report.json")
    report = json.loads(report_path.read_text())
    report["timeline"]["output_origin_sample"] = 1
    report_path.write_text(json.dumps(report))
    result = run_audit(output, "--force", expected="error")
    assert "timeline differs" in result["error"]


def test_missing_retrieved_read_is_not_a_pass(tmp_path, monkeypatch):
    output, _ = make_conversion(tmp_path)
    original = pod5.Reader

    class MissingReadReader:
        def __init__(self, path):
            self.reader = original(path)

        def __enter__(self):
            self.reader.__enter__()
            return self

        def __exit__(self, *args):
            return self.reader.__exit__(*args)

        def reads(self, selection=None):
            records = self.reader.reads(selection=selection) if selection is not None else self.reader.reads()
            return (r for r in records if selection is None or str(r.read_id) != str(uuid.UUID(int=1)))

    monkeypatch.setattr(pod5, "Reader", MissingReadReader)
    result = run_audit(output, expected="failed")
    assert result["comparison"]["missing_reads"] == 1
    assert result["comparison"]["checked_reads"] == 5


def test_source_size_change_is_detected(tmp_path):
    output, sources = make_conversion(tmp_path)
    with sources[0].open("ab") as handle:
        handle.write(b"changed")
    result = run_audit(output, expected="error")
    assert "Source size differs" in result["error"]


def test_background_only_is_inconclusive(tmp_path):
    output, _ = make_conversion(tmp_path, ["--max-duration-seconds", "0.001", "--complete-reads-only"])
    result = run_audit(output, expected="inconclusive")
    assert result["comparison"]["checked_samples"] == 0


def test_refuse_output_symlink(tmp_path):
    output, sources = make_conversion(tmp_path)
    sidecar = Path(str(output) + ".audit.json")
    sidecar.symlink_to(sources[0])
    before = sources[0].read_bytes()
    args = converter.parser().parse_args(["audit", str(output), "--force"])
    with pytest.raises(ValueError, match="symlink"):
        args.func(args)
    assert sources[0].read_bytes() == before


@pytest.mark.parametrize("flag", ["--reads", "--chunk-samples"])
@pytest.mark.parametrize("value", ["0", "-1"])
def test_positive_audit_arguments(flag, value):
    with pytest.raises(SystemExit):
        converter.parser().parse_args(["audit", "x.fast5", flag, value])


def test_overlap_sweep_nested_and_touching():
    with sqlite3.connect(":memory:") as con:
        con.execute("CREATE TABLE reads(file_id,read_id,channel,start,samples)")
        con.executemany("INSERT INTO reads VALUES (?,?,?,?,?)", [
            (0, "long", 1, 0, 100), (0, "nested1", 1, 10, 10), (0, "nested2", 1, 30, 10),
            (0, "touching", 1, 100, 10), (0, "another-channel", 2, 0, 100),
        ])
        settings = converter.parser().parse_args(["convert", "x", "--output", "x", "--background-cache", "x"])
        result = audit_module.prepare_index(con, converter, settings, 0, 110, 2)
        assert result["reads_involved"] == 3
        assert result["overlap_events"] == 2
