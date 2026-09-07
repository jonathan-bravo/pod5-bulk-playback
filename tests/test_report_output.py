"""Exercise report publication using a tiny background-only conversion."""
import importlib.util
import json
import sqlite3
from contextlib import nullcontext
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("report_converter", ROOT / "pod5_to_bulk_fast5.py")
converter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(converter)


def test_report_publication_and_failed_forced_rebuild(tmp_path, monkeypatch):
    source = tmp_path / "input.pod5"
    source.touch()
    output = tmp_path / "output.fast5"
    report_path = tmp_path / "output.fast5.report.json"
    info = {
        "acquisition_id": "test", "sample_rate": 5000,
        "adc_min": -32768, "adc_max": 32767,
        "context_tags": {}, "tracking_id": {},
        "sequencing_kit": "sqk-rbk114-24", "flow_cell_product_code": "FLO-MIN114",
    }

    def index(inputs, database, exclude_forced):
        database.unlink(missing_ok=True)
        with sqlite3.connect(database) as con:
            con.executescript(converter.INDEX_SCHEMA)
        return {
            "run_info": info, "total": 0, "kept": 0, "removed": 0,
            "end_reasons": {}, "source_min_start": 0, "source_max_end": 50,
            "min_start": 0, "max_end": 50, "max_channel": 1,
        }

    monkeypatch.setattr(converter, "build_index", index)
    monkeypatch.setattr(converter.pod5, "Reader", lambda path: nullcontext())
    args = converter.parser().parse_args([
        "convert", str(source), "--output", str(output),
        "--background-cache", str(ROOT / "background_r10_4_1_5khz.h5"),
        "--tmp-dir", str(tmp_path / "index"), "--channels", "1",
        "--compression", "gzip", "--auxiliary", "none",
    ])
    converter.convert(args)
    report = json.loads(report_path.read_text())
    assert report["status"] == "conversion_completed"
    assert report["timeline"]["duration_samples"] == 50
    assert report["output_bytes"] == output.stat().st_size
    assert not Path(str(report_path) + ".tmp").exists()
    with pytest.raises(FileExistsError):
        converter.convert(args)
    assert report_path.exists()

    def fail(*args):
        raise RuntimeError("simulated rebuild failure")

    args.force = True
    monkeypatch.setattr(converter, "build_index", fail)
    with pytest.raises(RuntimeError, match="simulated rebuild failure"):
        converter.convert(args)
    assert not report_path.exists()
