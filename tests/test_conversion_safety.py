"""Path-collision tests use disposable byte markers, never lab data."""
import os
from pathlib import Path

import pytest

import pod5_to_bulk_fast5 as converter


def snapshot(root):
    return {str(path.relative_to(root)): (path.read_bytes(), os.readlink(path) if path.is_symlink() else None)
            for path in root.rglob("*") if path.is_file()}


def link(alias, target, kind):
    alias.unlink(missing_ok=True)
    if kind == "symlink":
        alias.symlink_to(target)
    else:
        os.link(target, alias)


@pytest.mark.parametrize("protected", ["input", "cache"])
@pytest.mark.parametrize("artifact", ["output", "report", "pending", "index", "journal", "wal", "shm"])
@pytest.mark.parametrize("alias_kind", ["direct", "symlink", "hardlink"])
def test_convert_rejects_alias_before_any_mutation(tmp_path, monkeypatch, protected, artifact, alias_kind):
    source = tmp_path / "input.pod5"
    cache = tmp_path / "cache.h5"
    output = tmp_path / "output.fast5"
    temp = tmp_path / "tmp"
    temp.mkdir()
    database = temp / "pod5_bulk_index.sqlite"
    artifacts = {
        "output": output,
        "report": Path(str(output) + ".report.json"),
        "pending": Path(str(output) + ".report.json.tmp"),
        "index": database,
        "journal": Path(str(database) + "-journal"),
        "wal": Path(str(database) + "-wal"),
        "shm": Path(str(database) + "-shm"),
    }
    for path in (source, cache, *artifacts.values()):
        path.write_bytes(f"keep {path.name}".encode())
    target = source if protected == "input" else cache
    alias = artifacts[artifact]
    if alias_kind == "direct":
        # The protected input itself occupies an artifact path.
        if protected == "input":
            source = alias
        else:
            cache = alias
    else:
        link(alias, target, alias_kind)

    def indexing_must_not_start(*args):
        pytest.fail("Collision must be rejected before indexing or overwriting files")

    monkeypatch.setattr(converter, "build_index", indexing_must_not_start)
    args = converter.parser().parse_args([
        "convert", str(source), "--background-cache", str(cache), "--output", str(output),
        "--tmp-dir", str(temp), "--force",
    ])
    before = snapshot(tmp_path)
    with pytest.raises(ValueError, match="Path collision"):
        converter.convert(args)
    assert snapshot(tmp_path) == before


@pytest.mark.parametrize("alias_kind", ["direct", "symlink", "hardlink"])
def test_convert_artifacts_cannot_alias_each_other(tmp_path, alias_kind):
    source = tmp_path / "input.pod5"
    source.write_bytes(b"source")
    cache = tmp_path / "cache.h5"
    cache.write_bytes(b"cache")
    output = tmp_path / "output.fast5"
    output.write_bytes(b"old output")
    database = tmp_path / "pod5_bulk_index.sqlite"
    if alias_kind == "direct":
        output = database
        output.write_bytes(b"old index/output")
    else:
        link(database, output, alias_kind)
    args = converter.parser().parse_args([
        "convert", str(source), "--background-cache", str(cache), "--output", str(output),
        "--tmp-dir", str(tmp_path), "--force",
    ])
    before = snapshot(tmp_path)
    with pytest.raises(ValueError, match="Path collision"):
        converter.convert(args)
    assert snapshot(tmp_path) == before


@pytest.mark.parametrize("alias_kind", ["direct", "symlink", "hardlink"])
def test_background_output_cannot_alias_donor(tmp_path, alias_kind):
    donor = tmp_path / "donor.fast5"
    donor.write_bytes(b"do not modify donor")
    output = donor if alias_kind == "direct" else tmp_path / "cache.h5"
    if alias_kind != "direct":
        link(output, donor, alias_kind)
    args = converter.parser().parse_args(["make-background", str(donor), str(output), "--force"])
    before = snapshot(tmp_path)
    with pytest.raises(ValueError, match="Path collision"):
        converter.make_background(args)
    assert snapshot(tmp_path) == before


@pytest.mark.parametrize("flag", ["--channels", "--chunk-samples"])
@pytest.mark.parametrize("value", ["0", "-1", "1.5", "nan"])
def test_convert_rejects_invalid_positive_integers(flag, value):
    with pytest.raises(SystemExit):
        converter.parser().parse_args([
            "convert", "input.pod5", "--output", "out.fast5", "--background-cache", "cache.h5",
            flag, value,
        ])


def test_convert_valid_counts_and_defaults_unchanged():
    base = ["convert", "input.pod5", "--output", "out.fast5", "--background-cache", "cache.h5"]
    args = converter.parser().parse_args(base)
    assert args.channels is None
    assert args.chunk_samples == 180480
    args = converter.parser().parse_args(base + ["--channels", "1", "--chunk-samples", "1"])
    assert (args.channels, args.chunk_samples) == (1, 1)
