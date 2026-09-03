from pathlib import Path

import h5py
import numpy as np
import pytest
import vbz_h5py_plugin

vbz_h5py_plugin.register_plugin()

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
CACHE = ROOT / "background_r10_4_1_5khz.h5"


def test_background_cache_is_complete():
    with h5py.File(CACHE, "r") as handle:
        assert handle.attrs["format"] == "pod5-bulk-background-v1"
        assert int(handle.attrs["sample_rate"]) == 5000
        assert len(handle["background"]) == 512
        for channel in (1, 3, 256, 512):
            group = handle[f"background/Channel_{channel}"]
            assert any(name.startswith("segment_") for name in group)


def test_512_channel_playback_fixture_structure():
    path = DATA / "test_512ch_30s.fast5"
    with h5py.File(path, "r") as handle:
        assert int(handle["Meta"].attrs["sample_rate"]) == 5000
        assert int(handle["Meta"].attrs["duration_samples"]) == 150_000
        assert len(handle["Raw"]) == 512
        assert handle["Raw/Channel_1/Signal"].shape == (150_000,)
        assert handle["Raw/Channel_512/Signal"].shape == (150_000,)
        for group in ("IntermediateData", "MultiplexData", "StateData", "Device"):
            assert group in handle
        assert handle["Device/MetaData"].shape == (150_000,)
        # ASIC command history is not represented by POD5 and must not be invented.
        assert handle["Device/AsicCommands"].shape == (0,)


def test_small_gzip_and_vbz_fixtures_open():
    expected = {
        "test_8ch_30s.fast5": (8, 150_000),
        "test_vbz.fast5": (2, 5_000),
    }
    for filename, (channel_count, duration) in expected.items():
        with h5py.File(DATA / filename, "r") as handle:
            assert len(handle["Raw"]) == channel_count
            assert int(handle["Meta"].attrs["duration_samples"]) == duration
            # Read a tiny slice to verify the installed compression filter works.
            assert handle["Raw/Channel_1/Signal"][:16].dtype == np.int16


def test_embedded_signal_matches_local_pod5_when_available():
    """Optional integration check using the large source file kept outside Git."""
    pod5_path = ROOT.parent / "FAX54151_4f394667_45b25379_8.pod5"
    if not pod5_path.exists():
        pytest.skip("large local POD5 integration input is not available")

    import pod5

    read_id = "15bef6f0-6a08-4b12-92ae-d84a3602fcdd"
    channel = 5
    output_start = 13_436
    sample_count = 116_877
    with pod5.Reader(pod5_path) as reader:
        source = next(reader.reads(selection=[read_id])).signal
    with h5py.File(DATA / "test_512ch_30s.fast5", "r") as handle:
        reconstructed = handle[f"Raw/Channel_{channel}/Signal"][
            output_start : output_start + sample_count
        ]
    assert np.array_equal(source, reconstructed)
