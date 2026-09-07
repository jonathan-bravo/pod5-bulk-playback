# Setup, storage, and tests

[Documentation index](README.md) · [Project overview](../README.md)

Run command examples from the repository root unless stated otherwise.

## Environment

Docker is optional. Conversion, validation, and signal-integrity auditing can
run directly in the Conda environment below or inside a container with the same
dependencies and access to the required files.

```bash
conda env create -f environment.yml
conda activate pod5-bulk-playback
```

To update an existing environment after `environment.yml` changes:

```bash
conda env update -n pod5-bulk-playback -f environment.yml --prune
```

## Storage planning

A 5 kHz, 512-channel run contains 5,120,000 raw bytes per second before
compression. Output size depends strongly on run duration and compressibility.
Check free space before converting all 22 shards and retain the POD5 sources.

## Test suite

The repository includes generated fixtures under `tests/data/`. Because the
FAST5 fixtures are binary files, Git LFS tracks them.

```bash
git lfs install
git lfs pull
conda activate pod5-bulk-playback
pytest -q
```

The tests check the reusable cache, 512-channel bulk structure, gzip and VBZ
readability, reconstructed auxiliary groups, and the intentionally empty ASIC
command history. When the original POD5 remains in the parent directory, an
optional integration test also confirms that a complete 116,877-sample read is
bit-for-bit identical in the generated bulk FAST5.
