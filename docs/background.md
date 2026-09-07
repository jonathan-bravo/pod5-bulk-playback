# Reconstruction and background signal

[Documentation index](README.md) · [Project overview](../README.md)

Run command examples from the repository root unless stated otherwise.

## Reconstruction limitations

This is a **reconstruction**, not a lossless format conversion. POD5 contains
signal and metadata for selected reads; it does not contain the continuous
between-read signal, ASIC commands, or complete MinKNOW state history found in
a bulk FAST5. The tool therefore:

1. extracts genuine `pore`-state background from a compatible donor bulk FAST5;
2. converts donor ADC values through pA into each destination channel's
   calibration and repeats those background segments;
3. overlays each retained POD5 read at its recorded channel and sample position;
4. reconstructs read, mux, and state tables from POD5 metadata; and
5. leaves `/Device/AsicCommands` empty because those commands cannot be
   recovered from POD5.

MinKNOW's API states that `/Raw/Channel_*/Signal` is the dataset used for
playback. Auxiliary tables are retained/reconstructed for inspection, but they
must not be interpreted as the original instrument history.

## Use or rebuild the background cache

A reusable R10.4.1 5 kHz background cache is already included in this repository
as `background_r10_4_1_5khz.h5`, so normal conversions can skip this step and
use that file directly with `--background-cache`.

To rebuild the cache from the original open-access NA12878 human donor bulk
FAST5, download the R10.4 5 kHz bulk file from the readfish playback dataset:

```bash
curl -L -o GXB02001_20230509_1250_FAW79338_X3_sequencing_run_NA12878_B1_19382aa5_ef4362cd.fast5 \
  https://s3.amazonaws.com/nanopore-human-wgs/bulkfile/GXB02001_20230509_1250_FAW79338_X3_sequencing_run_NA12878_B1_19382aa5_ef4362cd.fast5
```

The cache output must be distinct from the donor bulk FAST5, including symlink
and hard-link aliases. `--force` cannot override this protection.

The donor bulk FAST5 is only required when rebuilding the cache. The command
reads the state tables and short signal slices classified as `pore`; it does not
read the complete 22 GB signal.

```bash
python pod5_to_bulk_fast5.py make-background \
  GXB02001_20230509_1250_FAW79338_X3_sequencing_run_NA12878_B1_19382aa5_ef4362cd.fast5 \
  background_r10_4_1_5khz.h5
```

The default cache contains up to ten seconds per channel. Preserve and share
this small cache with the script and environment file. A donor should match the
flow-cell chemistry and sample rate of the destination data. If a donor channel
has no qualifying `pore` interval, the cache records and uses the nearest donor
channel; `--strict-channels` disables this fallback.

For the files here, both runs are FLO-MIN114 at 5 kHz, but the donor used
SQK-LSK114 while the POD5 run used SQK-RBK114-24. Reusing the donor background
across those preparations follows the project decision that open-pore
background is acceptable across these runs; it is not evidence that the two
library preparations have identical background distributions.

## Sources

- ONT POD5 specification (read `start`, channel, calibration, and signal table):
  https://software-docs.nanoporetech.com/pod5/latest/specification/
- ONT POD5 tools (`pod5 convert to_fast5` produces read FAST5, not bulk FAST5):
  https://software-docs.nanoporetech.com/pod5/latest/tools/
- MinKNOW API bulk configuration (`Signal` is used in playback):
  https://github.com/nanoporetech/minknow_api/blob/master/proto/minknow_api/analysis_configuration.proto
- MinKNOW API playback source setting:
  https://github.com/nanoporetech/minknow_api/blob/master/proto/minknow_api/protocol.proto
- readfish playback instructions and donor dataset:
  https://github.com/LooseLab/readfish#configuring-bulk-fast5-file-playback
- POD5 Python API:
  https://software-docs.nanoporetech.com/pod5/latest/reference/api/reader/
