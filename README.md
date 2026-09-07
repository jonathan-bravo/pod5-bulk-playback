# POD5 → bulk FAST5 playback builder

This project builds a continuous, 512-channel bulk FAST5 for MinKNOW/readfish
playback from a complete set of read-oriented POD5 files.

## Important limitation

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

## Reproducible environment

```bash
conda env create -f environment.yml
conda activate pod5-bulk-playback
```

To update an existing environment after `environment.yml` changes:

```bash
conda env update -n pod5-bulk-playback -f environment.yml --prune
```

## 1. Use or rebuild the reusable background cache

A reusable R10.4.1 5 kHz background cache is already included in this repository
as `background_r10_4_1_5khz.h5`, so normal conversions can skip this step and
use that file directly with `--background-cache`.

To rebuild the cache from the original open-access NA12878 human donor bulk
FAST5, download the R10.4 5 kHz bulk file from the readfish playback dataset:

```bash
curl -L -o GXB02001_20230509_1250_FAW79338_X3_sequencing_run_NA12878_B1_19382aa5_ef4362cd.fast5 \
  https://s3.amazonaws.com/nanopore-human-wgs/bulkfile/GXB02001_20230509_1250_FAW79338_X3_sequencing_run_NA12878_B1_19382aa5_ef4362cd.fast5
```

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

## 2. Convert all POD5 shards from one run

```bash
python pod5_to_bulk_fast5.py convert '/data/run/*.pod5' \
  --background-cache background_r10_4_1_5khz.h5 \
  --output bacterial_playback.fast5 \
  --tmp-dir tmp
```

Defaults chosen for this project:

- all forced reads (including previous unblocks and mux-change terminations) are
  excluded;
- the earliest retained read is rebased to sample zero;
- original relative channel/read timing is preserved;
- 512 channels are emitted;
- VBZ compression is used; and
- reconstructed auxiliary tables and device metadata are written.

Multiple POD5 shards do **not** restore a read that MinKNOW already truncated.
A following shard contains later read records, not the missing continuation of
an unblocked read. Excluding forced reads replaces their intervals with donor
open-pore background.

Use `--no-exclude-forced` to retain them. Use `--time-origin absolute` only when
leading time from acquisition start is important; it can make the output much
larger. Do not combine POD5 files with different acquisition IDs.

The metadata index is retained as `tmp/pod5_bulk_index.sqlite` for auditing.
Only one channel's read signals are held in memory at a time.

## Automatic conversion report and playback flags

Every successful conversion writes an indented JSON report beside the FAST5:
`bacterial_playback.fast5.report.json`. No extra flag is required. The report
contains input paths and sizes, run metadata, total/retained/excluded read
counts, counts by end reason, source and retained timing, output duration,
conversion settings, and recommended MinKNOW simulation flags.

Retained index reads and reads intersecting the actual output are reported
separately, including boundary-clipped reads and reads outside the requested
channel/time window. These are metadata counts, not proof of signal integrity
or successful playback. Source timing ends at the latest supplied read, not
necessarily the original experiment stop time.

The `recommended_simulation.shell_flags` field is ready to append to your
installed `start_protocol.py` command. For the 21.03-hour bacterial example:

```bash
--position MS00000 \
--product-code FLO-MIN114 \
--kit SQK-RBK114-24 \
--experiment-duration 22 \
--fastq --bam --pod5 --verbose --basecalling \
--simulation /tmp/.dorado/bacterial_playback.fast5
```

Kit and product code are taken from POD5 metadata and uppercased. Missing
values become explicit placeholders. Duration is computed from the **generated
FAST5**, rounded up to whole hours (minimum one hour). Thus a 30-second test
recommends one hour, not 22. This flag sets a protocol time limit; it does not
guarantee how playback behaves at end-of-file. In particular,
`--experiment-duration 1` stops an overnight replay after one hour.

`MS00000` and `/tmp/.dorado/<output filename>` are editable simulation defaults,
not values inferred from the original device. Verify those settings and that
the installed MinKNOW version supports the source kit/product combination.
The flags are also printed when conversion finishes.

The report is published only after the FAST5 has been written and closed.
Existing output/report files require `--force`; a forced rebuild removes the
old report first so a failed rebuild cannot leave a stale success report.
Reports are generated for new conversions only; old indexes cannot recover
previously discarded read counts. Reports include source metadata and paths,
so review them before sharing.

## Small verification run

Before creating a many-hour output, test 512 channels and 30 seconds:

```bash
python pod5_to_bulk_fast5.py convert ../FAX54151_4f394667_45b25379_8.pod5 \
  --background-cache background_r10_4_1_5khz.h5 \
  --output tmp/test_512ch_30s.fast5 \
  --tmp-dir tmp \
  --channels 512 \
  --max-duration-seconds 30

python pod5_to_bulk_fast5.py validate tmp/test_512ch_30s.fast5
```

`validate` checks structure and dimensions without scanning signal values. The
final acceptance test must be playback in the same Ubuntu/MinKNOW container and
protocol configuration used by readfish.

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

## Storage planning

A 5 kHz, 512-channel run contains 5,120,000 raw bytes per second before
compression. Output size depends strongly on run duration and compressibility.
Check free space before converting all 22 shards and retain the POD5 sources.

## License

This project is licensed under the GNU General Public License v3.0. See [LICENSE](LICENSE).

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
