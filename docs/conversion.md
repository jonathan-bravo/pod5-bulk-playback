# Conversion, timing, and reports

[Documentation index](README.md) · [Project overview](../README.md)

Run command examples from the repository root unless stated otherwise.

## Convert all POD5 shards from one run

```bash
python pod5_to_bulk_fast5.py convert '/data/run/*.pod5' \
--background-cache background_r10_4_1_5khz.h5 \
--output bacterial_playback.fast5 \
--tmp-dir tmp
```

Defaults chosen for this project:

- all recorded reads, including forced-ended reads from previous unblocks and
  mux-change terminations, are included;
- the timeline starts at acquisition sample zero and ends at the latest end of
  **all supplied reads**, including excluded reads;
- original acquisition-relative channel/read timing is preserved;
- 512 channels are emitted;
- VBZ compression is used; and
- reconstructed auxiliary tables and device metadata are written.

Multiple POD5 shards do **not** restore a read that MinKNOW already truncated.
A following shard contains later read records, not the missing continuation of
an unblocked read. Excluding forced reads replaces their intervals with donor
open-pore background.

Forced-ended reads are valid captured signal and are now retained by default
(`--no-exclude-forced` makes this explicit). Use `--exclude-forced` only if you
want the previous exclusion behavior. Filtering changes which signals are
overlaid, not the default timeline. Do not combine POD5 files with different
acquisition IDs, even when all reads from one acquisition would be filtered out.

### Full-run versus compact timing

The default `--timeline source --time-origin absolute` preserves acquisition
start through the latest supplied read end. For a partial dataset recorded late
in a run, this can create a large leading background interval. It does **not**
recover missing reads or prove the inputs cover the full experiment.

For compact tests, use the previous timing behavior explicitly:

```bash
--timeline retained --time-origin rebase
```

`--timeline` chooses which reads establish the time bounds (`source`: all input
reads; `retained`: only reads surviving filtering). `--time-origin` independently
chooses sample zero: acquisition start (`absolute`) or the earliest read in the
selected bounds (`rebase`). All relative read/channel timing remains unchanged.

POD5 does not establish the experiment stop time after the last read. If you
know the full duration, supply it in seconds; for example, with the default
absolute origin, a known 22-hour experiment can use:

```bash
--duration-seconds 79200
```

This specifies total output length from the chosen origin, not extra padding.
If a window-start offset is supplied, length is measured from that window start.
It may extend the remaining selected timeline with background but cannot shorten
it. `--max-duration-seconds` instead caps output length for tests. These options
are mutually exclusive. Conversions still fail if no reads remain after
forced-read filtering.

**Compatibility:** the default used to be retained/rebased timing. Add
`--timeline retained --time-origin rebase` to reproduce that behavior.

The metadata index is retained as `tmp/pod5_bulk_index.sqlite` for auditing.
Only one channel's read signals are held in memory at a time.

### Complete-read test windows

Use `--complete-reads-only` to omit reads crossing either output boundary rather
than clipping their signal. For a 30-second window starting ten minutes after
acquisition start, append:

```bash
--window-start-seconds 600 \
--max-duration-seconds 30 \
--complete-reads-only
```

The start offset defaults to zero and is measured from the chosen time origin:
acquisition start for `--time-origin absolute`, or the earliest read in the
selected timeline for `--time-origin rebase`. Seconds are rounded to the nearest
sample. The window must start before the selected timeline ends, and its length
is capped by the remaining timeline unless explicitly extended with
`--duration-seconds`. Source window start becomes sample zero in the output.

A read is included only if it starts at or after the window start and ends at
or before the window end. Exact boundary fits are included. Reads ending at the
window start or starting at its end are outside the window. Skipped reads leave
background behind; reads are not moved together or extended. The same selection
is used for raw signal and reconstructed read/mux/state tables. If no reads fit,
the output is background-only and a warning is printed.

This option is independent of forced-read filtering: recorded unblock/mux-change
reads that fit the window are included by default. Use `--exclude-forced` only
if you explicitly want to omit those reads as well.
"Complete" here means the entire **recorded signal**, not necessarily an entire
molecule. Without `--complete-reads-only`, boundary-crossing reads are still
clipped, preserving the previous test-window behavior.

## Small verification run

Before creating a many-hour output, test 512 channels and 30 seconds:

```bash
python pod5_to_bulk_fast5.py convert ../FAX54151_4f394667_45b25379_8.pod5 \
--background-cache background_r10_4_1_5khz.h5 \
--output tmp/test_512ch_30s.fast5 \
--tmp-dir tmp \
--channels 512 \
--timeline retained \
--time-origin rebase \
--max-duration-seconds 30 \
--complete-reads-only

python pod5_to_bulk_fast5.py validate tmp/test_512ch_30s.fast5
```

`validate` checks structure and dimensions without scanning signal values. Use
the [signal-integrity audit](auditing.md) to check recorded read samples. The
final acceptance test must be [playback](playback.md) in the intended
MinKNOW/readfish environment and protocol configuration, whether containerized
or not.

## Automatic conversion report and playback flags

Every successful conversion writes an indented JSON report beside the FAST5:
`bacterial_playback.fast5.report.json`. No extra flag is required. The report
contains input paths and sizes, run metadata, total/retained/excluded read
counts, counts by end reason, source and retained timing, output duration,
conversion settings, and recommended MinKNOW simulation flags.

Retained index reads and reads intersecting the actual output are reported
separately, including reads selected for overlay, boundary-crossing reads,
reads excluded at the boundary, reads actually clipped, and reads outside the
requested channel/time window. Boundary exclusion counts are separate from
forced-read exclusions. These are metadata counts, not proof of signal integrity
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
The flags are also printed when conversion finishes. See the
[playback guide](playback.md) for output persistence, logs, and yield comparisons.

The report is published only after the FAST5 has been written and closed.
Existing output/report files require `--force`; a forced rebuild removes the
old report first so a failed rebuild cannot leave a stale success report.
Reports are generated for new conversions only; old indexes cannot recover
previously discarded read counts. Reports include source metadata and paths,
so review them before sharing.
