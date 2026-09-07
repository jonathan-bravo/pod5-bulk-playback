# Post-conversion signal-integrity audit

[Documentation index](README.md) · [Project overview](../README.md)

Run command examples from the repository root unless stated otherwise.

Run this **after conversion**, in the conversion environment on the machine
where the source POD5 paths in the report are accessible. It does not require
MinKNOW, Dorado, or Docker.

Start with a reproducible sample (defaults: 1,000 reads, seed 42):

```bash
python pod5_to_bulk_fast5.py audit bacterial_playback.fast5 \
--mode sampled \
--reads 1000 \
--seed 42
```

Then, for exhaustive comparison of every selected read interval:

```bash
python pod5_to_bulk_fast5.py audit bacterial_playback.fast5 \
--mode full \
--force
```

`--force` replaces an existing **audit report only**, not the FAST5 or conversion
report. Results are written to `bacterial_playback.fast5.audit.json`. The audit
uses `<fast5>.report.json` by default; `--report PATH` selects another conversion
report. Old conversions without a report are not supported by this command.

### What is checked

- Rebuild source metadata in a private temporary index and verify it against the
  conversion report: source sizes, run metadata, counts, filtering, and timeline.
  The reusable `tmp/pod5_bulk_index.sqlite` is deliberately not trusted or changed.
- Check output channel dimensions, sample rate, and signed 16-bit signal type.
- Recover the expected retained intervals, respecting forced-read filtering,
  channel limits, clipping, and complete-read windows.
- Compare source ADC values **exactly** against their output locations. Calibration
  metadata is checked for channels with audited reads. Clipped reads are compared
  only over their expected output interval.
- Scan **all selected read intervals** for same-channel overlaps, even in sampled
  mode. Overlapping reads are reported as conflicts, including when their shared
  samples agree. A conflicting overlap is not silently counted as successful
  preservation of every original read. The report includes the number and
  fraction of selected reads involved, plus example overlaps.

Sampled mode prioritizes representatives across input files, channels, and 24
relative time bins, as well as clipped boundary reads. Small budgets cannot
necessarily cover every group. The seed and selected-read digest make a sample
reproducible for the same inputs. Full mode ignores `--reads`.

The report includes reads/samples checked, mismatch counts, and up to 100 detailed
issues with source path, read ID, channel, and first differing sample/value.
Samples are counted per read, so overlapping samples can be counted more than
once. Overlap events are not a count of every overlapping pair.

### Results and limitations

| Status | Meaning | Exit code |
|---|---|---|
| `passed` | Every eligible read interval checked; no mismatch or overlap conflict | 0 |
| `passed_sampled` | Selected sample passed; unchecked reads are not certified | 0 |
| `failed` | Signal mismatch, missing retrieved read, or calibration mismatch | 1 |
| `overlap_conflict` | Selected comparisons matched, but overlapping intervals exist | 1 |
| `inconclusive` | No eligible read signal to compare (e.g. background-only test) | 2 |
| `error` | Audit could not complete, e.g. missing source or inconsistent report | 2 |

Failures before audit setup, such as refusing to overwrite an audit report, use
the CLI's general error exit code 1. When mismatches and overlaps coexist, the
status is `failed` and overlap details are still retained.

Both modes scan all source **metadata**, but sampled mode loads only selected
read signals. Full mode can require substantial I/O. Memory is bounded by one
source read plus small batches and FAST5 slices; `--chunk-samples` controls slice
size. The temporary metadata index uses disk space and is cleaned up afterward.
Use `--tmp-dir /existing/scratch/directory` to choose its parent directory.

The FAST5, POD5 sources, conversion report, and original index are read-only.
Keep source files unchanged during auditing. Source identity is checked using
paths, sizes and metadata, **not** conversion-time cryptographic hashes. This
checks read-signal preservation, not synthetic background, auxiliary state
history, MinKNOW read recovery, or biological/basecalling accuracy.
