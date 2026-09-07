# POD5 → bulk FAST5 playback builder

Build a continuous, 512-channel bulk FAST5 for MinKNOW/readfish playback from
read-oriented POD5 files. The tool overlays recorded read signals on reusable
open-pore background and generates a conversion report with recommended playback
flags. A separate audit checks read-signal preservation.

**This is a reconstruction, not a lossless format conversion:** between-read
signal and instrument history cannot be recovered from POD5. See
[reconstruction limits and background compatibility](docs/background.md).
Docker is optional.

## Quick start

Run these commands from the repository root. Replace `/data/run/*.pod5` with
files from **one acquisition**. The included background cache is for compatible
R10.4.1, 5 kHz data.

```bash
conda env create -f environment.yml
conda activate pod5-bulk-playback

python pod5_to_bulk_fast5.py convert '/data/run/*.pod5' \
--background-cache background_r10_4_1_5khz.h5 \
--output bacterial_playback.fast5 \
--tmp-dir tmp
```

By default, all recorded reads—including forced-ended reads—are included, and
timing runs from acquisition start through the latest supplied read's end.
For large runs, [check storage](docs/setup.md#storage-planning) and try a short
window first.

### Short test with complete recorded reads

```bash
python pod5_to_bulk_fast5.py convert '/data/run/*.pod5' \
--background-cache background_r10_4_1_5khz.h5 \
--output test_30s.fast5 \
--timeline retained \
--time-origin rebase \
--max-duration-seconds 30 \
--complete-reads-only
```

### Validate and audit after conversion

```bash
python pod5_to_bulk_fast5.py validate bacterial_playback.fast5
python pod5_to_bulk_fast5.py audit bacterial_playback.fast5 --mode sampled --reads 1000
```

Validation checks structure; the audit compares original POD5 samples against
the output. Keep the source files and `bacterial_playback.fast5.report.json`.
For playback, use the report's `recommended_simulation.shell_flags` with your
installed MinKNOW `start_protocol.py`; see the [playback guide](docs/playback.md).

## Documentation

- [Setup, storage, and tests](docs/setup.md)
- [Reconstruction and background cache](docs/background.md)
- [Conversion defaults, time windows, and reports](docs/conversion.md)
- [Signal-integrity audit and interpreting results](docs/auditing.md)
- [Playback, Docker persistence, and yield troubleshooting](docs/playback.md)

Licensed under [GPL-3.0](LICENSE).
