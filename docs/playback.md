# Playback and troubleshooting

[Documentation index](README.md) · [Project overview](../README.md)

Conversion and signal auditing do not require Docker, MinKNOW, or Dorado.
Playback does require an appropriate MinKNOW simulated-device environment and,
when basecalling is enabled, a configured basecaller. Docker is one way to host
that environment, not a requirement. These examples are guidance for an already
configured installation, not a universal container startup recipe.

## Before playback

1. Run `validate` and preferably the [signal-integrity audit](auditing.md).
2. Check the conversion report's source kit, product code, and generated duration.
3. Confirm the installed MinKNOW protocol supports that kit/product combination.
4. Make the FAST5 accessible at the path used by the playback process.
5. Choose persistent output and log locations, with enough free space for any
   requested FASTQ, BAM, and POD5 outputs.

A running Dorado basecall server is not itself a playback run: it receives signal
from a client. MinKNOW manager services can also remain running after sequencing
has stopped. If readfish is making adaptive-sampling decisions, record its
configuration; intentional unblocks can change recovered yield.

## Use the generated simulation flags

Each conversion writes `<output>.fast5.report.json` and prints recommended flags.
For example, from the repository root in the conversion environment:

```bash
python -c 'import json; print(json.load(open("bacterial_playback.fast5.report.json"))["recommended_simulation"]["shell_flags"])'
```

Append those flags to the command for your **installed** MinKNOW API
`start_protocol.py`. Script paths and setup steps vary by installation/version.
For the bacterial example, the flags are:

```bash
--position MS00000 \
--product-code FLO-MIN114 \
--kit SQK-RBK114-24 \
--experiment-duration 22 \
--fastq --bam --pod5 --verbose --basecalling \
--simulation /tmp/.dorado/bacterial_playback.fast5
```

This block is a set of arguments, not a standalone shell command. The kit and
product come from source metadata. `MS00000` and `/tmp/.dorado/...` are editable
assumptions: use the actual simulated position and a path visible to MinKNOW.
Do not substitute a native barcoding kit for a rapid barcoding kit merely because
both support the same number of barcodes.

**`--experiment-duration` is in hours**, not a request to play an entire file.
The recommendation rounds generated signal duration up to a whole hour (minimum
one hour). A 21.03-hour file recommends 22; a 45-minute file recommends 1.
A one-hour protocol limit will stop an overnight replay early. Verify the help
for your installed script and use the logs to confirm the stop reason.
End-of-file behavior is playback-version dependent; the time limit is not a
guarantee that playback will end automatically at EOF.

The requested flags enable output formats but do not by themselves select the
output directory or guarantee demultiplexing. Configure output paths and barcode
processing through the protocol options supported by your installation.
See [conversion report details](conversion.md#automatic-conversion-report-and-playback-flags).

## Preserve Docker outputs and logs

Skip this section if you are not using Docker. Run the following commands on the
Docker host, including when connected to that host over SSH—not inside the
container.

**A container started with `--rm` loses its writable filesystem when it exits.**
Without that flag, removing the container still deletes that filesystem.
Do not stop/remove a disposable container until important data is copied out.

### Prefer bind mounts for new runs

Create persistent host directories before starting the container:

```bash
mkdir -p "$PWD/playback-output" "$PWD/playback-logs/minknow" "$PWD/playback-logs/dorado"
```

Add bind mounts to your existing `docker run` command, adapting the input path:

```bash
--mount "type=bind,src=/absolute/path/to/input-directory,dst=/playback-input,readonly" \
--mount "type=bind,src=$PWD/playback-output,dst=/playback-output" \
--mount "type=bind,src=$PWD/playback-logs/minknow,dst=/var/log/minknow" \
--mount "type=bind,src=$PWD/playback-logs/dorado,dst=/var/log/dorado"
```

These are arguments to `docker run`, not a standalone command. The source input
directory must exist on the Docker host. Use
`--simulation /playback-input/bacterial_playback.fast5` in this layout. Configure
the protocol to write results under `/playback-output`; **mounting that directory
does not automatically redirect outputs**. Check the actual log paths in your
image before mounting them and give the container's service user appropriate
write permissions on the host output/log directories.

A read-only input mount also avoids copying a large FAST5 into the container's
writable layer. Conversion and auditing inside a container likewise need access
to their input files; auditing resolves the POD5 paths stored in its conversion
report, so keep those paths accessible in the audit environment.

### Rescue data from an existing container

Before stopping/removing it:

```bash
C=your_container_name
mkdir -p playback-backup

docker logs --timestamps "$C" > playback-backup/docker.log 2>&1
docker cp "$C":/var/log/minknow playback-backup/minknow-logs
docker cp "$C":/var/log/dorado playback-backup/dorado-logs
docker cp "$C":/actual/protocol/output/directory playback-backup/sequencing-output
```

Replace the last path with the real protocol output directory and verify the
copied files before deleting anything. Log directories vary by image. Copy the
FAST5 and conversion/audit reports too if their only copies are in the container.

## Diagnose short playback or unexpectedly low yield

### Check timing and termination first

From the repository root in the conversion environment:

```bash
python pod5_to_bulk_fast5.py validate /path/to/bacterial_playback.fast5
```

Compare `duration_seconds` against the protocol's configured runtime and actual
start/stop timestamps. Wall-clock runtime, acquisition time, and time spent
finishing basecalls are different quantities.

For Docker, inspect current processes and capture logs on the host:

```bash
C=your_container_name
docker top "$C" -eo pid,ppid,etime,args
docker logs --timestamps "$C" > playback-docker.log 2>&1
docker exec "$C" df -h / /tmp /var/log
docker exec "$C" sh -lc 'grep -Ein "runtime|starting_playback|protocol_stop|exit_reason|stop_reason|finished|error" /var/log/minknow/MS00000/control_server_log-*.txt | tail -100'
```

Adjust the position/log path. Outside Docker, inspect the corresponding local
processes, log files, and filesystems directly. Check the filesystem holding the
actual sequencing output too, especially when it is separately mounted.

`docker logs` captures the container main process's stdout/stderr; commands
started in interactive `docker exec` sessions may log elsewhere. Useful paths
in the tested image were:

- `/var/log/minknow/MS00000/control_server_log-0.txt`
- `/var/log/minknow/basecall_manager_log-0.txt`
- `/var/log/minknow/mk_manager_svc_log-0.txt`
- `/var/log/dorado/`

Look for the runtime stop criterion, playback start, stop request, and final
reason. For example, `runtime=3600` with a stop approximately one hour after
playback starts explains a short replay even when the FAST5 contains 21 hours.
Inspect surrounding lines for warnings/errors rather than relying on grep alone;
a successful protocol exit does not certify every signal or output record.

### Compare reads and bases, not just file sizes

A bulk FAST5 stores continuous raw signal across channels, including reconstructed
background. FASTQ stores called bases and qualities. **Their byte sizes are not
a yield comparison.** Compressed and uncompressed FASTQ sizes are not directly
comparable either; the tested MinKNOW protocol enabled FASTQ compression.

Compare these quantities instead:

- Passing and failing read counts and total bases.
- Read-length distributions and, where available, quality scores.
- All output shards/barcodes combined—not one file or barcode directory.
- Equivalent source/replay time windows, with the same basecalling model,
  quality threshold, trimming, and demultiplexing settings where possible.

For FASTQ files, an optional external tool is `seqkit` (not installed by this
project's environment). On a Bash/Linux host, summarize every FASTQ beneath an
output directory, including compressed files:

```bash
find /path/to/fastq-output -type f \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) -print0 |
xargs -0 -r seqkit stats -a -T > fastq-stats.tsv
```

This produces per-file statistics; sum read/base counts across the relevant
files. Do not include a pre-combined FASTQ alongside its constituent shards,
since that double-counts reads. Keep pass/fail summaries separate when comparing
quality outcomes. Sequencing summaries and final basecaller/writer log counters
are alternatives when available; verify actual output-file counts if writer
errors were reported.

Replay re-detects read boundaries and may split or merge records. A count mismatch
is not automatically proof of conversion damage. Conversely, a high sequence
mapping percentage does not prove exact preservation of each original read—many
bacterial reads overlap the same genomic regions. Use the
[signal-integrity audit](auditing.md) to test conversion directly, then investigate
protocol settings and downstream processing separately.
