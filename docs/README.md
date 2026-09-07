# Documentation

[Project overview and quick start](../README.md)

Command examples use the repository root as the working directory unless stated
otherwise. Docker is optional; container-specific examples are clearly labeled.

| Guide | Contents |
|---|---|
| [Setup](setup.md) | Conda environment, storage planning, and tests |
| [Background](background.md) | Reconstruction limits, donor compatibility, cache rebuilding, and sources |
| [Conversion](conversion.md) | Read filtering, full-run timing, complete-read windows, and conversion reports |
| [Signal audit](auditing.md) | Post-conversion sampled/full comparison, overlaps, and exit codes |
| [Playback](playback.md) | Simulation settings, logs, Docker persistence, and yield comparisons |

Suggested workflow: set up → verify background compatibility → convert a small
window → validate and audit → test playback → convert and replay the full run.
