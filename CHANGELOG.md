# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- Runs with [Cromwell](https://github.com/broadinstitute/cromwell) 92 (Slurm + Apptainer backend,
  `backends.conf/cromwell.conf`) instead of sprocket; the WDL is now version 1.1.
- No `Directory` type in WDL 1.1: `inputs.run` is a `String` (the launchers bind it in the
  containers), and `ta03_report` returns `viz.tar.gz` (reports + `samples.tsv`), extracted by
  `msi_detection_post_p.py index`, which now reads Cromwell's metadata JSON.
- Images are named Docker-style in `params/install.json` (`msisensor-pro:1.3.0`...) and given to
  the tasks as the `apptainer` runtime attribute (Cromwell 92 bug with `container`, issue #7846).
- The launchers run their own Python steps in `python-MSI-3.12.sif`: no Python needed on the host.

### Added
- `companion-scripts/BASH.msi_detection.cromwell_test.sh`: one params JSON by hand, with a summary.

### Removed
- `backends.conf/sprocket.cron.toml`, `backends.conf/sprocket.batch.toml`.

## [1.0.0] - 2026-10-28

### Added
- Initial release.
