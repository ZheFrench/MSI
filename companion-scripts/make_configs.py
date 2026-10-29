#!/usr/bin/env python3
"""
make_configs.py

Writes the params JSON of the given runs for BASH.msi_detection.batch.sh: for each
run, <output-dir>/<run>_params.json = copy of --template with
"msi_detection.inputs".run = <bams-dir>/<run> (which must hold at least one
*.bam) and, when sample IDs are given for the run, "msi_detection.inputs".samples
= those IDs (the workflow then only analyses the BAMs named <ID>.bam, exact name,
every library). The same run name can be in several groups
(input/colon-msi, input/endo-mss...): bams-dir picks the group.

    make_configs.py BAMS_DIR [RUN ...] [--runs-file TXT] --output-dir DIR [--template JSON]

    make_configs.py /media/jp/SSD2/CHU/input/colon-msi Run141_Juno Run123_Juno \\
        --output-dir /media/jp/SSD2/CHU/config/colon-msi
        -> /media/jp/SSD2/CHU/config/colon-msi/Run141_Juno_params.json, ...

--runs-file: a text file, one line per run or per sample (blank lines and #
comments ignored), added to the runs given on the command line:

    Run141_Juno                  every BAM of the run (but the controls)
    Run171_Juno  C22D03556T      only these samples of the run (one line each,
    Run171_Juno  C22D03587T      columns separated by tabs or spaces)

A run with at least one sample ID only gets its IDs; each ID must match a BAM
of the run (<ID>.bam, exact name), checked here before the workflow.

--output-dir (required) is created if needed. Existing params JSON are
overwritten. Exits 1 when a run has no BAM or a sample ID no BAM (the other
runs are still written); the lines of those runs are written, in the
--runs-file format, in <output-dir>/failed-runs.txt (for a rerun), removed when
no run fails.
"""

import argparse
import json
import logging
import sys
from pathlib import Path

TEMPLATE = Path(__file__).resolve().parent.parent / "params" / "template.json"

logger = logging.getLogger("make-configs")


def read_runs_file(path: Path) -> list:
    """(run, sample ID or None) per line of a runs file."""
    entries = []
    for line in path.read_text().splitlines():
        fields = line.split("#")[0].split()
        if len(fields) > 2:
            sys.exit(f"{path}: more than 2 columns (run, sample): {line!r}")
        if fields:
            entries.append((fields[0], fields[1] if len(fields) == 2 else None))
    return entries


def main():
    parser = argparse.ArgumentParser(description="Params JSON of the given runs of a BAM group folder.")
    parser.add_argument("bams_dir", type=Path, help="Group folder holding one BAM folder per run.")
    parser.add_argument("runs", nargs="*", help="Run names (BAM folders of bams_dir), every sample.")
    parser.add_argument("--runs-file", type=Path, help="Text file, one line per run (every sample) or run + sample ID.")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Where <run>_params.json is written (created if needed).")
    parser.add_argument("--template", type=Path, default=TEMPLATE,
                        help="Params JSON whose run input is replaced per run.")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")

    entries = [(run, None) for run in args.runs]
    if args.runs_file:
        entries += read_runs_file(args.runs_file)
    if not entries:
        parser.error("no run (give RUN names or --runs-file)")
    # run -> its sample IDs (empty: every sample), in the order of the lines
    samples_by_run = {}
    for run, sample in entries:
        samples = samples_by_run.setdefault(run, [])
        if sample and sample not in samples:
            samples.append(sample)

    bams_dir = args.bams_dir.resolve()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    template = json.loads(args.template.read_text())

    failed = []
    for run, samples in samples_by_run.items():
        bam_dir = bams_dir / run
        bam_names = [p.name[:-len(".bam")] for p in bam_dir.rglob("*.bam")] if bam_dir.is_dir() else []
        if not bam_names:
            logger.error(f"{run}: no BAM in {bam_dir}")
            failed.append(run)
            continue
        missing = [s for s in samples if s not in bam_names]
        if missing:
            logger.error(f"{run}: no BAM for sample(s) {', '.join(missing)} in {bam_dir}")
            failed.append(run)
            continue
        params = json.loads(json.dumps(template))
        params["msi_detection.inputs"]["run"] = str(bam_dir)
        params["msi_detection.inputs"]["samples"] = samples or [""]
        output = output_dir / f"{run}_params.json"
        output.write_text(json.dumps(params, indent=4) + "\n")
        selected = f"samples {', '.join(samples)}" if samples else "every sample"
        logger.info(f"{run}: {len(bam_names)} BAM(s), {selected} -> {output}")

    failed_file = output_dir / "failed-runs.txt"
    if failed:
        lines = [f"{run}\t{sample}" if sample else run for run, sample in entries if run in failed]
        failed_file.write_text("\n".join(dict.fromkeys(lines)) + "\n")
        logger.error(f"{len(failed)} run(s) failed -> {failed_file}")
    else:
        failed_file.unlink(missing_ok=True)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
