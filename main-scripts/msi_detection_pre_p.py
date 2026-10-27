#!/usr/bin/env python3
"""
msi_detection_pre_p.py

Pre-processing of the msi_detection-0.1.0.wdl workflow: picks the run(s) left
to analyse and writes, for each, the params JSON (config) sprocket runs it
with. Steps:

    1-2. The runs to analyse: those ADIVAR has finished (status 3 in
       --locker-adivar) and not yet in --locker-msi (a text file, one run
       name per line: the runs already analysed): one run (the new one) or
       none; more is an error (exit 1, nothing written). --runs-dir is not
       listed.
    3. For each run of the difference, the [Header] of its
       <runs-dir>/<run>/SampleSheet.csv
       must have Assay = Juno and Description = ADIVaR_Juno_INCa_v3
       (otherwise the run is skipped); its Experiment Name is the only link
       to its BAMs: the folder <bams-dir>/<Experiment Name>/, which must
       hold at least one *.bam (the workflow's find_bams lists them all).
    4. <output-dir>/<date>-<run>-<Experiment Name>.json (date:
       YYYY-MM-DD_HH-MM-SS): copy of --template with
       "msi_detection.inputs".run = <bams-dir>/<Experiment Name>.
    5. The run is added to --locker-msi as soon as its params JSON is
       written, before it runs: a later launch (cron) never takes it again,
       nor is a failed run retried (remove its line to rerun it).

    msi_detection_pre_p.py --runs-dir RUNS --locker-adivar adivar-testNS.json \\
        --locker-msi msi-sensor.locker.txt --bams-dir BAMS_DIR --template params/template.json \\
        --output-dir CONFIG_DIR

Prints the written params JSON paths (one per line); exits 1 when a run of
the difference could not get its params JSON.

    msi_detection_pre_p.py --print-run --locker-adivar adivar-testNS.json \\
        --locker-msi msi-sensor.locker.txt

Steps 1-2 only (reads the two lockers, writes nothing): prints the run to
analyse (nothing: none), so that the caller can give the container only its
<runs-dir>/<run>/SampleSheet.csv.

    msi_detection_pre_p.py --print-experiment --runs-dir RUNS \\
        --locker-adivar adivar-testNS.json --locker-msi msi-sensor.locker.txt

Steps 1-3 without the BAM folder (reads the lockers and the run's sample
sheet, writes nothing): prints the run's Experiment Name (nothing: no run),
so that the caller can give the container only <bams-dir>/<Experiment Name>/.
Exits 1 when the sample sheet is missing or rejected.
"""

import argparse
import csv
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

# status of a run ADIVAR has finished, in its locker
ADIVAR_DONE = 3

# [Header] values a run must have to be analysed
EXPECTED_HEADER = {"Assay": "Juno", "Description": "ADIVaR_Juno_INCa_v3"}

logger = logging.getLogger("MSI-pre-processing")


def read_locker(locker_path: Path) -> dict:
    """{run: {"status": N, ...}} of ADIVAR's locker JSON; {} when it does not exist yet."""
    if not locker_path.exists():
        return {}
    return json.loads(locker_path.read_text())


def list_runs_done(locker_adivar_path: Path) -> list:
    """Runs ADIVAR has finished (status 3 in its locker), sorted."""
    locker = read_locker(locker_adivar_path)
    return sorted(run for run, state in locker.items() if state.get("status") == ADIVAR_DONE)


def list_runs_msi(locker_msi_path: Path) -> list:
    """Runs of the MSI locker (one per line, already analysed), sorted; []
    when it does not exist yet."""
    if not locker_msi_path.exists():
        return []
    return sorted({line.strip() for line in locker_msi_path.read_text().splitlines() if line.strip()})


def lock(locker_msi_path: Path, run: str) -> None:
    """Appends run to the MSI locker (one run name per line, created if
    missing); nothing when it is already there."""
    text = locker_msi_path.read_text() if locker_msi_path.exists() else ""
    if run in {line.strip() for line in text.splitlines()}:
        logger.info(f"{locker_msi_path}: {run} already there")
        return
    with open(locker_msi_path, "a") as fh:
        if text and not text.endswith("\n"):
            fh.write("\n")
        fh.write(f"{run}\n")
    logger.info(f"{locker_msi_path}: {run} added")


def read_samplesheet_header(path: Path) -> dict:
    """{key: value} of the [Header] section of an Illumina SampleSheet.csv."""
    header, section = {}, None
    with open(path, newline="") as fh:
        for cols in csv.reader(fh):
            cols = [c.strip() for c in cols]
            if not any(cols):
                continue
            if cols[0].startswith("["):
                section = cols[0].strip("[]")
            elif section == "Header":
                header[cols[0]] = cols[1] if len(cols) > 1 else ""
    return header


def read_experiment(run_dir: Path) -> str:
    """Checks run_dir's sample sheet header and returns its Experiment Name
    (the name of its BAM folder); None (error logged) when the run cannot be
    analysed."""
    sheet = run_dir / "SampleSheet.csv"
    if not sheet.is_file():
        logger.error(f"{run_dir.name}: no {sheet}")
        return None
    header = read_samplesheet_header(sheet)

    for key, expected in EXPECTED_HEADER.items():
        if header.get(key) != expected:
            logger.error(f"{run_dir.name}: {key} = {header.get(key)!r} in {sheet}, expected {expected!r}")
            return None
    experiment = header.get("Experiment Name")
    if not experiment:
        logger.error(f"{run_dir.name}: no Experiment Name in {sheet}")
        return None
    return experiment


def make_params(run_dir: Path, bams_dir: Path, template: dict, output_dir: Path) -> Path:
    """Checks run_dir's sample sheet header and BAM folder, writes its params
    JSON and returns its path; None (error logged) when the run cannot be
    analysed."""
    experiment = read_experiment(run_dir)
    if not experiment:
        return None
    sheet = run_dir / "SampleSheet.csv"

    bam_dir = bams_dir / experiment
    bam_count = sum(1 for _ in bam_dir.rglob("*.bam")) if bam_dir.is_dir() else 0
    if not bam_count:
        logger.error(f"{run_dir.name}: no BAM in {bam_dir} (Experiment Name of {sheet})")
        return None
    logger.info(f"{run_dir.name}: {experiment}, {bam_count} BAM(s) in {bam_dir}")

    params = json.loads(json.dumps(template))
    params["msi_detection.inputs"]["run"] = str(bam_dir.resolve())
    date = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    output = output_dir / f"{date}-{run_dir.name}-{experiment}.json"
    output.write_text(json.dumps(params, indent=4) + "\n")
    return output


def main():
    parser = argparse.ArgumentParser(description="Params JSON of the runs finished by ADIVAR and not yet in the MSI locker.")
    parser.add_argument("--print-run", action="store_true",
                        help="Only print the run to analyse (steps 1-2): the other options but the lockers are not needed.")
    parser.add_argument("--print-experiment", action="store_true",
                        help="Only print the run's Experiment Name (its BAM folder): needs --runs-dir and the lockers.")
    parser.add_argument("--runs-dir", help="Folder holding one sub-folder per run (its SampleSheet.csv).")
    parser.add_argument("--locker-adivar", required=True,
                        help="ADIVAR's locker ({run: {\"status\": N}}); status 3 = finished by ADIVAR. Read only.")
    parser.add_argument("--locker-msi", required=True,
                        help="MSI locker: text file, one run name per line, the runs already taken (skipped); the run is appended to it.")
    parser.add_argument("--bams-dir", help="Holds one BAM folder per Experiment Name.")
    parser.add_argument("--template", help="Params JSON whose run input is replaced per run.")
    parser.add_argument("--output-dir", help="Where <date>-<run>-<Experiment Name>.json is written.")
    args = parser.parse_args()
    if args.print_run and args.print_experiment:
        parser.error("--print-run and --print-experiment are exclusive")
    needed = [] if args.print_run else ["runs_dir"] if args.print_experiment else ["runs_dir", "bams_dir", "template", "output_dir"]
    missing = [opt for opt in needed if not getattr(args, opt)]
    if missing:
        parser.error(f"required here: {', '.join('--' + m.replace('_', '-') for m in missing)}")

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")

    # 1-2. Runs finished by ADIVAR and not in the MSI locker
    taken = set(list_runs_msi(Path(args.locker_msi)))
    todo = [run for run in list_runs_done(Path(args.locker_adivar)) if run not in taken]
    if not todo:
        logger.info(f"No run to analyse (check {args.locker_msi})")
        return
    if len(todo) > 1:
        logger.error(f"{len(todo)} runs to analyse, one expected (nothing written): {', '.join(todo)}")
        sys.exit(1)
    logger.info(f"Run to analyse: {todo[0]}")
    if args.print_run:
        print(todo[0])
        return
    if args.print_experiment:
        experiment = read_experiment(Path(args.runs_dir) / todo[0])
        if not experiment:
            sys.exit(1)
        print(experiment)
        return

    # 3-5. Each run's sample sheet header and BAM folder -> its params JSON,
    # then the run into the MSI locker
    runs_dir, output_dir = Path(args.runs_dir), Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    template = json.loads(Path(args.template).read_text())
    failed = []
    for run in todo:
        output = make_params(runs_dir / run, Path(args.bams_dir), template, output_dir)
        if output:
            lock(Path(args.locker_msi), run)
            print(output)
        else:
            failed.append(run)
    if failed:
        logger.error(f"No params JSON for {', '.join(failed)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
