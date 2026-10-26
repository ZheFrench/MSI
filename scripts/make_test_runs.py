#!/usr/bin/env python3
"""
make_test_runs.py

Builds sequencer-like test runs for `BASH.sprocket.sh pending` from the
BAM folders of input/colon-msi (no real sample sheets yet): each run folder
only holds a SampleSheet.csv, the BAMs stay where they are.

    <runs-dir>/<run>/SampleSheet.csv - Experiment Name = source run (its BAM
        folder under input/colon-msi), one [Data] line per BAM: Sample_Name =
        BAM name, Manifest A = alignments_libA, B = alignments_libB, plus
        the controls of a real sample sheet (TN01-03, TPHD730-731, no BAM)

and records each run in the ADIVAR locker: status = last-status = 3
(finished, so listed by pending) or 2 for the last one (still running, so
never listed: checks the filter). Re-runnable: existing test runs are rebuilt.

Usage:
    python3 make_test_runs.py [--source DIR] [--runs-dir DIR] [--adivar-locker JSON]
"""

import argparse
import json
import shutil
from pathlib import Path

CHU = Path("/media/jp/SSD2/CHU")

# (test run name, source run under --source, ADIVAR status)
TEST_RUNS = [
    ("261001_NB501631_0541_ATESTR141", "Run141_Juno", 3),
    ("261002_NB501631_0542_ATESTR144", "Run144_Juno", 3),
    ("261003_NB501631_0543_ATESTR157", "Run157_Juno", 3),
    ("261004_NB501631_0544_ATESTR158", "Run158_Juno", 2),
]

SAMPLESHEET_HEADER = """[Header]
IEMFileVersion,4
Experiment Name,{experiment}
Date,{date}
Investigator Name,TEST
Workflow,GenerateFASTQ
Application,NextSeq FASTQ Only
Assay,Juno
Description,ADIVaR_Juno_INCa_v3
Chemistry,Default

[Manifests]
A,libA.txt
B,libB.txt

[Reads]
151
151

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description,Manifest
"""

# Control lines (negative TN0x, positive TPHD73x) of a real sample sheet
# (260930_NB501631_0540_AHKTHMAFXC), listed before each manifest's samples.
CONTROL_LINES = {
    "A": """TN01_A,TN01,,A01,TSP0049,GTATTGCTGC,,,,ADIVaR_JunoINCaV3,A
TN02_A,TN02,,G02,TSP0063,TTGGTAGATG,,,,ADIVaR_JunoINCaV3,A
TN03_A,TN03,,D03,TSP0068,GTGTGACTAA,,,,ADIVaR_JunoINCaV3,A
TPHD730_A,TPHD730,,G03,TSP0071,CATGCGGAGA,,,,ADIVaR_JunoINCaV3,A
TPHD731_A,TPHD731,,H03,TSP0072,CACTCGATAA,,,,ADIVaR_JunoINCaV3,A
""",
    "B": """TN01_B,TN01,,A04,TSP0073,GAGTAGCCTG,,,,ADIVaR_JunoINCaV3,B
TN02_B,TN02,,G05,TSP0087,GCTTCACTAG,,,,ADIVaR_JunoINCaV3,B
TN03_B,TN03,,D06,TSP0092,CTAGGCTACA,,,,ADIVaR_JunoINCaV3,B
TPHD730_B,TPHD730,,G06,TSP0095,GAACACTCTG,,,,ADIVaR_JunoINCaV3,B
TPHD731_B,TPHD731,,H06,TSP0096,TCACACTTAG,,,,ADIVaR_JunoINCaV3,B
""",
}


def build_run(run_dir: Path, source_dir: Path) -> int:
    """SampleSheet.csv of one test run (one line per BAM of source_dir);
    returns the BAM count."""
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    lines = []
    count = 0
    for manifest in CONTROL_LINES:
        lines.append(CONTROL_LINES[manifest])
        for bam in sorted(source_dir.glob(f"*lib{manifest}/*.bam")):
            lines.append(f"{bam.stem}_{manifest},{bam.stem},,,,,,,,ADIVaR_JunoINCaV3,{manifest}\n")
            count += 1
    header = SAMPLESHEET_HEADER.format(experiment=source_dir.name, date=run_dir.name.split("_")[0])
    (run_dir / "SampleSheet.csv").write_text(header + "".join(lines))
    return count


def main():
    parser = argparse.ArgumentParser(description="Test runs for BASH.sprocket.sh pending.")
    parser.add_argument("--source", default=CHU / "input/colon-msi", type=Path)
    parser.add_argument("--runs-dir", default=CHU / "input/runs", type=Path)
    parser.add_argument("--adivar-locker", default=CHU / "config/adivar-testNS.json", type=Path)
    args = parser.parse_args()

    locker = json.loads(args.adivar_locker.read_text())
    for run, source, status in TEST_RUNS:
        count = build_run(args.runs_dir / run, args.source / source)
        locker[run] = {"status": status, "last-status": status}
        print(f"{args.runs_dir / run}: {count} BAM(s) from {source}, ADIVAR status {status}")
    args.adivar_locker.write_text(json.dumps(locker, indent=4) + "\n")


if __name__ == "__main__":
    main()
