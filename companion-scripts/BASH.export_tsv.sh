#!/bin/sh
# export_tsv.py as a Slurm job (srun, job msi-export_tsv) in the python image:
# the results of a reports folder as a TSV, one line per sample.
#   ./companion-scripts/BASH.export_tsv.sh VIZ_DIR OUTPUT_DIR
# VIZ_DIR: RESULTS_DIR of BASH.msi_detection.batch.sh or BASH.msi_detection.cron.sh (reports folder),
# read only; the TSV is written in OUTPUT_DIR, OUTPUT_DIR/results_msi_<date>.tsv.
set -eu

PIPELINE_DIR=/media/jp/SSD2/CHU/code/pipeline-jp
# Apptainer images (*.sif; their recipes: $PIPELINE_DIR/apptainer-recipes/*.def)
CONTAINERS_DIR=/media/jp/SSD2/CHU/containers
# Scratch: the script runs from it (once its arguments are resolved)
SCRATCH_DIR=/media/jp/SSD2/CHU/Scratch/MSI

VIZ_DIR=$(realpath "${1:?usage: $0 VIZ_DIR OUTPUT_DIR}")
OUTPUT_DIR=$(realpath "${2:?usage: $0 VIZ_DIR OUTPUT_DIR}")
mkdir -p "$SCRATCH_DIR"
cd "$SCRATCH_DIR"

# Apptainer bind of a path the container only reads.
ro() {
	echo "$1:$1:ro"
}

# OUTPUT_DIR (read-write) first: the read-only binds after it stay read-only
# when they are inside it.
srun --partition=prod --ntasks=1 --cpus-per-task=1 --mem=1G --job-name=msi-export_tsv \
	apptainer exec --cleanenv --bind "$OUTPUT_DIR,$(ro "$PIPELINE_DIR"),$(ro "$VIZ_DIR")" \
	"$CONTAINERS_DIR/python-MSI-3.12.sif" \
	python3 -I "$PIPELINE_DIR/companion-scripts/export_tsv.py" --viz "$VIZ_DIR" --output-dir "$OUTPUT_DIR"
