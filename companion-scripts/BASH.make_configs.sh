#!/bin/sh
# make_configs.py as a Slurm job (srun, job msi-make_configs) in the python
# image: the params JSON of the runs of RUNS_FILE, for BASH.msi_detection.batch.sh.
# RUNS_FILE: one line per run (every sample) or per run + sample ID (only those
# samples), columns separated by tabs or spaces.
#   ./companion-scripts/BASH.make_configs.sh BAMS_DIR RUNS_FILE OUTPUT_DIR
# BAMS_DIR: group folder of BAM folders (input/colon-msi...); OUTPUT_DIR:
# where <run>_params.json is written (created if needed). Exits 1 when a run
# has no BAM (listed in OUTPUT_DIR/failed-runs.txt).
set -eu

PIPELINE_DIR=/media/jp/SSD2/CHU/code/pipeline-jp
# Apptainer images (*.sif; their recipes: $PIPELINE_DIR/apptainer-recipes/*.def)
CONTAINERS_DIR=/media/jp/SSD2/CHU/containers
# Scratch: the script runs from it (once its arguments are resolved)
SCRATCH_DIR=/media/jp/SSD2/CHU/Scratch/MSI

BAMS_DIR=$(realpath "${1:?usage: $0 BAMS_DIR RUNS_FILE OUTPUT_DIR}")
RUNS_FILE=$(realpath "${2:?usage: $0 BAMS_DIR RUNS_FILE OUTPUT_DIR}")
# -m: OUTPUT_DIR may not exist yet.
OUTPUT_DIR=$(realpath -m "${3:?usage: $0 BAMS_DIR RUNS_FILE OUTPUT_DIR}")
mkdir -p "$SCRATCH_DIR"
cd "$SCRATCH_DIR"

# Bound by apptainer: must exist before.
mkdir -p "$OUTPUT_DIR"

# Apptainer bind of a path the container only reads.
ro() {
	echo "$1:$1:ro"
}

srun --partition=prod --ntasks=1 --cpus-per-task=1 --mem=1G --job-name=msi-make_configs \
	apptainer exec --cleanenv --bind "$(ro "$PIPELINE_DIR"),$(ro "$BAMS_DIR"),$(ro "$RUNS_FILE"),$OUTPUT_DIR" \
	"$CONTAINERS_DIR/python-MSI-3.12.sif" \
	python3 -I "$PIPELINE_DIR/companion-scripts/make_configs.py" "$BAMS_DIR" \
	--runs-file "$RUNS_FILE" --output-dir "$OUTPUT_DIR"
