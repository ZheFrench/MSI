#!/bin/sh
# Batch mode: runs the workflow on every params JSON of a folder, one after
# the other (no pre_process, no MSI locker). Cron mode: BASH.msi_detection.cron.sh.
#   ./companion-scripts/BASH.msi_detection.batch.sh CONFIGS_DIR RESULTS_DIR
# RESULTS_DIR: final reports folder (index copies each run's reports into it and
# rebuilds the dashboard), required; one RESULTS_DIR per group keeps their dashboards
# apart. sprocket's runs are computed in the scratch, $SCRATCH_DIR/Sprocket/,
# as in cron mode.
# CONFIGS_DIR from BASH.make_configs.sh; results as TSV after: BASH.export_tsv.sh RESULTS_DIR.
set -eu

PIPELINE_DIR=/media/jp/SSD2/CHU/code/pipeline-jp
# Apptainer images (*.sif; their recipes: $PIPELINE_DIR/apptainer-recipes/*.def)
CONTAINERS_DIR=/media/jp/SSD2/CHU/containers
# Scratch: the script runs from it, sprocket's runs are in it
SCRATCH_DIR=/media/jp/SSD2/CHU/Scratch/MSI/Batch-1

WDL=$PIPELINE_DIR/msi_detection-1.0.0.wdl
SPROCKET_CONFIG=$PIPELINE_DIR/backends.conf/sprocket.batch.toml
# Images and python helpers (paths of this machine), before each run's params
# JSON: the last file wins, so a params JSON can still override them.
INSTALL=$PIPELINE_DIR/params/install.json

# Folder of the params JSON to run.
CONFIGS_DIR=$(realpath "${1:?usage: $0 CONFIGS_DIR RESULTS_DIR}")
# Absolute (apptainer binds, and resolved before the cd below); -m: it may
# not exist yet.
RESULTS_DIR=$(realpath -m "${2:?usage: $0 CONFIGS_DIR RESULTS_DIR}")

mkdir -p "$RESULTS_DIR" "$SCRATCH_DIR"
cd "$SCRATCH_DIR"

# Apptainer bind of a path the container only reads.
ro() {
	echo "$1:$1:ro"
}

# Script $3 of main-scripts/ as a Slurm job (srun, job msi-$1) in the python image,
# which sees $PIPELINE_DIR (read-only) and the binds $2.
python_slurm_job() {
	name=$1; binds=$2; script=$3; shift 3
	srun --partition=prod --ntasks=1 --cpus-per-task=1 --mem=1G --job-name="msi-$name" \
		apptainer exec --cleanenv --bind "$(ro "$PIPELINE_DIR"),$binds" \
		"$CONTAINERS_DIR/python-MSI-3.12.sif" \
		python3 "$PIPELINE_DIR/main-scripts/$script" "$@"
}

# post_p subcommand $2 (job msi-<subcommand>), with the binds $1.
post_process() {
	binds=$1; shift
	python_slurm_job "$1" "$binds" msi_detection_post_p.py "$@"
}

# After run $1 (sprocket --index-on name): index adds its viz/runs/<run_id>/ to
# $RESULTS_DIR (copied) and rebuilds the dashboard (index_msi_*.js; index_msi_*.html
# written once). Not _latest: a concurrent sprocket run would repoint it.
index() {
	post_process "$(ro "$SCRATCH_DIR/Sprocket"),$RESULTS_DIR" index --viz "$RESULTS_DIR" \
		--run-outputs "$(realpath "$SCRATCH_DIR/Sprocket/index/$1/outputs.json")"
}

FAILED=""
# Every *.json of the folder, in name order.
LIST=$(find "$CONFIGS_DIR" -maxdepth 1 -name '*.json' -type f | sort)
[ -n "$LIST" ] || { echo "no *.json in $CONFIGS_DIR" >&2; exit 1; }

for p in $LIST; do
	echo "=== $p"
	# Named after its params JSON: $SCRATCH_DIR/Sprocket/index/msi_detection/<name>/
	# (links to its outputs and its outputs.json).
	RUN=msi_detection/$(basename "$p" .json)
	if sprocket run -v -s -c "$SPROCKET_CONFIG" -o "$SCRATCH_DIR/Sprocket" --index-on "$RUN" "$WDL" "@$INSTALL" "@$p" && index "$RUN"; then
		echo "    ok"
	else
		echo "    FAILED"
		FAILED="$FAILED $p"
	fi
done

if [ -n "$FAILED" ]; then
	echo "Failed:"; for p in $FAILED; do echo "  $p"; done
	exit 1
fi
echo "$(date '+%Y-%m-%d %H:%M:%S') OK: every params JSON of $CONFIGS_DIR analysed, reports in $RESULTS_DIR/index_msi_run.html"
