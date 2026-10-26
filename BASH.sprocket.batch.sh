#!/bin/sh
# Batch mode: runs the workflow on every params JSON of a folder, one after
# the other (no pre_process, no MSI locker). Cron mode: BASH.sprocket.sh.
#   ./BASH.sprocket.batch.sh CONFIG_DIR [OUT_DIR]
# OUT_DIR (default /media/jp/SSD2/CHU/test-batch) gets viz/ (dashboard) and
# work/ (sprocket runs): one OUT_DIR per group keeps their outputs apart.
set -eu

WDL_DIR=/media/jp/SSD2/CHU/code/pipeline-jp/wdl

WDL=$WDL_DIR/msi_detection-0.1.0.wdl
SPROCKET_CONFIG=$WDL_DIR/backends.conf/sprocket.prod.toml

# Folder of the params JSON to run.
CONFIGS=${1:?usage: $0 CONFIG_DIR [OUT_DIR]}
OUT=${2:-/media/jp/SSD2/CHU/test-batch}
VIZ=$OUT/viz
WORK=$OUT/work

mkdir -p "$VIZ" "$WORK"

# Apptainer bind of a path the container only reads.
ro() {
	echo "$1:$1:ro"
}

# Script $3 of scripts/ as a Slurm job (srun, job msi-$1) in the python image,
# which sees $WDL_DIR (read-only) and the binds $2.
python_slurm_job() {
	name=$1; binds=$2; script=$3; shift 3
	srun --partition=prod --ntasks=1 --cpus-per-task=1 --mem=1G --job-name="msi-$name" \
		apptainer exec --bind "$(ro "$WDL_DIR"),$binds" \
		"$WDL_DIR/apptainer-recipes/python-3.12.sif" \
		python3 "$WDL_DIR/scripts/$script" "$@"
}

# post_p subcommand $2 (job msi-<subcommand>), with the binds $1.
post_process() {
	binds=$1; shift
	python_slurm_job "$1" "$binds" msi_detection_post_p.py "$@"
}

# After each run: index adds its viz/runs/<run_id>/ to $VIZ (copied) and
# rebuilds the dashboard (index_msi_*.js; index_msi_*.html written once).
index() {
	post_process "$(ro "$WORK/sprocket"),$VIZ" index --viz "$VIZ" \
		--run-outputs "$(realpath "$WORK/sprocket/runs/msi_detection/_latest")/outputs.json"
}

FAILED=""
# Every *.json of the folder, in name order.
LIST=$(find "$CONFIGS" -maxdepth 1 -name '*.json' -type f | sort)
[ -n "$LIST" ] || { echo "no *.json in $CONFIGS" >&2; exit 1; }

for p in $LIST; do
	echo "=== $p"
	if sprocket run -v -s -c "$SPROCKET_CONFIG" -o "$WORK/sprocket" "$WDL" "@$p" && index; then
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

