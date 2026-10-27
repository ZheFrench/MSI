#!/bin/sh
# Batch mode: runs the workflow on every params JSON of a folder, one after
# the other (no pre_process, no MSI locker). Cron mode: BASH.sprocket.cron.sh.
#   ./companion-scripts/BASH.sprocket.batch.sh CONFIG_DIR VIZ WORK
# VIZ (reports and dashboard) and WORK (sprocket runs), required, as in cron
# mode; one VIZ per group keeps their dashboards apart.
# CONFIG_DIR from BASH.make_configs.sh; results as TSV after: BASH.export_tsv.sh VIZ.
set -eu

WDL_DIR=/media/jp/SSD2/CHU/code/pipeline-jp

WDL=$WDL_DIR/msi_detection-0.1.0.wdl
SPROCKET_CONFIG=$WDL_DIR/backends.conf/sprocket.batch.toml
# Images and python helpers (paths of this machine), before each run's params
# JSON: the last file wins, so a params JSON can still override them.
INSTALL=$WDL_DIR/params/install.json

# Folder of the params JSON to run.
CONFIGS=${1:?usage: $0 CONFIG_DIR VIZ WORK}
# Absolute (apptainer binds); -m: they may not exist yet.
VIZ=$(realpath -m "${2:?usage: $0 CONFIG_DIR VIZ WORK}")
WORK=$(realpath -m "${3:?usage: $0 CONFIG_DIR VIZ WORK}")

mkdir -p "$VIZ" "$WORK"

# Apptainer bind of a path the container only reads.
ro() {
	echo "$1:$1:ro"
}

# Script $3 of main-scripts/ as a Slurm job (srun, job msi-$1) in the python image,
# which sees $WDL_DIR (read-only) and the binds $2.
python_slurm_job() {
	name=$1; binds=$2; script=$3; shift 3
	srun --partition=prod --ntasks=1 --cpus-per-task=1 --mem=1G --job-name="msi-$name" \
		apptainer exec --cleanenv --bind "$(ro "$WDL_DIR"),$binds" \
		"$WDL_DIR/apptainer-recipes/python-3.12.sif" \
		python3 "$WDL_DIR/main-scripts/$script" "$@"
}

# post_p subcommand $2 (job msi-<subcommand>), with the binds $1.
post_process() {
	binds=$1; shift
	python_slurm_job "$1" "$binds" msi_detection_post_p.py "$@"
}

# After run $1 (sprocket --index-on name): index adds its viz/runs/<run_id>/ to
# $VIZ (copied) and rebuilds the dashboard (index_msi_*.js; index_msi_*.html
# written once). Not _latest: a concurrent sprocket run would repoint it.
index() {
	post_process "$(ro "$WORK/sprocket"),$VIZ" index --viz "$VIZ" \
		--run-outputs "$(realpath "$WORK/sprocket/index/$1/outputs.json")"
}

FAILED=""
# Every *.json of the folder, in name order.
LIST=$(find "$CONFIGS" -maxdepth 1 -name '*.json' -type f | sort)
[ -n "$LIST" ] || { echo "no *.json in $CONFIGS" >&2; exit 1; }

for p in $LIST; do
	echo "=== $p"
	# Named after its params JSON: $WORK/sprocket/index/msi_detection/<name>/
	# (links to its outputs and its outputs.json).
	RUN=msi_detection/$(basename "$p" .json)
	if sprocket run -v -s -c "$SPROCKET_CONFIG" -o "$WORK/sprocket" --index-on "$RUN" "$WDL" "@$INSTALL" "@$p" && index "$RUN"; then
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

