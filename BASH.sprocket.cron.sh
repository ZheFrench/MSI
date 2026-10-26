#!/bin/sh


set -eu
# cron's PATH is minimal: sprocket is in ~/.cargo/bin.
#export PATH=/home/jp/.cargo/bin:$PATH

WDL_DIR=/media/jp/SSD2/CHU/code/pipeline-jp/wdl

WDL=$WDL_DIR/msi_detection-0.1.0.wdl
SPROCKET_CONFIG=$WDL_DIR/backends.conf/sprocket.prod.toml
TEMPLATE=$WDL_DIR/params/template.json

VIZ=/media/jp/SSD2/CHU/test/viz
WORK=/media/jp/SSD2/CHU/test/work
RUNS_DIR=/media/jp/SSD2/CHU/input/fake-data/runs
BAM_ROOT=/media/jp/SSD2/CHU/input/fake-data/fake_bam
LOCKER_ADIVAR=/media/jp/SSD2/CHU/config/adivar-testNS.json
LOCKER_MSI=/media/jp/SSD2/CHU/test/msi-sensor.txt

mkdir -p "$VIZ" "$WORK"
# Bound into the containers, so it must exist (pre_process appends to it).
touch "$LOCKER_MSI"

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

# After run $1 (sprocket --index-on name): index adds its viz/runs/<run_id>/ to
# $VIZ (copied) and rebuilds the dashboard (index_msi_*.js; index_msi_*.html
# written once). Not _latest: a concurrent sprocket run would repoint it.
index() {
	post_process "$(ro "$WORK/sprocket"),$VIZ" index --viz "$VIZ" \
		--run-outputs "$(realpath "$WORK/sprocket/index/$1/outputs.json")"
}

# Params JSON of the run left to analyse (one at most), from $TEMPLATE:
# pre_process writes $WORK/config/<date>-<run>-<Experiment Name>.json (Experiment
# Name: its BAM folder, e.g. Run141_Juno) and prints its path (nothing: no run
# to analyse). It also adds the run to the MSI locker before it runs: a later
# launch (cron) never takes it again, nor is a failed run retried (remove its
# line to rerun it). pre_process fails (see its ERROR lines above) when the run is
# rejected or more than one run is left.
mkdir -p "$WORK/config"
CONFIG=$(python_slurm_job pre_process "$(ro "$RUNS_DIR"),$(ro "$BAM_ROOT"),$(ro "$LOCKER_ADIVAR"),$LOCKER_MSI,$WORK/config" \
	msi_detection_pre_p.py --runs-dir "$RUNS_DIR" \
	--locker-adivar "$LOCKER_ADIVAR" --locker-msi "$LOCKER_MSI" \
	--bam-root "$BAM_ROOT" --template "$TEMPLATE" --output-dir "$WORK/config") \
	|| { echo "msi_detection_pre_p.py failed" >&2; exit 1; }
# No run to analyse (pre_process logged why): stop here, not an error.
[ -n "$CONFIG" ] || exit 0

# Named after its params JSON (<date>-<run>-<Experiment Name>, one per run).
RUN=msi_detection/$(basename "$CONFIG" .json)
echo "=== $CONFIG"
# Its outputs JSON (stdout) stays out of the cron log: it is the run's outputs.json.
sprocket run -v -s -c "$SPROCKET_CONFIG" -o "$WORK/sprocket" --index-on "$RUN" "$WDL" "@$CONFIG" >/dev/null
index "$RUN"
