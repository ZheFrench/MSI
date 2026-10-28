#!/bin/sh

set -eu

WDL_DIR=/media/jp/SSD2/CHU/code/pipeline-jp
SCRATCH_DIR=/media/jp/SSD2/CHU/code/pipeline-jp

cd "$SCRATCH_DIR"

WDL=$WDL_DIR/msi_detection-1.0.0.wdl
SPROCKET_CONFIG=$WDL_DIR/backends.conf/sprocket.cron.toml
# Images and python helpers (paths of this machine), before each run's params
# JSON: the last file wins, so a params JSON can still override them.
INSTALL=$WDL_DIR/params/install.json
TEMPLATE=$WDL_DIR/params/template.json

VIZ=/media/jp/SSD2/CHU/test/viz
WORK=/media/jp/SSD2/CHU/test/work
# Parent folders: $RUNS_DIR/<run>/SampleSheet.csv (<run>: name in ADIVAR's
# locker), $BAMS_DIR/<Experiment Name>/ (the run's BAMs, the workflow's input).
RUNS_DIR=/media/jp/SSD2/CHU/input/fake-data/runs
BAMS_DIR=/media/jp/SSD2/CHU/input/fake-data/fake_bam
LOCKER_ADIVAR=/media/jp/SSD2/CHU/config/locker/adivar-testNS.json
LOCKER_MSI=/media/jp/SSD2/CHU/config/locker/msi-sensor.locker.txt

mkdir -p "$VIZ" "$WORK"
# Bound into the containers, so it must exist (pre_process appends to it).
touch "$LOCKER_MSI"

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

# Params JSON of the run left to analyse (one at most), from $TEMPLATE:
# pre_process writes $WORK/config/<date>-<run>-<Experiment Name>.json (Experiment
# Name: its BAM folder, e.g. Run141_Juno) and prints its path (nothing: no run
# to analyse). It also adds the run to the MSI locker before it runs: a later
# launch (cron) never takes it again, nor is a failed run retried (remove its
# line to rerun it). pre_process fails (see its ERROR lines above) when the run
# has no Experiment Name or BAM, or more than one run is left.
# The run is first read from the lockers (job msi-pick_run, --print-run), then
# its Experiment Name from its SampleSheet.csv (job msi-experiment,
# --print-experiment), so that pre_process gets only that SampleSheet.csv and
# that run's BAM folder, not the whole $RUNS_DIR and $BAMS_DIR.
RUN_ID=$(python_slurm_job pick_run "$(ro "$LOCKER_ADIVAR"),$(ro "$LOCKER_MSI")" \
	msi_detection_pre_p.py --print-run \
	--locker-adivar "$LOCKER_ADIVAR" --locker-msi "$LOCKER_MSI") \
	|| { echo "msi_detection_pre_p.py --print-run failed" >&2; exit 1; }
# No run to analyse (logged why): stop here, not an error.
[ -n "$RUN_ID" ] || exit 0
SAMPLE_SHEET=$RUNS_DIR/$RUN_ID/SampleSheet.csv
[ -f "$SAMPLE_SHEET" ] || { echo "$RUN_ID: no $SAMPLE_SHEET" >&2; exit 1; }

# Fails (ERROR lines above) when its sample sheet has no Experiment Name: the
# run stays out of the MSI locker, retried at the next launch.
EXPERIMENT=$(python_slurm_job experiment "$(ro "$LOCKER_ADIVAR"),$(ro "$LOCKER_MSI"),$(ro "$SAMPLE_SHEET")" \
	msi_detection_pre_p.py --print-experiment --runs-dir "$RUNS_DIR" \
	--locker-adivar "$LOCKER_ADIVAR" --locker-msi "$LOCKER_MSI") \
	|| { echo "msi_detection_pre_p.py --print-experiment failed" >&2; exit 1; }
RUN_BAMS=$BAMS_DIR/$EXPERIMENT
[ -d "$RUN_BAMS" ] || { echo "$RUN_ID: no $RUN_BAMS (Experiment Name of $SAMPLE_SHEET)" >&2; exit 1; }

mkdir -p "$WORK/config"
CONFIG=$(python_slurm_job pre_process "$(ro "$SAMPLE_SHEET"),$(ro "$RUN_BAMS"),$(ro "$LOCKER_ADIVAR"),$LOCKER_MSI,$WORK/config" \
	msi_detection_pre_p.py --runs-dir "$RUNS_DIR" \
	--locker-adivar "$LOCKER_ADIVAR" --locker-msi "$LOCKER_MSI" \
	--bams-dir "$BAMS_DIR" --template "$TEMPLATE" --output-dir "$WORK/config") \
	|| { echo "msi_detection_pre_p.py failed" >&2; exit 1; }
# No run to analyse (pre_process logged why): stop here, not an error.
[ -n "$CONFIG" ] || exit 0

# Named after its params JSON (<date>-<run>-<Experiment Name>, one per run).
RUN=msi_detection/$(basename "$CONFIG" .json)
echo "=== $CONFIG"
# Its outputs JSON (stdout) stays out of the cron log: it is the run's outputs.json.
sprocket run -v -s -c "$SPROCKET_CONFIG" -o "$WORK/sprocket" --index-on "$RUN" "$WDL" "@$INSTALL" "@$CONFIG" >/dev/null
index "$RUN"
