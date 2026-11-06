#!/bin/sh

set -eu

PIPELINE_DIR=/media/jp/SSD2/CHU/code/pipeline-jp
# Apptainer images (*.sif; their recipes: $PIPELINE_DIR/apptainer-recipes/*.def)
CONTAINERS_DIR=/media/jp/SSD2/CHU/containers
# References (genome, panel and ROI BEDs of params/template.json), bound
# read-only in the task containers with $PIPELINE_DIR (python helpers): across
# filesystems, Cromwell soft-links its inputs, whose targets must be visible.
REFDATA_DIR=/media/jp/SSD2/CHU/refData
SCRATCH_DIR=/media/jp/SSD2/CHU/Scratch/MSI/Cron

mkdir -p "$SCRATCH_DIR"
cd "$SCRATCH_DIR"

WDL=$PIPELINE_DIR/msi_detection-1.0.0.wdl
# Cromwell (WDL 1.1) needs Java 17 or later.
JAVA=/usr/lib/jvm/java-17-openjdk-amd64/bin/java
CROMWELL_JAR=$PIPELINE_DIR/soft/cromwell/cromwell-92.jar
CROMWELL_CONFIG=$PIPELINE_DIR/backends.conf/cromwell.conf
# Cromwell runs from here: its executions (cromwell-executions/), its call
# cache database (cromwell-db/), each run's merged inputs, metadata and log.
CROMWELL_DIR=$SCRATCH_DIR/Cromwell
# Images and python helpers (paths of this machine), before each run's params
# JSON: the last file wins, so a params JSON can still override them.
INSTALL=$PIPELINE_DIR/params/install.json
TEMPLATE=$PIPELINE_DIR/params/template.json

# Final reports folder: index copies each run's runs/<run_id>/ into it and
# rebuilds the dashboard. Everything else is computed in $SCRATCH_DIR:
# Cromwell's runs (Cromwell/), the params JSON (Config/) and the copies of
# ADIVAR's locker and of the run's sample sheet (Temp/).
RESULTS_DIR=/media/jp/SSD2/CHU/test/Cron-Results
# Parent folders: $RUNS_DIR/<run>/bcl2fastq/SampleSheet.csv (<run>: name in ADIVAR's
# locker), $BAMS_DIR/<Experiment Name>/ (the run's BAMs, the workflow's input).
RUNS_DIR=/media/jp/SSD2/CHU/input/fake-data/runs
BAMS_DIR=/media/jp/SSD2/CHU/input/fake-data/fake_bam
LOCKER_ADIVAR=/media/jp/SSD2/CHU/config/locker/adivar-testNS.json
LOCKER_MSI=/media/jp/SSD2/CHU/config/locker/msi-sensor.locker.txt

mkdir -p "$RESULTS_DIR"
# Bound into the containers, so it must exist (pre_process appends to it).
touch "$LOCKER_MSI"

# The jobs below read copies, made here in $TEMP_DIR, of ADIVAR's locker
# ($ADIVAR_COPY, once: the three jobs see the same state even if ADIVAR
# updates it meanwhile) and of the run's sample sheet ($SAMPLE_SHEET, below):
# the containers never see ADIVAR's file nor the sequencer's run folder.
TEMP_DIR=$SCRATCH_DIR/Temp
mkdir -p "$TEMP_DIR"
ADIVAR_COPY=$TEMP_DIR/$(basename "$LOCKER_ADIVAR")
cp "$LOCKER_ADIVAR" "$ADIVAR_COPY"

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

# python3 of the python image (none needed on the host), run here (no Slurm
# job: a few lines of JSON), seeing only the binds $1.
py() {
	binds=$1; shift
	apptainer exec --cleanenv --bind "$binds" "$CONTAINERS_DIR/python-MSI-3.12.sif" python3 "$@"
}

# post_p subcommand $2 (job msi-<subcommand>), with the binds $1.
post_process() {
	binds=$1; shift
	python_slurm_job "$1" "$binds" msi_detection_post_p.py "$@"
}

# Run $1 (name of its params JSON $2) with Cromwell, from $CROMWELL_DIR: Cromwell
# takes a single inputs JSON, so $INSTALL and $2 are merged into inputs/$1.json
# (the last file wins). Its log goes to logs/$1.log, its outputs (paths in
# cromwell-executions/) to metadata/$1.json. $3: data folders bound in the
# task containers (APPTAINER_BIND: the run directory, read by path, with
# $PIPELINE_DIR and $REFDATA_DIR).
cromwell_run() {
	mkdir -p "$CROMWELL_DIR/inputs" "$CROMWELL_DIR/metadata" "$CROMWELL_DIR/logs"
	py "$CROMWELL_DIR/inputs,$(ro "$INSTALL"),$(ro "$2")" -c 'import json, sys
merged = {}
for path in sys.argv[2:]:
    with open(path) as f:
        merged.update(json.load(f))
with open(sys.argv[1], "w") as f:
    json.dump(merged, f, indent=4)' "$CROMWELL_DIR/inputs/$1.json" "$INSTALL" "$2"
	(cd "$CROMWELL_DIR" && MSI_CONTAINERS_DIR=$CONTAINERS_DIR APPTAINER_BIND=$3 \
		"$JAVA" -Dconfig.file="$CROMWELL_CONFIG" -jar "$CROMWELL_JAR" run "$WDL" \
		-i "inputs/$1.json" -m "metadata/$1.json" > "logs/$1.log" 2>&1) \
		|| { echo "Cromwell failed on $1: see $CROMWELL_DIR/logs/$1.log" >&2; return 1; }
}

# After run $1: index adds its viz/runs/<run_id>/ (from its viz.tar.gz) to
# $RESULTS_DIR (copied) and rebuilds the dashboard (index_msi_*.js; index_msi_*.html
# written once).
index() {
	post_process "$(ro "$CROMWELL_DIR"),$RESULTS_DIR" index --viz "$RESULTS_DIR" \
		--run-outputs "$CROMWELL_DIR/metadata/$1.json"
}

# Params JSON of the run left to analyse (one at most), from $TEMPLATE:
# pre_process writes $SCRATCH_DIR/Config/<date>-<run>-<Experiment Name>.json (Experiment
# Name: its BAM folder, e.g. Run141_Juno) and prints its path (nothing: no run
# to analyse). It also adds the run to the MSI locker before it runs: a later
# launch (cron) never takes it again, nor is a failed run retried (remove its
# line to rerun it). pre_process fails (see its ERROR lines above) when the run
# has no Experiment Name or BAM, or more than one run is left.
# The run is first read from the lockers (job msi-pick_run, --print-run), then
# its Experiment Name from its SampleSheet.csv (job msi-experiment,
# --print-experiment), so that pre_process gets only that SampleSheet.csv and
# that run's BAM folder, not the whole $RUNS_DIR and $BAMS_DIR.
RUN_ID=$(python_slurm_job pick_run "$(ro "$ADIVAR_COPY"),$(ro "$LOCKER_MSI")" \
	msi_detection_pre_p.py --print-run \
	--locker-adivar "$ADIVAR_COPY" --locker-msi "$LOCKER_MSI") \
	|| { echo "msi_detection_pre_p.py --print-run failed" >&2; exit 1; }
# No run to analyse (logged why): stop here, not an error.
[ -n "$RUN_ID" ] || exit 0
RUN_SHEET=$RUNS_DIR/$RUN_ID/bcl2fastq/SampleSheet.csv
[ -f "$RUN_SHEET" ] || { echo "$RUN_ID: no $RUN_SHEET" >&2; exit 1; }

# Its copy, given to the jobs below as --samplesheet. Copied again at each
# launch (a run retried takes the current one).
SAMPLE_SHEET=$TEMP_DIR/$RUN_ID-SampleSheet.csv
cp "$RUN_SHEET" "$SAMPLE_SHEET"

# Fails (ERROR lines above) when its sample sheet has no Experiment Name: the
# run stays out of the MSI locker, retried at the next launch.
EXPERIMENT=$(python_slurm_job experiment "$(ro "$ADIVAR_COPY"),$(ro "$LOCKER_MSI"),$(ro "$SAMPLE_SHEET")" \
	msi_detection_pre_p.py --print-experiment --samplesheet "$SAMPLE_SHEET" \
	--locker-adivar "$ADIVAR_COPY" --locker-msi "$LOCKER_MSI") \
	|| { echo "msi_detection_pre_p.py --print-experiment failed" >&2; exit 1; }
RUN_BAMS=$BAMS_DIR/$EXPERIMENT
[ -d "$RUN_BAMS" ] || { echo "$RUN_ID: no $RUN_BAMS (Experiment Name of $RUN_SHEET)" >&2; exit 1; }

mkdir -p "$SCRATCH_DIR/Config"
CONFIG=$(python_slurm_job pre_process "$(ro "$SAMPLE_SHEET"),$(ro "$RUN_BAMS"),$(ro "$ADIVAR_COPY"),$LOCKER_MSI,$SCRATCH_DIR/Config" \
	msi_detection_pre_p.py --samplesheet "$SAMPLE_SHEET" \
	--locker-adivar "$ADIVAR_COPY" --locker-msi "$LOCKER_MSI" \
	--bams-dir "$BAMS_DIR" --template "$TEMPLATE" --output-dir "$SCRATCH_DIR/Config") \
	|| { echo "msi_detection_pre_p.py failed" >&2; exit 1; }
# No run to analyse (pre_process logged why): stop here, not an error.
[ -n "$CONFIG" ] || exit 0

# Named after its params JSON (<date>-<run>-<Experiment Name>, one per run).
RUN=$(basename "$CONFIG" .json)
echo "=== $CONFIG"
cromwell_run "$RUN" "$CONFIG" "$(ro "$RUN_BAMS"),$(ro "$PIPELINE_DIR"),$(ro "$REFDATA_DIR")"
index "$RUN"
# Only reached when every step above succeeded (set -e).
echo "$(date '+%Y-%m-%d %H:%M:%S') [DONE]: $RUN_ID ($EXPERIMENT) analysed, reports in $RESULTS_DIR/index_msi_run.html"
