#!/bin/sh
# Batch mode: runs the workflow on every params JSON of a folder, one after
# the other (no pre_process, no MSI locker). Cron mode: BASH.msi_detection.cron.sh.
#   ./companion-scripts/BASH.msi_detection.batch.sh CONFIGS_DIR RESULTS_DIR
# RESULTS_DIR: final reports folder (index copies each run's reports into it and
# rebuilds the dashboard), required; one RESULTS_DIR per group keeps their dashboards
# apart. Cromwell's runs are computed in the scratch, $SCRATCH_DIR/Cromwell/,
# as in cron mode.
# CONFIGS_DIR from BASH.make_configs.sh; results as TSV after: BASH.export_tsv.sh RESULTS_DIR.
set -eu

PIPELINE_DIR=/media/jp/SSD2/CHU/code/pipeline-jp
# Apptainer images (*.sif; their recipes: $PIPELINE_DIR/apptainer-recipes/*.def)
CONTAINERS_DIR=/media/jp/SSD2/CHU/containers
# Read by backends.conf/cromwell.conf (submit) to find each task's .sif.
export MSI_CONTAINERS_DIR=$CONTAINERS_DIR
# References (genome, panel and ROI BEDs of params/template.json), bound
# read-only in the task containers with $PIPELINE_DIR (python helpers): across
# filesystems, Cromwell soft-links its inputs, whose targets must be visible.
REFDATA_DIR=/media/jp/SSD2/CHU/refData
# Scratch: the script runs from it, Cromwell's runs are in it
SCRATCH_DIR=/media/jp/SSD2/CHU/Scratch/MSI/Batch-1

WDL=$PIPELINE_DIR/msi_detection-1.0.0.wdl
# Cromwell (WDL 1.1) needs Java 17 or later.
JAVA=/usr/lib/jvm/java-17-openjdk-amd64/bin/java
CROMWELL_JAR=$PIPELINE_DIR/soft/cromwell/cromwell-92.jar
# Same settings as the cron mode; its call cache is this scratch's own
# (cromwell-db/, relative to $CROMWELL_DIR).
CROMWELL_CONFIG=$PIPELINE_DIR/backends.conf/cromwell.conf
# Cromwell runs from here: its executions (cromwell-executions/), its call
# cache database (cromwell-db/), each run's merged inputs, metadata and log.
CROMWELL_DIR=$SCRATCH_DIR/Cromwell
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
# cromwell-executions/) to metadata/$1.json. The run directory of $2, read by
# path, is bound in the task containers (APPTAINER_BIND), with $PIPELINE_DIR and
# $REFDATA_DIR.
cromwell_run() {
	mkdir -p "$CROMWELL_DIR/inputs" "$CROMWELL_DIR/metadata" "$CROMWELL_DIR/logs"
	py "$CROMWELL_DIR/inputs,$(ro "$INSTALL"),$(ro "$2")" -c 'import json, sys
merged = {}
for path in sys.argv[2:]:
    with open(path) as f:
        merged.update(json.load(f))
with open(sys.argv[1], "w") as f:
    json.dump(merged, f, indent=4)' "$CROMWELL_DIR/inputs/$1.json" "$INSTALL" "$2"
	run_dir=$(py "$(ro "$CROMWELL_DIR/inputs")" -c 'import json, sys
with open(sys.argv[1]) as f:
    print(json.load(f)["msi_detection.inputs"]["run"])' "$CROMWELL_DIR/inputs/$1.json")
	(cd "$CROMWELL_DIR" && APPTAINER_BIND=$(ro "$run_dir"),$(ro "$PIPELINE_DIR"),$(ro "$REFDATA_DIR") \
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

FAILED=""
# Every *.json of the folder, in name order.
LIST=$(find "$CONFIGS_DIR" -maxdepth 1 -name '*.json' -type f | sort)
[ -n "$LIST" ] || { echo "no *.json in $CONFIGS_DIR" >&2; exit 1; }

for p in $LIST; do
	echo "=== $p"
	# Named after its params JSON: $CROMWELL_DIR/{inputs,metadata,logs}/<name>.*
	RUN=$(basename "$p" .json)
	if cromwell_run "$RUN" "$p" && index "$RUN"; then
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
