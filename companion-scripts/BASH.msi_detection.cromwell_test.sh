#!/bin/sh
# Test run: the workflow with Cromwell on one params JSON, by hand (no MSI
# locker, no reports folder). Checks the WDL with womtool, merges install.json
# with the params JSON, runs Cromwell and prints the status and outputs.
#   ./companion-scripts/BASH.msi_detection.cromwell_test.sh [PARAMS_JSON] [WORK_DIR]
# PARAMS_JSON: params JSON of the run [default: Run141_Juno, fake data].
# WORK_DIR: Cromwell's folder (cromwell-executions/, cromwell-db/ call cache,
# inputs.json, metadata.json, cromwell.log) [default: $SCRATCH_DIR/Cromwell-jp];
# kept from one test to the next, so a rerun reuses the call cache.
set -eu

PIPELINE_DIR=/media/jp/SSD2/CHU/code/pipeline-jp
# Apptainer images (*.sif; their recipes: $PIPELINE_DIR/apptainer-recipes/*.def)
CONTAINERS_DIR=/media/jp/SSD2/CHU/containers
# References (genome, panel and ROI BEDs of params/template.json), bound
# read-only in the task containers with $PIPELINE_DIR (python helpers): across
# filesystems, Cromwell soft-links its inputs, whose targets must be visible.
REFDATA_DIR=/media/jp/SSD2/CHU/refData
SCRATCH_DIR=/media/jp/SSD2/CHU/Scratch/MSI

WDL=$PIPELINE_DIR/msi_detection-1.0.0.wdl
INSTALL=$PIPELINE_DIR/params/install.json
# Cromwell (WDL 1.1) needs Java 17 or later.
JAVA=/usr/lib/jvm/java-17-openjdk-amd64/bin/java
CROMWELL_JAR=$PIPELINE_DIR/soft/cromwell/cromwell-92.jar
WOMTOOL_JAR=$PIPELINE_DIR/soft/cromwell/womtool-92.jar
CROMWELL_CONFIG=$PIPELINE_DIR/backends.conf/cromwell.conf

PARAMS=$(realpath "${1:-$SCRATCH_DIR/Cron/Config/2026-10-29_17-24-03-261001_NB501631_0541_ATESTR141-Run141_Juno.json}")
WORK_DIR=$(realpath -m "${2:-$SCRATCH_DIR/Cromwell-jp}")

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# Apptainer bind of a path the container only reads.
ro() {
	echo "$1:$1:ro"
}

# python3 of the python image (none needed on the host), run here (no Slurm
# job: a few lines of JSON), seeing only the binds $1.
py() {
	binds=$1; shift
	apptainer exec --cleanenv --bind "$binds" "$CONTAINERS_DIR/python-MSI-3.12.sif" python3 "$@"
}

echo "=== womtool validate"
"$JAVA" -jar "$WOMTOOL_JAR" validate "$WDL"

# Cromwell takes a single inputs JSON: install.json then the params JSON (the
# last file wins).
py "$WORK_DIR,$(ro "$INSTALL"),$(ro "$PARAMS")" -c 'import json, sys
merged = {}
for path in sys.argv[2:]:
    with open(path) as f:
        merged.update(json.load(f))
with open(sys.argv[1], "w") as f:
    json.dump(merged, f, indent=4)' inputs.json "$INSTALL" "$PARAMS"

# The run directory, read by path by the tasks: bound in their containers.
RUN_DIR=$(py "$WORK_DIR" -c 'import json, sys
with open(sys.argv[1]) as f:
    print(json.load(f)["msi_detection.inputs"]["run"])' inputs.json)
[ -d "$RUN_DIR" ] || { echo "No run directory $RUN_DIR (inputs.run of $PARAMS)" >&2; exit 1; }

echo "=== cromwell run: $PARAMS"
echo "    folder: $WORK_DIR (log of this run: cromwell-executions/workflow-logs/workflow.<id>.log; all of Cromwell's: cromwell.log; jobs: squeue)"
status=0
MSI_CONTAINERS_DIR=$CONTAINERS_DIR APPTAINER_BIND=$(ro "$RUN_DIR"),$(ro "$PIPELINE_DIR"),$(ro "$REFDATA_DIR") \
	"$JAVA" -Dconfig.file="$CROMWELL_CONFIG" -jar "$CROMWELL_JAR" \
	run "$WDL" -i inputs.json -m metadata.json > cromwell.log 2>&1 || status=$?

# Status, failures, call cache hits and outputs, from Cromwell's metadata.
[ -s metadata.json ] || { echo "No metadata.json: see $WORK_DIR/cromwell.log" >&2; exit 1; }
py "$WORK_DIR" -c 'import json
with open("metadata.json") as f:
    m = json.load(f)
print("=== status:", m["status"])
def failures(items, indent="    "):
    for item in items:
        print(indent + item.get("message", ""))
        failures(item.get("causedBy", []), indent + "  ")
failures(m.get("failures", []))
for name, calls in sorted(m.get("calls", {}).items()):
    for call in calls:
        task = name.split(".")[-1]
        shard = call.get("shardIndex", -1)
        if shard >= 0:
            task += "[%d]" % shard
        cached = " (cache)" if call.get("callCaching", {}).get("hit") else ""
        print("    %s: %s%s" % (task, call.get("executionStatus"), cached))
print("=== outputs")
for name, value in m.get("outputs", {}).items():
    print("    %s: %s" % (name, value))'
exit $status
