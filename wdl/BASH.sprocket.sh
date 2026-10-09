#!/bin/sh
# Runs msi_detection on Slurm (every task in its Apptainer image). Results stay
# in sprocket's run folder, $WORK/sprocket/runs/msi_detection/<timestamp>/
# (paths in its outputs.json).
# Usage: [CONFIG=backends.conf/sprocket.prod.toml] BASH.sprocket.sh [params.json | params_dir]
#   params_dir: every *.json in it (sub-folders included) is run in turn; a
#   failed run is skipped (no index) and listed at the end.
#   defaults: backends.conf/sprocket.local.toml (this workstation),
#             params/example_params_1.json
set -eu
WDL_DIR=$(dirname "$(realpath "$0")")
WDL=$WDL_DIR/msi_detection-0.1.0.wdl
PARAMS=${1:-$WDL_DIR/params/example_params_1.json}
CONFIG=${CONFIG:-$WDL_DIR/backends.conf/sprocket.local.toml}

# 1. Destinations (not in the params JSON): VIZ, the reports folder, filled
#    after each run by index (below); WORK, sprocket's runs ($WORK/sprocket/).
VIZ=/media/jp/SSD2/CHU/test/viz
WORK=/media/jp/SSD2/CHU/test/work


mkdir -p "$VIZ" "$WORK" 

# After each run: index adds its viz/runs/<run_id>/ to $VIZ (moved) and
# rebuilds the dashboard (index_msi_*.js; index_msi_*.html written once). A
# Slurm job, in the python image of the workflow's own python tasks.
index() {
	srun --partition=prod --ntasks=1 --cpus-per-task=1 --mem=1G --job-name=msi-index \
		apptainer exec --bind /media/jp/SSD2/CHU \
		"$WDL_DIR/apptainer-recipes/python-3.12.sif" \
		python3 "$WDL_DIR/scripts/msi_detection_pp.py" index --viz "$VIZ" \
		--run-outputs "$(realpath "$WORK/sprocket/runs/msi_detection/_latest")/outputs.json"
}

# 2. Run. sprocket works in $WORK/sprocket/runs/msi_detection/<timestamp>/:
#    one folder per task (command, stdout, stderr, work/ with its files).
#    index only after a successful run (_latest is the run just made).
if [ -d "$PARAMS" ]; then
	LIST=$(find "$PARAMS" -name '*.json' -type f | sort)
else
	LIST=$PARAMS
fi
[ -n "$LIST" ] || { echo "no *.json in $PARAMS" >&2; exit 1; }

FAILED=""
for p in $LIST; do
	echo "=== $p"
	if sprocket run -v -s -c "$CONFIG" -o "$WORK/sprocket" "$WDL" "@$p" && index; then
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

# Run 1 : Run141_Juno
#LOG_FILE=$WORK/msi_sprocket_$(date +%Y-%m-%d_%H-%M-%S).log
# sprocket run -v -s \
#   -c /media/jp/SSD2/CHU/code/pipeline-jp/wdl/backends.conf/sprocket.local.toml \
#   -o /media/jp/SSD2/CHU/test/work/sprocket \
#   /media/jp/SSD2/CHU/code/pipeline-jp/wdl/msi_detection-0.1.0.wdl \
#   @/media/jp/SSD2/CHU/code/pipeline-jp/wdl/params/example_params_1.json

# index


# Run 2 : Run144_Juno
# LOG_FILE=$WORK/msi_sprocket_$(date +%Y-%m-%d_%H-%M-%S).log
# sprocket run -v -s \
#   -c /media/jp/SSD2/CHU/code/pipeline-jp/wdl/backends.conf/sprocket.local.toml \
#   -o /media/jp/SSD2/CHU/test/work/sprocket \
#   /media/jp/SSD2/CHU/code/pipeline-jp/wdl/msi_detection-0.1.0.wdl \
#   @/media/jp/SSD2/CHU/code/pipeline-jp/wdl/params/example_params_2.json \
# index