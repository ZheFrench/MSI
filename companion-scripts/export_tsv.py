#!/usr/bin/env python3
"""
export_tsv.py

Exports the results of a reports folder (VIZ of BASH.sprocket.batch.sh or
BASH.sprocket.cron.sh) as a TSV table, one line per sample, from the dashboard data
the workflow and msi_detection_post_p.py index already computed (nothing is
recomputed):

    run_id, sample_name, call, consensus_roi, consensus_total, pct_region_unstable
        - VIZ/index_msi_sample.js
    <lib>_unstable_region_count, <lib>_unstable_roi_count (per library)
        - VIZ/index_msi_sample.js
    <marker>, <lib>_<marker>_pro_p, <lib>_<marker>_cov (per ROI marker: BAT25...)
        - VIZ/runs/<run_id>/<sample>_index_sample.js: the marker's site_status
          (same for all libraries), then its pro_p and CovReads in each library

Empty cell: no value (no consensus, marker not covered...).

    export_tsv.py --viz VIZ [--output-dir DIR | --output results.tsv]

Written in --output-dir (default: VIZ) as results_msi_<date>.tsv (date:
YYYY-MM-DD_HH-MM-SS, as the params JSON of msi_detection_pre_p.py): each
export is a new file. --output gives the file itself instead.
"""

import argparse
import csv
import logging
import sys
from datetime import datetime
from pathlib import Path

# msi_detection_post_p.py sits in main-scripts/, next to this script's folder.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "main-scripts"))
from msi_detection_post_p import SAMPLES_DATA_FILENAME, read_data_js  # noqa: E402
from msi_detection_report_p import RUNS_DIRNAME, extract_library_suffix  # noqa: E402

logger = logging.getLogger("export-tsv")

SUMMARY_COLUMNS = ["run_id", "sample_name", "call", "consensus_roi", "consensus_total", "pct_region_unstable"]
LIBRARY_COLUMNS = ["unstable_region_count", "unstable_roi_count"]


def read_markers(viz_dir: Path, run_id: str, sample_name: str) -> dict:
    """ROI marker -> its columns ({marker: status, <lib>_<marker>_pro_p,
    <lib>_<marker>_cov}) for a sample, from its <sample>_index_sample.js;
    {} (warning logged) when it is missing."""
    js_path = viz_dir / RUNS_DIRNAME / run_id / f"{sample_name}_index_sample.js"
    if not js_path.is_file():
        logger.warning(f"{run_id}/{sample_name}: no {js_path}, marker columns left empty")
        return {}
    markers = {}
    for row in read_data_js(js_path):
        marker = row.get("roi_marker")
        if not marker:
            continue
        lib = extract_library_suffix(row["library"])
        columns = markers.setdefault(marker, {})
        columns[marker] = row.get("site_status")
        columns[f"{lib}_{marker}_pro_p"] = row.get("pro_p")
        columns[f"{lib}_{marker}_cov"] = row.get("CovReads")
    return markers


def main():
    parser = argparse.ArgumentParser(description="Results of a reports folder as a TSV table, one line per sample.")
    parser.add_argument("--viz", type=Path, required=True, help="Reports folder (holds index_msi_sample.js).")
    parser.add_argument("--output-dir", type=Path, help="Where results_msi_<date>.tsv is written (default: VIZ).")
    parser.add_argument("--output", type=Path, help="TSV written, instead of --output-dir/results_msi_<date>.tsv.")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")

    data = read_data_js(args.viz / SAMPLES_DATA_FILENAME)
    libraries = data["library_labels"]
    samples = data["samples"]
    markers_by_sample = [read_markers(args.viz, s["run_id"], s["sample_name"]) for s in samples]
    markers = sorted({m for sample_markers in markers_by_sample for m in sample_markers})

    columns = (SUMMARY_COLUMNS
               + [f"{lib}_{col}" for lib in libraries for col in LIBRARY_COLUMNS]
               + [col for m in markers
                  for col in [m] + [f"{lib}_{m}_{v}" for lib in libraries for v in ("pro_p", "cov")]])
    date = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    output = args.output or (args.output_dir or args.viz) / f"results_msi_{date}.tsv"
    with output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sample, sample_markers in zip(samples, markers_by_sample):
            line = {col: sample.get(col) for col in SUMMARY_COLUMNS}
            for lib in libraries:
                for col in LIBRARY_COLUMNS:
                    line[f"{lib}_{col}"] = sample.get("libraries", {}).get(lib, {}).get(col)
            for marker_columns in sample_markers.values():
                line.update(marker_columns)
            writer.writerow(line)
    logger.info(f"{len(samples)} sample(s), {len(markers)} ROI marker(s) -> {output}")


if __name__ == "__main__":
    main()
