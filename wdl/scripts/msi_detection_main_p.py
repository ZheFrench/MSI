#!/usr/bin/env python3
"""
msi_detection_main_p.py

Python steps run inside the msi_detection-0.1.0.wdl workflow, before the
msisensor-pro analysis: they prepare the reference microsatellite list. One
subcommand each (the report and the dashboard are in msi_detection_post_p.py):

    list-to-bed  LIST OUT_BED
        MSIsensor-pro .list -> BED6 (one interval per microsatellite site).
    mark-list    LIST PANEL_BED OUT_LIST
        Full .list with "filter" = PASS for sites in PANEL_BED, NO otherwise
        (the list msisensor-pro pro analyses when a panel is given).

Usage (normally only from msi_detection-0.1.0.wdl):
    sprocket run msi_detection-0.1.0.wdl @params/example_params_1.json
"""

import argparse
import logging
from pathlib import Path

logger = logging.getLogger("MSI-detection")


# ---------------------------------------------------------------------------
# Reference list preparation (list-to-bed / mark-list)
# ---------------------------------------------------------------------------

def convert_list_to_bed(list_path: Path, bed_path: Path) -> Path:
    """Convert a MSIsensor-pro *.msisensor_pro.list reference file into a BED6
    file (chrom, start, end, name, score=0, strand=+).

    "location" in the .list is the 1-based start of the repeat; BED is
    0-based half-open, so start = location - 1 and end = start +
    repeat_unit_length * repeat_times.
    """
    logger.info(f"Generating BED file from microsatellite list: {bed_path}")
    count = 0
    with open(list_path) as fh_in, open(bed_path, "w") as fh_out:
        header = fh_in.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh_in:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            chromosome = cols[idx["chromosome"]]
            location = int(cols[idx["location"]])
            repeat_unit_length = int(cols[idx["repeat_unit_length"]])
            repeat_times = int(cols[idx["repeat_times"]])
            repeat_unit_bases = cols[idx["repeat_unit_bases"]]
            left_flank_bases = cols[idx["left_flank_bases"]]
            right_flank_bases = cols[idx["right_flank_bases"]]

            start = location - 1
            end = start + repeat_unit_length * repeat_times
            name = f"{repeat_times}-{repeat_unit_bases}-{left_flank_bases}-{right_flank_bases}"
            fh_out.write(f"{chromosome}\t{start}\t{end}\t{name}\t0\t+\n")
            count += 1
    logger.info(f"BED file written: {bed_path} ({count:,} sites)")
    return bed_path


def read_bed_keys(bed_path: Path) -> set:
    """(chrom, start) of every interval of a correspondence BED built by the
    bedtools_intersect task from this same list -- matching a .list row on
    (chromosome, location - 1) is then exact, no fuzzy matching needed."""
    keys = set()
    with open(bed_path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2:
                keys.add((cols[0], cols[1]))
    return keys


def apply_panel_filter_to_list(list_path: Path, panel_bed_path: Path, output_path: Path) -> Path:
    """Write output_path with every row of the full genome-wide list kept (same
    row count/order as list_path), but the "filter" column set to "PASS" for
    sites selected in panel_bed_path and "NO" for every other site.

    This keeps the panel restriction visible/inspectable directly in the
    list -- every genome-wide site is still there, just tagged.
    """
    keys = read_bed_keys(panel_bed_path)
    pass_count = 0
    total = 0
    with open(list_path) as fh_in, open(output_path, "w") as fh_out:
        header = fh_in.readline()
        fh_out.write(header)
        idx = {name: i for i, name in enumerate(header.rstrip("\n").split("\t"))}
        filter_idx = idx["filter"]
        for line in fh_in:
            cols = line.rstrip("\n").split("\t")
            if not cols or not cols[0]:
                continue
            total += 1
            start = str(int(cols[idx["location"]]) - 1)
            is_panel = (cols[idx["chromosome"]], start) in keys
            cols[filter_idx] = "PASS" if is_panel else "NO"
            if is_panel:
                pass_count += 1
            fh_out.write("\t".join(cols) + "\n")

    logger.info(f"Panel-marked list written: {output_path} ({pass_count:,}/{total:,} sites PASS)")
    return output_path


def main():
    parser = argparse.ArgumentParser(description="Reference list steps of the msi_detection-0.1.0.wdl workflow.")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("list-to-bed", help="MSIsensor-pro .list -> BED6.")
    p.add_argument("list")
    p.add_argument("output")

    p = sub.add_parser("mark-list", help="Full .list, filter=PASS for sites in the panel BED, NO otherwise.")
    p.add_argument("list")
    p.add_argument("bed")
    p.add_argument("output")

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    if args.command == "list-to-bed":
        convert_list_to_bed(Path(args.list), Path(args.output))
    elif args.command == "mark-list":
        apply_panel_filter_to_list(Path(args.list), Path(args.bed), Path(args.output))


if __name__ == "__main__":
    main()
