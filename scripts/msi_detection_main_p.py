#!/usr/bin/env python3
"""
msi_detection_main_p.py

Python steps run inside the msi_detection-0.1.0.wdl workflow (msi_main_helper,
a File given to its tasks), one subcommand each. The dashboard of every run is
built after the workflow, outside it, by msi_detection_post_p.py index.

    list-to-bed  LIST OUT_BED
        MSIsensor-pro .list -> BED6 (one interval per microsatellite site).
    mark-list    LIST PANEL_BED OUT_LIST
        Full .list with "filter" = PASS for sites in PANEL_BED, NO otherwise
        (the list msisensor-pro pro analyses when a panel is given).

    report       --run-id RUN --samples samples.tsv --genome FASTA
                 --data-output DIR --viz-output DIR --log-output DIR
                 --min-reads N --pval-cutoff P --sample-na-fraction F
                 --consensus-roi-cutoff C
                 [--roi-bed ROI_BED]
        Copies every sample's raw msisensor-pro output under
        <data-output>/<run_id>/<library>/, then generates this run's HTML
        reports (never changed afterwards):

    <viz-output>/runs/<run_id>/<run_id>_index_run.html            - one row per sample, one column-group per library
    <viz-output>/runs/<run_id>/<sample>_index_sample.html         - that sample's result pooled across its libraries

samples.tsv (written by the workflow, no header) has one line per BAM:
    sample_name  library  bam_path  all  unstable  dis  summary
library is the BAM's folder relative to input_dir ("." when the BAM sits
directly in it). A failed BAM fails the workflow before report runs.

The script must stay self-contained (no import of the other scripts): in a
WDL task it is localized on its own.

Usage (normally only from msi_detection-0.1.0.wdl):
    sprocket run msi_detection-0.1.0.wdl @<params JSON>
"""

import argparse
import csv
import json
import logging
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import NoReturn
from urllib.parse import quote

# Set by report from the workflow's params.p_process (min_reads, pval_cutoff, sample_na_fraction,
# consensus_roi_cutoff).
MIN_READS = 200  # a site counts (% unstable, site status) only with >= this many reads
MSI_PVAL_CUTOFF = 0.1  # a site is unstable when pro_p > this
SAMPLE_NA_FRACTION = 0.5  # sample status NA when more than this fraction of its listed loci are NA or MSU
MSI_CONSENSUS_ROI_CUTOFF = 3  # consensus ROI unstable sites >= this -> MSI

# Columns of msisensor-pro's *_all / *_unstable files
UNSTABLE_COLUMNS = [
    "chromosome", "location", "left_flank_bases", "repeat_times",
    "repeat_unit_bases", "right_flank_bases", "pro_p", "pro_q",
    "CovReads",
]
# Columns of the workflow's samples TSV (one row per BAM)
SAMPLES_TSV_COLUMNS = ["sample_name", "library", "bam_path", "all", "unstable", "dis", "summary"]

# Reports folder: viz/runs/<run_id>/ per run, and the all-runs dashboard
# (its data files: msi_detection_post_p.py)
RUNS_DIRNAME = "runs"
OUTPUT_INDEX_FILENAME = "index_msi_run.html"

logger = logging.getLogger("MSI-detection")


def setup_logging(run_id: str, log_dir: Path) -> Path:
    """Configure the module logger: console + a per-run file under log_dir."""
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_path = log_dir / f"{run_id}_MSI_{timestamp}.log"

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")

    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return log_path


def fatal(msg: str) -> NoReturn:
    """Log an error and exit(1). Use instead of sys.exit() once logging is set up."""
    logger.error(msg)
    sys.exit(1)


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


# ---------------------------------------------------------------------------
# HTML reporting: RUN index -> SAMPLE index -> per-sample RESULT report
# ---------------------------------------------------------------------------

def parse_roi_bed(roi_bed_path: Path) -> dict:
    """Map (chromosome, location) -> ROI marker name from a *.roi.bed /
    *.panel.roi.bed correspondence file, for highlighting a sample's own
    unstable sites that land on a classic MSI marker.

    location = start + 1 (BED is 0-based; .list/._unstable "location" is
    1-based), matching how generate_panel_correspondence_bed() built it. The
    marker name is the last "#"-separated segment of the combined name field
    (it may already carry a panel region name ahead of it).
    """
    markers = {}
    if not roi_bed_path or not roi_bed_path.exists():
        return markers
    with open(roi_bed_path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            markers[(cols[0], int(cols[1]) + 1)] = cols[3].rsplit("#", 1)[-1]
    return markers


def parse_msisensor_tsv(path: Path) -> list:
    """Parse a MSIsensor-pro *_all or *_unstable TSV file into a list of typed
    dicts -- both share the same column layout, only *_all also includes the
    sites MSIsensor-pro didn't call unstable."""
    if not path.exists():
        return []
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}

        def get(cols, name):
            i = idx.get(name)
            return cols[i] if i is not None and i < len(cols) else ""

        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            rows.append({
                "chromosome": get(cols, "chromosome"),
                "location": int(get(cols, "location") or 0),
                "left_flank_bases": get(cols, "left_flank_bases"),
                "repeat_times": int(get(cols, "repeat_times") or 0),
                "repeat_unit_bases": get(cols, "repeat_unit_bases"),
                "right_flank_bases": get(cols, "right_flank_bases"),
                "pro_p": float(get(cols, "pro_p") or 0),
                "pro_q": float(get(cols, "pro_q") or 0),
                "CovReads": int(get(cols, "CovReads") or 0),
            })
    return rows


def row_is_unstable(row: dict) -> bool:
    """Unstable call for one *_all row: pro_p > MSI_PVAL_CUTOFF (params
    pval_cutoff). At 0.1 it is msisensor-pro's own call (its threshold
    column, *_unstable holds exactly the sites with pro_p > threshold)."""
    return row["pro_p"] > MSI_PVAL_CUTOFF


def row_is_counted(row: dict) -> bool:
    """A site that counts towards % unstable: read >= MIN_READS times (and
    actually evaluated: msisensor-pro reports pro_p as NaN at 0 reads)."""
    return row["CovReads"] >= MIN_READS and row["pro_p"] == row["pro_p"]


def counted_site_keys(rows) -> set:
    return {(r["chromosome"], r["location"]) for r in rows if row_is_counted(r)}


def counted_unstable_keys(rows) -> set:
    return {(r["chromosome"], r["location"]) for r in rows if row_is_counted(r) and row_is_unstable(r)}


def pct(part: int, total: int):
    return round(100 * part / total, 1) if total else None


REPORT_CSS = """
:root{
  --bg:#f5f6f8; --surface:#ffffff; --surface-alt:#eef1f4; --border:#dbe0e6;
  --text:#16202b; --text-muted:#5c6773;
  --accent:#0f7d6c; --accent-soft:#e3f3ef;
  --sig-high:#0f7d6c; --sig-high-soft:#e3f3ef;
  --sig-mid:#a8790a; --sig-mid-soft:#f7eed9;
  --sig-low:#6b7680; --sig-low-soft:#eceef0;
  --sig-danger:#b3261e; --sig-danger-soft:#fbe9e7;
  --lib-color-0:#2f6fed; --lib-color-1:#e2792e; --lib-color-2:#0f7d6c; --lib-color-3:#a8790a;
  --shadow: 0 1px 2px rgba(20,30,40,.06), 0 8px 24px -12px rgba(20,30,40,.12);
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --bg:#0f151b; --surface:#161e26; --surface-alt:#1c2530; --border:#2a3641;
    --text:#e7ecf0; --text-muted:#94a2ad;
    --accent:#3bc0a4; --accent-soft:rgba(59,192,164,.14);
    --sig-high:#3bc0a4; --sig-high-soft:rgba(59,192,164,.14);
    --sig-mid:#dcab3f; --sig-mid-soft:rgba(220,171,63,.14);
    --sig-low:#8b98a3; --sig-low-soft:rgba(139,152,163,.14);
    --sig-danger:#e5726b; --sig-danger-soft:rgba(229,114,107,.16);
    --lib-color-0:#6d9bff; --lib-color-1:#f0a35f; --lib-color-2:#3bc0a4; --lib-color-3:#dcab3f;
    --shadow: 0 1px 2px rgba(0,0,0,.3), 0 8px 24px -12px rgba(0,0,0,.5);
  }
}
:root[data-theme="dark"]{
  --bg:#0f151b; --surface:#161e26; --surface-alt:#1c2530; --border:#2a3641;
  --text:#e7ecf0; --text-muted:#94a2ad;
  --accent:#3bc0a4; --accent-soft:rgba(59,192,164,.14);
  --sig-high:#3bc0a4; --sig-high-soft:rgba(59,192,164,.14);
  --sig-mid:#dcab3f; --sig-mid-soft:rgba(220,171,63,.14);
  --sig-low:#8b98a3; --sig-low-soft:rgba(139,152,163,.14);
  --sig-danger:#e5726b; --sig-danger-soft:rgba(229,114,107,.16);
  --lib-color-0:#6d9bff; --lib-color-1:#f0a35f; --lib-color-2:#3bc0a4; --lib-color-3:#dcab3f;
  --shadow: 0 1px 2px rgba(0,0,0,.3), 0 8px 24px -12px rgba(0,0,0,.5);
}
*{box-sizing:border-box;}
body{margin:0;background:var(--bg);color:var(--text);font-family:"IBM Plex Sans",ui-sans-serif,system-ui,sans-serif;
  font-size:14px;line-height:1.5;padding-inline:20px;padding-block:28px 60px;}
.wrap{max-width:1120px;margin:0 auto;display:flex;flex-direction:column;gap:20px;}
code,.mono{font-family:"IBM Plex Mono",ui-monospace,SFMono-Regular,monospace;}
.tabular{font-variant-numeric:tabular-nums;}
a{color:var(--accent);}
h1{text-wrap:balance;margin:0 0 4px 0;font-size:1.5rem;font-weight:700;letter-spacing:-.01em;}
.breadcrumb{font-size:.8rem;color:var(--text-muted);margin-bottom:2px;}
.breadcrumb a{color:var(--text-muted);text-decoration:none;}
.breadcrumb a:hover{color:var(--accent);}
header{border-bottom:1px solid var(--border);padding-bottom:16px;display:flex;flex-direction:column;gap:6px;}
header p{margin:0;color:var(--text-muted);max-width:70ch;}
.meta-row{display:flex;flex-wrap:wrap;gap:8px;margin-top:6px;}
.meta-chip{display:inline-flex;align-items:center;gap:5px;background:var(--surface-alt);border:1px solid var(--border);
  border-radius:6px;padding:3px 9px;font-size:.78rem;color:var(--text-muted);}
.meta-chip strong{color:var(--text);font-weight:600;}
.meta-chip.consensus{border-color:var(--text);}
.meta-note{font-size:12px;color:var(--text-muted);line-height:1.5;}
.panel-row{display:grid;grid-template-columns:1.3fr 1fr;gap:16px;}
@media (max-width:820px){.panel-row{grid-template-columns:1fr;}}
.panel{background:var(--surface);border:1px solid var(--border);border-radius:12px;padding:16px 18px;box-shadow:var(--shadow);}
.panel h2{margin:0 0 12px 0;font-size:.82rem;text-transform:uppercase;letter-spacing:.06em;color:var(--text-muted);font-weight:600;}
#chart{width:100%;height:170px;overflow:visible;}
#chart .bar{fill:var(--accent);}
#chart .bar:hover{fill:var(--sig-mid);}
#chart .axis{stroke:var(--border);}
#chart text{fill:var(--text-muted);font-size:10px;font-family:"IBM Plex Mono",monospace;}
.field{display:flex;flex-direction:column;gap:6px;margin-bottom:14px;}
.field:last-child{margin-bottom:0;}
.field label{font-size:.78rem;color:var(--text-muted);display:flex;justify-content:space-between;}
.field label .val{font-family:"IBM Plex Mono",monospace;color:var(--text);}
.field input[type="text"]{border:1px solid var(--border);border-radius:8px;background:var(--bg);color:var(--text);
  padding:8px 10px;font:inherit;width:100%;}
input[type="range"]{width:100%;accent-color:var(--accent);}
.reset-row{display:flex;justify-content:flex-end;margin-top:4px;}
.reset-row button{background:none;border:none;color:var(--accent);font:inherit;font-weight:500;cursor:pointer;padding:4px 0;}
.reset-row button:hover{text-decoration:underline;}
.legend{display:flex;gap:14px;flex-wrap:wrap;margin-top:14px;padding-top:12px;border-top:1px solid var(--border);}
.legend .item{display:flex;align-items:center;gap:6px;font-size:.76rem;color:var(--text-muted);}
.legend .dot{width:9px;height:9px;border-radius:50%;flex:none;}
.table-panel{background:var(--surface);border:1px solid var(--border);border-radius:12px;box-shadow:var(--shadow);overflow:hidden;}
.table-scroll{overflow-x:auto;}
table{border-collapse:collapse;width:100%;min-width:760px;}
thead th{position:sticky;top:0;background:var(--surface-alt);text-align:left;font-size:.72rem;text-transform:uppercase;
  letter-spacing:.05em;color:var(--text-muted);padding:10px 12px;border-bottom:1px solid var(--border);white-space:nowrap;
  cursor:pointer;user-select:none;}
thead th:hover{color:var(--text);}
thead th .arrow{opacity:.5;margin-left:3px;font-size:.7em;}
thead th.active{color:var(--accent);}
thead th.active .arrow{opacity:1;}
tbody td{padding:9px 12px;border-bottom:1px solid var(--border);vertical-align:middle;white-space:nowrap;}
tbody tr:last-child td{border-bottom:none;}
tbody tr:hover{background:var(--surface-alt);}
.cell-mono{font-family:"IBM Plex Mono",monospace;font-size:.78rem;}
.lib-c0{color:var(--lib-color-0);}
.lib-c1{color:var(--lib-color-1);}
.lib-c2{color:var(--lib-color-2);}
.lib-c3{color:var(--lib-color-3);}
.locus{font-family:"IBM Plex Mono",monospace;font-size:.82rem;color:var(--text-muted);}
.locus b{color:var(--text);font-weight:700;background:var(--accent-soft);border-radius:3px;padding:0 2px;}
.chip{display:inline-flex;align-items:center;gap:5px;border-radius:999px;padding:3px 10px;font-size:.74rem;font-weight:600;}
.chip.high{background:var(--sig-high-soft);color:var(--sig-high);}
.chip.mid{background:var(--sig-mid-soft);color:var(--sig-mid);}
.chip.low{background:var(--sig-low-soft);color:var(--sig-low);}
.chip.danger{background:var(--sig-danger-soft);color:var(--sig-danger);}
.chip.roi{background:var(--accent-soft);color:var(--accent);}
tr.roi-hit td{background:var(--accent-soft);}
tr.roi-hit:hover td{background:var(--accent-soft);filter:brightness(0.97);}
tr.roi-hit td:first-child{box-shadow:inset 3px 0 0 var(--accent);}
.lib-cov{display:block;font-size:.68rem;color:var(--text-muted);margin-top:2px;}
.lib-empty{color:var(--text-muted);}
tr.site-group td{border-top:2px solid var(--border);}
tr.site-group:first-child td{border-top:none;}
.num{font-family:"IBM Plex Mono",monospace;text-align:right;}
.table-foot{display:flex;align-items:center;justify-content:space-between;gap:12px;flex-wrap:wrap;padding:12px 16px;
  border-top:1px solid var(--border);background:var(--surface-alt);font-size:.82rem;color:var(--text-muted);}
.pager{display:flex;gap:6px;align-items:center;}
.pager button{border:1px solid var(--border);background:var(--surface);color:var(--text);border-radius:6px;
  padding:5px 10px;font:inherit;cursor:pointer;}
.pager button:disabled{opacity:.4;cursor:default;}
.pager button:not(:disabled):hover{border-color:var(--accent);color:var(--accent);}
.empty-row td{text-align:center;color:var(--text-muted);padding:32px 12px;}
.status-pill{display:inline-flex;align-items:center;gap:5px;border-radius:999px;padding:3px 10px;font-size:.76rem;font-weight:600;}
.status-pill.ok{background:var(--sig-high-soft);color:var(--sig-high);}
.status-pill.failed{background:var(--sig-mid-soft);color:var(--sig-mid);}
.index-table{width:100%;border-collapse:collapse;}
.index-table th{text-align:left;font-size:.72rem;text-transform:uppercase;letter-spacing:.05em;color:var(--text-muted);
  padding:10px 12px;border-bottom:1px solid var(--border);}
.index-table th.num{text-align:right;}
.index-table th.consensus{color:var(--text);}
.index-table td.consensus{font-weight:700;}
.index-table td{padding:11px 12px;border-bottom:1px solid var(--border);}
.index-table tr:last-child td{border-bottom:none;}
.index-table tr:hover td{background:var(--surface-alt);}
footer{color:var(--text-muted);font-size:.78rem;text-align:center;padding-top:6px;}
"""

# The page icon, embedded here rather than read from a file next to the script:
# in a WDL task the script is localized on its own (msi_helper is a File).
FAVICON_SVG = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64">
  <defs>
    <linearGradient id="g" x1="0" y1="0" x2="1" y2="1">
      <stop offset="0" stop-color="#4F46E5"/>
      <stop offset="1" stop-color="#0EA5E9"/>
    </linearGradient>
  </defs>
  <rect width="64" height="64" rx="14" fill="url(#g)"/>
  <text x="32" y="41" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" font-weight="800" font-size="24" letter-spacing="-0.5" fill="#FFFFFF">MSI</text>
</svg>"""


def favicon_link() -> str:
    """<link rel="icon"> embedding FAVICON_SVG as a data URI, so every page
    carries its icon wherever viz/ ends up and even opened on its own."""
    return f'<link rel="icon" type="image/svg+xml" href="data:image/svg+xml,{quote(FAVICON_SVG.strip())}">'


# Shared <head> links of every page: web fonts + favicon.
HEAD_LINKS = ('<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;500;600&display=swap">'
              + favicon_link())


def extract_library_suffix(lib: str) -> str:
    """"alignements_libA" -> "libA": strips everything before the trailing
    "lib<x>" segment (e.g. the "alignements_" folder-naming prefix), falling
    back to the original string if it doesn't end in "lib<x>"."""
    m = re.search(r"(lib[A-Za-z0-9]+)$", lib, re.IGNORECASE)
    return m.group(1) if m else lib


def sample_status_rule() -> str:
    """Sample status rule, shown on the run pages and the dashboard."""
    return (f"Sample status: NA when more than {round(100 * SAMPLE_NA_FRACTION)}% of the loci listed on its page are NA or MSU; "
            f"otherwise MSI when &ge;{MSI_CONSENSUS_ROI_CUTOFF} ROI markers are unstable in consensus, MSS if not.")


def site_status(libs: list, ok_lib_count: int) -> str:
    """One locus's status from its rows in the sample's ok libraries (same
    value the sample page's Status column shows): NA when no library read it
    >= MIN_READS times, MSU when only some did; otherwise MSI when every
    library calls it unstable, MSS when none does, MSU when they disagree."""
    counted = [r for r in libs if row_is_counted(r)]
    if not counted:
        return "NA"
    if len(counted) < ok_lib_count:
        return "MSU"
    passing = sum(1 for r in counted if row_is_unstable(r))
    if passing == len(counted):
        return "MSI"
    return "MSS" if passing == 0 else "MSU"


def build_run_table(samples: list, rows_by_sample: dict = None) -> dict:
    """Precompute (Phase 2) everything <run_id>_index_run.html displays: one row
    per sample, one column-group per library (e.g. alignements_libA,
    alignements_libB) side by side instead of a row per (sample, library) pair.

    Per library: unstable-site count over the whole analysed region
    (unstable_region_count), and how many of those land on a classic MSI marker
    (unstable_roi_count); a failed library is just {"failed": true}. Site and
    total counts only include sites with >= MIN_READS reads (ROI counts don't).
    When a sample has 2+ ok libraries, consensus counts give how many unstable
    sites (total, and ROI-only) every one of its libraries shares, "% unstable"
    (pct_region_unstable) is consensus_total over the sites every library
    covers (the lone library's own % otherwise), and "call" is MSI when the
    consensus ROI count reaches MSI_CONSENSUS_ROI_CUTOFF."""
    rows_by_sample = rows_by_sample or {}
    libraries = sorted({s.get("library") or "-" for s in samples})

    by_sample = {}
    for s in samples:
        by_sample.setdefault(s["sample_name"], {})[s.get("library") or "-"] = s

    # rows come from *_all (every evaluated site, not just the unstable ones),
    # so these need their own unstable filter.
    def roi_site_keys(rows):
        return {(r["chromosome"], r["location"]) for r in rows if row_is_unstable(r) and r.get("roi_marker")}

    def roi_covered(rows):
        # msisensor-pro reports pro_p as NaN for a site with 0 reads
        return any(r.get("roi_marker") and r["CovReads"] > 0 and r["pro_p"] == r["pro_p"] for r in rows)

    table = []
    for sname, entries in by_sample.items():
        lib_rows = rows_by_sample.get(sname, {})
        libs = {}
        for lib in libraries:
            entry = entries.get(lib)
            if entry is None:
                continue
            if entry.get("failed"):
                libs[lib] = {"failed": True}
            else:
                libs[lib] = {
                    "unstable_region_count": entry["unstable_region_count"],
                    "unstable_roi_count": len(roi_site_keys(lib_rows.get(lib, []))),
                }

        ok_libs = [lib for lib in libraries if lib in entries and not entries[lib].get("failed")]
        ok_pcts = [entries[lib]["pct_region_unstable"] for lib in ok_libs
                   if entries[lib].get("pct_region_unstable") is not None]
        row = {
            "sample_name": sname,
            "pct_region_unstable": ok_pcts[0] if len(ok_pcts) == 1 else None,
            "libraries": libs,
            "consensus_total": None,
            "consensus_roi": None,
            "call": None,
        }
        if len(ok_libs) >= 2:
            # Consensus % unstable: sites every library read >= MIN_READS times
            # (common), and how many of those all of them call unstable.
            common = set.intersection(*(counted_site_keys(lib_rows.get(lib, [])) for lib in ok_libs))
            row["consensus_total"] = len(set.intersection(*(counted_unstable_keys(lib_rows.get(lib, [])) for lib in ok_libs)))
            row["pct_region_unstable"] = pct(row["consensus_total"], len(common))
            row["consensus_roi"] = len(set.intersection(*(roi_site_keys(lib_rows.get(lib, [])) for lib in ok_libs)))
            loci = {}
            for lib in ok_libs:
                for r in lib_rows.get(lib, []):
                    loci.setdefault((r["chromosome"], r["location"]), []).append(r)
            # Only the loci its sample page lists (whatever the p-value
            # filter): an ROI marker, or a p-value in at least one library.
            # The others (no p-value anywhere: 0 reads, e.g. not covered by
            # an older panel) say nothing about how interpretable it is.
            loci = {k: libs_at for k, libs_at in loci.items()
                    if libs_at[0]["roi_marker"] or any(r["pro_p"] == r["pro_p"] for r in libs_at)}
            uninterpretable = sum(1 for libs_at in loci.values() if site_status(libs_at, len(ok_libs)) in ("NA", "MSU"))
            if loci and uninterpretable > SAMPLE_NA_FRACTION * len(loci):
                row["call"] = "NA"  # most loci NA or MSU: no call
            elif not any(roi_covered(lib_rows.get(lib, [])) for lib in ok_libs):
                row["call"] = "NA"  # no ROI marker read at all (0 reads everywhere): no call possible
            else:
                row["call"] = "MSI" if row["consensus_roi"] >= MSI_CONSENSUS_ROI_CUTOFF else "MSS"
        table.append(row)

    return {
        "libraries": libraries,
        "library_labels": {lib: extract_library_suffix(lib) for lib in libraries},
        "table": table,
    }


def render_run_index(run_id: str) -> str:
    """Static shell for <run_id>_index_run.html: everything it shows comes from
    its one data file, <run_id>_index_run.js (window.RUN_DATA), loaded via
    <script src>."""
    html = """<title>__TITLE__</title>
__HEAD_LINKS__
<style>__CSS__</style>
<div class="wrap">
  <div class="breadcrumb"><a href="../../__OUTPUT_INDEX_FILE__">&larr; All runs</a></div>
  <header>
    <h1>__RUN_ID__</h1>
    <p><span id="sampleCount">0</span> sample(s) processed successfully.</p>
    <div class="meta-note">__SAMPLE_STATUS_RULE__</div>
  </header>
  <div class="table-panel">
    <div class="table-scroll">
      <table class="index-table">
        <thead id="thead"></thead>
        <tbody id="tbody"></tbody>
      </table>
    </div>
  </div>
  <footer>MSIsensor-pro pipeline report</footer>
</div>
<script src="./__DATA_JS_FILE__"></script>
<script>
(function(){
  "use strict";
  var D = window.RUN_DATA || {libraries: [], library_labels: {}, table: []};
  function n(v){ return (v === null || v === undefined) ? "-" : Number(v).toLocaleString("en-US"); }
  function pct(v){ return (v === null || v === undefined) ? "-" : v + "%"; }
  function num(v){ return '<td class="num tabular">' + v + '</td>'; }
  var FAILED = '<span class="status-pill failed">failed</span>';

  var groups = D.libraries.map(function(lib){
    return '<th class="num" colspan="2">UNSTABLE_SITE_' + D.library_labels[lib].toUpperCase() + '</th>';
  }).join("");
  // UNSTABLE_CONSENSUS (headers in black, values in bold) stands out from the per-library groups.
  var sub = '';
  for (var i = 0; i < D.libraries.length; i++) sub += '<th class="num">Total</th><th class="num">ROI</th>';
  sub += '<th class="num consensus">Total</th><th class="num consensus">ROI</th>';
  var thead = '<tr><th rowspan="2">Sample</th><th rowspan="2">% unstable</th>' + groups
        + '<th class="num consensus" colspan="2">UNSTABLE_CONSENSUS</th><th rowspan="2">Status</th></tr><tr>' + sub + '</tr>';
  var colspan = 1 + 1 + 2 * (D.libraries.length + 1) + 1;
  var rows = D.table.map(function(r){
    var cells = D.libraries.map(function(lib){
      var e = r.libraries[lib];
      if (!e) return num("-") + num("-");
      if (e.failed) return '<td colspan="2">' + FAILED + '</td>';
      return num(n(e.unstable_region_count)) + num(e.unstable_roi_count);
    }).join("");
    var status = (r.call && r.call !== "NA") ? '<span class="chip ' + (r.call === "MSI" ? "mid" : "low") + '">' + r.call + '</span>'
               : r.call === "NA" ? '<span class="lib-empty">NA</span>'
                        : '<span class="lib-empty">—</span>';
    return '<tr><td><a href="./' + r.sample_name + '_index_sample.html">' + r.sample_name + '</a></td>'
      + num(pct(r.pct_region_unstable)) + cells
      + '<td class="num tabular consensus">' + n(r.consensus_total) + '</td><td class="num tabular consensus">' + n(r.consensus_roi) + '</td>'
      + '<td>' + status + '</td></tr>';
  });
  document.getElementById("sampleCount").textContent = D.table.length;
  document.getElementById("thead").innerHTML = thead;
  document.getElementById("tbody").innerHTML = rows.join("")
    || '<tr class="empty-row"><td colspan="' + colspan + '">No sample.</td></tr>';
})();
</script>
"""
    return (html
            .replace("__TITLE__", f"Run {run_id}")
            .replace("__HEAD_LINKS__", HEAD_LINKS)
            .replace("__CSS__", REPORT_CSS)
            .replace("__RUN_ID__", run_id)
            .replace("__OUTPUT_INDEX_FILE__", OUTPUT_INDEX_FILENAME)
            .replace("__SAMPLE_STATUS_RULE__", sample_status_rule())
            .replace("__DATA_JS_FILE__", f"{run_id}_index_run.js"))


def render_sample_group_report(run_id, sample_name, generated_at, rows, entries_by_library) -> str:
    """The only report a sample gets (run_dir/<sample>/index_sample.html):
    every unstable site from every library's *_unstable file, pooled into one
    table keyed on genomic position -- one row per (chromosome, location),
    with the Library/Reads/p-value columns concatenating each library's own
    value in library order (e.g. "libA | libB"). Repeats shows a single value
    (the reference locus's own repeat count from the genome scan, identical
    across libraries for the same site) rather than repeating it per library.
    A p-value filter (default: keep a site if at least one library's p-value
    is above 0.1) is applied client-side. Sites on a classic MSI marker (ROI)
    are highlighted and sorted first by default."""
    genome = next((e.get("genome") for e in entries_by_library.values() if e.get("genome")), "-")
    # Each library gets its own color (cycling through 4) instead of a boxed
    # badge, so alignements_libA/alignements_libB are visually distinguished
    # by color rather than by an outline.
    lib_color_json = json.dumps({lib: i % 4 for i, lib in enumerate(sorted(entries_by_library))})

    # % unstable per library, over its sites with >= MIN_READS reads only
    # (numerator and denominator alike). Consensus: sites every ok library read
    # >= MIN_READS times (common), and how many of those all of them call
    # unstable (common MSI) -- the run page's "% unstable".
    def pct_chip(label, unstable, total, cls="meta-chip"):
        value = pct(unstable, total)
        return (f'<span class="{cls}">{label}: <strong>{"-" if value is None else f"{value}%"} unstable</strong> '
                f'({unstable}/{total})</span> ')

    covered_by_lib, unstable_by_lib = {}, {}
    library_chips = []
    for library in sorted(entries_by_library):
        entry = entries_by_library[library]
        lib_label = extract_library_suffix(library)
        if not entry.get("failed"):
            lib_rows = [r for r in rows if r.get("library") == library]
            covered_by_lib[library] = counted_site_keys(lib_rows)
            unstable_by_lib[library] = counted_unstable_keys(lib_rows)
            library_chips.append(pct_chip(lib_label, len(unstable_by_lib[library]), len(covered_by_lib[library])))
        else:
            library_chips.append(
                f'<span class="meta-chip">{lib_label}: <span class="status-pill failed">failed</span></span>'
            )
    # Consensus on its own line, under the libraries, outlined.
    consensus_row = ""
    if len(covered_by_lib) >= 2:
        common = set.intersection(*covered_by_lib.values())
        common_msi = set.intersection(*unstable_by_lib.values())
        consensus_row = ('<div class="meta-row">'
                         + pct_chip("Consensus", len(common_msi), len(common), "meta-chip consensus") + '</div>')
    consensus_row += (f'<div class="meta-row"><div class="meta-note">'
                      f'Instability defined by &ge;{MIN_READS} reads &amp; p&gt;{MSI_PVAL_CUTOFF}.<br>'
                      f'MSU when one library is uninterpretable (&lt;{MIN_READS} reads) or the libraries '
                      f'disagree (p on either side of {MSI_PVAL_CUTOFF}).<br>'
                      f'NA when both libraries have &lt;{MIN_READS} reads.</div></div>')

    # rows come from *_all now (every evaluated site, per library), so these
    # footer counts need their own unstable filter to keep meaning "unstable
    # site(s)" rather than "every site MSIsensor-pro looked at".
    roi_hits = sum(1 for r in rows if row_is_unstable(r) and r.get("roi_marker"))
    lib_count = len(entries_by_library)
    site_count = len({(r["chromosome"], r["location"]) for r in rows if row_is_unstable(r)})

    html = """<title>__TITLE__</title>
__HEAD_LINKS__
<style>__CSS__</style>
<div class="wrap">
  <div class="breadcrumb"><a href="./__RUN_INDEX_FILE__">&larr; Run __RUN_ID__</a></div>
  <header>
    <h1>__SAMPLE__</h1>
    <div class="meta-row">
      __LIBRARY_CHIPS__
      <span class="meta-chip">Generated: <strong>__GENERATED__</strong></span>
    </div>
    __CONSENSUS_ROW__
  </header>

  <div class="panel">
    <h2>Filters</h2>
    <div class="field">
      <label for="pvalInput">Min p-value (kept if any library clears it; ROI always shown) <span class="val" id="pvalVal">0.1</span></label>
      <input type="text" id="pvalInput" value="0.1" inputmode="decimal">
    </div>
    <div class="reset-row"><button id="resetBtn" type="button">Reset filters</button></div>
    <div class="legend">
      <span class="item"><span class="dot" style="background:var(--accent)"></span>ROI marker</span>
    </div>
  </div>

  <div class="table-panel">
    <div class="table-scroll">
      <table><thead><tr id="headRow"></tr></thead><tbody id="tbody"></tbody></table>
    </div>
    <div class="table-foot">
      <span id="rangeLabel">-</span>
      <div class="pager">
        <button id="prevBtn" type="button">&larr; Prev</button>
        <span id="pageLabel">Page 1</span>
        <button id="nextBtn" type="button">Next &rarr;</button>
      </div>
    </div>
  </div>
  <footer>__COUNT__ site(s) across __LIBCOUNT__ librar__LIBPLURAL__</footer>
</div>

<script src="./__DATA_JS_FILE__"></script>
<script>
(function(){
  "use strict";
  var ROWS = window.ROWS || [];
  var LIB_COLOR_INDEX = __LIB_COLOR_JSON__;
  var COLS = [
    {key:"roi_marker", label:"ROI", type:"chip"},
    {key:"chromosome", label:"Position", type:"range"},
    {key:"library", label:"Library", type:"text"},
    {key:"locus", label:"Locus", type:"locus", noSort:true},
    {key:"repeat_times", label:"Repeats", type:"num"},
    {key:"CovReads", label:"Reads", type:"num"},
    {key:"pro_p", label:"p-value", type:"sci"},
    {key:"status", label:"Status", type:"chip"}
  ];
  var CHROM_ORDER = (function(){ var m={}; for(var i=1;i<=22;i++) m["chr"+i]=i; m.chrX=23;m.chrY=24;m.chrM=25;m.chrMT=25; return m; })();

  function chromRank(c){ if (CHROM_ORDER[c]!=null) return CHROM_ORDER[c]; var n=parseInt(c.replace(/^chr/i,""),10); return isNaN(n)?99:n; }
  function toSubscript(n){ var map={"0":"\\u2080","1":"\\u2081","2":"\\u2082","3":"\\u2083","4":"\\u2084","5":"\\u2085","6":"\\u2086","7":"\\u2087","8":"\\u2088","9":"\\u2089"}; return String(n).split("").map(function(d){return map[d]||d;}).join(""); }
  function fmtSci(v){ if (v===null || isNaN(v)) return "-"; if (v===0) return "0"; if (Math.abs(v)<0.0001) return v.toExponential(2); return v.toFixed(4); }
  function shortLib(lib){ var m=/(lib[A-Za-z0-9]+)$/i.exec(lib); return m?m[1]:lib; }

  // Position is the reference key: every *_unstable row for the same
  // (chromosome, location) across libA/libB/... becomes ONE row, with
  // Library/Repeats/Reads/p-value each holding every library's own value,
  // concatenated in the same library order (sorted by short name, e.g.
  // libA before libB).
  function groupBySite(rows){
    var map={}, order=[];
    rows.forEach(function(r){
      var key=r.chromosome+":"+r.location;
      if (!map[key]){
        map[key]={
          chromosome:r.chromosome, location:r.location,
          repeat_unit_bases:r.repeat_unit_bases,
          left_flank_bases:r.left_flank_bases, right_flank_bases:r.right_flank_bases,
          roi_marker:r.roi_marker, status:r.site_status||"NA", libs:[]
        };
        order.push(key);
      }
      map[key].libs.push(r);
    });
    order.forEach(function(key){
      map[key].libs.sort(function(a,b){ var la=shortLib(a.library), lb=shortLib(b.library); return la<lb?-1:(la>lb?1:0); });
    });
    return order.map(function(key){ return map[key]; });
  }
  var GROUPED = groupBySite(ROWS);

  // 1-based inclusive end, matching how the .bed correspondence files are built
  // (start = location-1, end = start + unit_length*repeat_times) -- so pasting
  // "chr:start-end" straight into IGV's locus box lands exactly on the repeat.
  // Uses the largest repeat_times seen across this site's libraries, so the
  // range covers every library's observed repeat span.
  function maxRepeatTimes(g){ return g.libs.reduce(function(m,l){ return Math.max(m,l.repeat_times); },0); }
  function igvRange(g){ return g.chromosome+":"+g.location+"-"+(g.location+g.repeat_unit_bases.length*maxRepeatTimes(g)-1); }

  // Site status (MSI/MSS/MSU/NA) is precomputed by site_status() in Python
  // and shipped in the data file as each row's site_status -- fixed at
  // MSI_PVAL_CUTOFF regardless of state.minPval, which only controls which
  // rows are shown, not this call.
  var STATUS_CHIP = { MSI:'<span class="chip mid">MSI</span>', MSU:'<span class="chip danger">MSU</span>',
                      MSS:'<span class="chip low">MSS</span>', NA:'<span class="lib-empty">NA</span>' };
  var STATUS_RANK = { MSI:2, MSU:1, MSS:0, NA:-1 };

  // ROI sites (roi_marker set) sort first by default (roi_marker descending
  // puts non-empty marker names ahead of "").
  var state = { filtered: GROUPED.slice(), sortKey:"roi_marker", sortDir:-1, page:0, pageSize:50, minPval:0.1 };

  function applyFilters(){
    state.filtered = GROUPED.filter(function(g){
      // The classic MSI markers (ROI) are always shown, whatever their
      // p-value or coverage. Any other site is kept if at least one library's
      // p-value clears the threshold -- a site only one library flags is
      // still worth surfacing.
      if (g.roi_marker) return true;
      return g.libs.some(function(l){ return l.pro_p > state.minPval; });
    });
    sortRows(); state.page=0; renderTable();
  }
  function sortValue(g, key){
    if (key==="chromosome") return chromRank(g.chromosome);
    if (key==="library") return g.libs.map(function(l){return shortLib(l.library);}).join("|");
    if (key==="repeat_times" || key==="CovReads" || key==="pro_p") return g.libs.length ? g.libs[0][key] : 0;
    if (key==="status") return STATUS_RANK[g.status];
    return g[key];
  }
  function sortRows(){
    var key=state.sortKey, dir=state.sortDir;
    state.filtered.sort(function(a,b){
      var av=sortValue(a,key), bv=sortValue(b,key);
      if (av<bv) return -1*dir; if (av>bv) return 1*dir;
      var ca=chromRank(a.chromosome), cb=chromRank(b.chromosome);
      if (ca!==cb) return ca<cb?-1:1;
      return a.location<b.location?-1:(a.location>b.location?1:0);
    });
  }
  function renderTable(){
    var head=document.getElementById("headRow");
    head.innerHTML=COLS.map(function(c){
      var active=c.key===state.sortKey, arrow=active?(state.sortDir===1?"\\u25B2":"\\u25BC"):"";
      var cls=(c.type==="num"||c.type==="sci"?"num ":"")+(active?"active":"");
      return '<th data-key="'+c.key+'" class="'+cls+'"'+(c.noSort?' data-nosort="1"':'')+'>'+c.label+' <span class="arrow">'+arrow+'</span></th>';
    }).join("");
    var start=state.page*state.pageSize, pageRows=state.filtered.slice(start,start+state.pageSize);
    var body=document.getElementById("tbody");
    if (!pageRows.length){
      body.innerHTML='<tr class="empty-row"><td colspan="'+COLS.length+'">No site matches the current filters.</td></tr>';
    } else {
      body.innerHTML=pageRows.map(function(g){
        var repeatTimes=maxRepeatTimes(g);
        var motif="("+g.repeat_unit_bases+")"+toSubscript(repeatTimes);
        var locus=g.left_flank_bases.slice(-6)+"<b>"+motif+"</b>"+g.right_flank_bases.slice(0,6);
        var roiCell = g.roi_marker ? '<span class="chip roi">'+g.roi_marker+'</span>' : '<span class="lib-empty">\\u2014</span>';
        var libCell = g.libs.map(function(l){
          return '<span class="cell-mono lib-c'+((LIB_COLOR_INDEX[l.library]||0))+'">'+shortLib(l.library)+'</span>';
        }).join(' <span class="lib-empty">|</span> ');
        // repeat_times is the reference locus's own repeat count (from the
        // genome scan), not a per-library measurement -- identical across
        // libA/libB/... for the same site, so showing every library's copy
        // side by side would just repeat the same number. One value (the
        // first library, in sorted order) is enough.
        var repeatsCell = g.libs.length ? g.libs[0].repeat_times : "";
        var readsCell = g.libs.map(function(l){return l.CovReads;}).join(' | ');
        var pvalCell = g.libs.map(function(l){return fmtSci(l.pro_p);}).join(' | ');
        var statusCell = STATUS_CHIP[g.status];
        return '<tr'+(g.roi_marker?' class="roi-hit"':'')+'>'
          +'<td>'+roiCell+'</td>'
          +'<td><span class="cell-mono">'+igvRange(g)+'</span></td>'
          +'<td>'+libCell+'</td>'
          +'<td class="locus">'+locus+'</td>'
          +'<td class="num tabular">'+repeatsCell+'</td>'
          +'<td class="num tabular">'+readsCell+'</td>'
          +'<td class="num tabular">'+pvalCell+'</td>'
          +'<td>'+statusCell+'</td>'
          +'</tr>';
      }).join("");
    }
    var total=state.filtered.length, last=Math.max(0,Math.ceil(total/state.pageSize)-1);
    document.getElementById("rangeLabel").textContent = total ? ("Rows "+(start+1)+"\\u2013"+Math.min(start+state.pageSize,total)+" of "+total.toLocaleString()) : "No rows";
    document.getElementById("pageLabel").textContent="Page "+(state.page+1)+" / "+(last+1);
    document.getElementById("prevBtn").disabled = state.page<=0;
    document.getElementById("nextBtn").disabled = state.page>=last;
  }
  document.getElementById("headRow").addEventListener("click", function(e){
    var th=e.target.closest("th"); if (!th || th.dataset.nosort) return;
    var key=th.dataset.key;
    if (state.sortKey===key) state.sortDir*=-1; else { state.sortKey=key; state.sortDir=1; }
    sortRows(); state.page=0; renderTable();
  });
  document.getElementById("pvalInput").addEventListener("input", function(e){
    var v=parseFloat(e.target.value);
    state.minPval = isNaN(v) ? 0 : v;
    document.getElementById("pvalVal").textContent = e.target.value;
    applyFilters();
  });
  document.getElementById("resetBtn").addEventListener("click", function(){
    document.getElementById("pvalInput").value="0.1"; state.minPval=0.1;
    document.getElementById("pvalVal").textContent="0.1";
    applyFilters();
  });
  document.getElementById("prevBtn").addEventListener("click", function(){ if (state.page>0){ state.page--; renderTable(); } });
  document.getElementById("nextBtn").addEventListener("click", function(){ var last=Math.max(0,Math.ceil(state.filtered.length/state.pageSize)-1); if (state.page<last){ state.page++; renderTable(); } });

  applyFilters();
})();
</script>
"""
    return (html
            .replace("__TITLE__", f"MSI Report - {sample_name}")
            .replace("__HEAD_LINKS__", HEAD_LINKS)
            .replace("__CSS__", REPORT_CSS)
            .replace("__RUN_ID__", run_id)
            .replace("__RUN_INDEX_FILE__", f"{run_id}_index_run.html")
            .replace("__LIBRARY_CHIPS__", "".join(library_chips))
            .replace("__CONSENSUS_ROW__", consensus_row)
            .replace("__SAMPLE__", sample_name)
            .replace("__GENOME__", genome)
            .replace("__GENERATED__", generated_at)
            .replace("__COUNT__", str(site_count))
            .replace("__LIBCOUNT__", str(lib_count))
            .replace("__LIBPLURAL__", "y" if lib_count == 1 else "ies")
            .replace("__DATA_JS_FILE__", f"{sample_name}_index_sample.js")
            .replace("__LIB_COLOR_JSON__", lib_color_json))


def write_data_js(js_path: Path, var_name: str, data) -> None:
    """Write data as a .js file assigning it to window.<var_name>. The HTML pages
    load it via <script src>: fetch() of a local .json is blocked by browsers when
    the page is opened from file://, a <script src> is not."""
    js_path.write_text(f"window.{var_name} = {json.dumps(data, indent=2)};\n")


def compute_sample_meta(data_library_dir: Path, sample_name: str, bam_path: Path,
                        genome: str, run_id: str, library: str = None, sample_key: str = None) -> dict:
    """Compute this (sample, library)'s stats from its raw msisensor-pro output
    (read from data/, where every sample of this library sits directly, named
    with the sample as a prefix). Nothing is written: the returned dict is one
    entry of the run's "samples" list, serialized once into <run_id>_index_run.js.
    """
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # msisensor-pro's own summary output is empty in tumor-only mode, so the
    # "% unstable" MSI score is recomputed here from *_all: unstable sites
    # over evaluated sites, both counting only sites with >= MIN_READS reads.
    all_rows = parse_msisensor_tsv(data_library_dir / f"{sample_name}_all")
    unstable_region_count = len(counted_unstable_keys(all_rows))
    region_count = len(counted_site_keys(all_rows))
    pct_region_unstable = pct(unstable_region_count, region_count)

    meta = {
        "sample_name": sample_name,
        "library": library,
        "sample_key": sample_key or sample_name,
        "bam_path": str(bam_path),
        "genome": genome,
        "unstable_region_count": unstable_region_count,
        "region_count": region_count,
        "pct_region_unstable": pct_region_unstable,
        "generated_at": generated_at,
    }
    return meta


def serialize_sample_rows_to_js(viz_run_dir: Path, sample_name: str, rows: list) -> Path:
    """Serialize sample rows (already parsed and enriched from *_all) to
    <sample>_index_sample.js (window.ROWS) in /viz/runs/<run_id>/, next to the
    <sample>_index_sample.html that loads it via <script src>."""
    data_file = viz_run_dir / f"{sample_name}_index_sample.js"
    write_data_js(data_file, "ROWS", rows)
    return data_file


def load_sample_library_rows(run_dir: Path, samples: list, roi_bed_path=None) -> dict:
    """Re-read every ok sample's *_all file (every site MSIsensor-pro
    evaluated for that library, not just the ones it called unstable) into
    sample_name -> library -> rows, tagging each row with roi_marker (the
    classic MSI marker name, e.g. "BAT26", or "" if the site isn't one) from
    the panel/ROI correspondence bed. Shared by write_run_index (site-level
    stats) and write_sample_group_indexes (the pooled per-sample report), so
    the raw files are only parsed once per run.

    Reading *_all instead of *_unstable matters once libraries are pooled by
    site: a site right at the p-value/threshold border can be unstable in one
    library and not the other, and *_unstable simply omits it for the library
    that didn't cross the line -- leaving that library with no row at all for
    the site instead of its (sub-threshold) value. *_all always has a row for
    every library that evaluated the site, so callers needing "unstable only"
    filter with row_is_unstable() themselves."""
    roi_markers = parse_roi_bed(Path(roi_bed_path)) if roi_bed_path else {}
    by_sample = {}
    for s in samples:
        if s.get("failed"):
            continue
        library = s.get("library") or "-"
        lib_dir = Path(s["library"]) if s.get("library") else Path(".")
        all_path = run_dir / lib_dir / f'{s["sample_name"]}_all'
        rows = parse_msisensor_tsv(all_path)
        for row in rows:
            row["library"] = library
            row["roi_marker"] = roi_markers.get((row["chromosome"], row["location"]), "")
        by_sample.setdefault(s["sample_name"], {})[library] = rows
    return by_sample


def run_data_js_path(run_dir: Path, run_id: str) -> Path:
    return run_dir / f"{run_id}_index_run.js"


def write_run_index(run_dir: Path, run_id: str, samples: list, rows_by_sample: dict = None) -> str:
    """Write this run's <run_id>_index_run.html and its one data file,
    <run_id>_index_run.js (window.RUN_DATA): run metadata (run_id, started_at,
    every (sample, library) entry in "samples") plus the precomputed table.
    Returns the timestamp used, so the caller can record the same one in the
    top-level runs manifest."""
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    run_data = {"run_id": run_id, "started_at": generated_at, "samples": samples,
                "status_rule": sample_status_rule()}
    run_data.update(build_run_table(samples, rows_by_sample))
    write_data_js(run_data_js_path(run_dir, run_id), "RUN_DATA", run_data)
    (run_dir / f"{run_id}_index_run.html").write_text(render_run_index(run_id))
    return generated_at


def write_sample_group_indexes(run_dir: Path, run_id: str, samples: list, rows_by_sample: dict) -> None:
    """Write run_dir/<sample_name>_index_sample.html for every sample,
    associating its results across all of its libraries (or just its single
    "-" entry for a run with no library subfolders) -- directly in run_dir,
    no per-sample subfolder.
    
    Each sample's rows go to <sample>_index_sample.js, which the HTML loads
    via <script src> (Phase 2 -> Phase 3 transition)."""
    by_sample = {}
    for s in samples:
        library = s.get("library") or "-"
        by_sample.setdefault(s["sample_name"], {})[library] = s

    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    for sample_name, entries in by_sample.items():
        combined_rows = []
        for library in entries:
            combined_rows.extend(rows_by_sample.get(sample_name, {}).get(library, []))

        # Site status (MSI/MSS/MSU/NA) is decided here, not in the page:
        # every row of a locus carries the same site_status.
        ok_lib_count = sum(1 for e in entries.values() if not e.get("failed"))
        loci = {}
        for r in combined_rows:
            loci.setdefault((r["chromosome"], r["location"]), []).append(r)
        for libs_at in loci.values():
            status = site_status(libs_at, ok_lib_count)
            for r in libs_at:
                r["site_status"] = status

        # Phase 2 -> Phase 3: serialize rows to the .js the HTML loads
        serialize_sample_rows_to_js(run_dir, sample_name, combined_rows)

        html = render_sample_group_report(run_id, sample_name, generated_at, combined_rows, entries)
        (run_dir / f"{sample_name}_index_sample.html").write_text(html)


# ---------------------------------------------------------------------------
# report: publish raw outputs to data/, (re)generate viz/
# ---------------------------------------------------------------------------

def read_samples_tsv(path: Path) -> list:
    """Rows of the samples.tsv the workflow writes (see module docstring)."""
    with open(path, newline="") as fh:
        rows = [dict(zip(SAMPLES_TSV_COLUMNS, cols)) for cols in csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE) if cols]
    for row in rows:
        if len(row) != len(SAMPLES_TSV_COLUMNS):
            sys.exit(f"Error: malformed line in {path}: {row}")
    return rows


def publish_raw_outputs(row: dict, data_library_dir: Path) -> None:
    """Copy one sample's msisensor-pro output files from its WDL task directory
    into data/<run_id>/<library>/, named <sample>, <sample>_all, ...
    exactly as msisensor-pro's own -o prefix names them."""
    data_library_dir.mkdir(parents=True, exist_ok=True)
    for column, suffix in (("summary", ""), ("all", "_all"), ("dis", "_dis"), ("unstable", "_unstable")):
        if row[column]:
            shutil.copyfile(row[column], data_library_dir / f'{row["sample_name"]}{suffix}')


def report(args) -> None:
    global MIN_READS, MSI_PVAL_CUTOFF, SAMPLE_NA_FRACTION, MSI_CONSENSUS_ROI_CUTOFF
    MIN_READS, MSI_PVAL_CUTOFF, SAMPLE_NA_FRACTION = args.min_reads, args.pval_cutoff, args.sample_na_fraction
    MSI_CONSENSUS_ROI_CUTOFF = args.consensus_roi_cutoff
    run_id = args.run_id

    viz_dir = Path(args.viz_output)
    data_dir = Path(args.data_output)
    viz_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    log_path = setup_logging(run_id, Path(args.log_output))
    logger.info(f"Log file: {log_path}")

    samples = read_samples_tsv(Path(args.samples))

    # Every generated .html/.js/.json report lives under viz/runs/<run_id>/;
    # every raw msisensor-pro output file under data/<run_id>/<library>/.
    viz_run_dir = viz_dir / RUNS_DIRNAME / run_id
    data_run_dir = data_dir / run_id
    viz_run_dir.mkdir(parents=True, exist_ok=True)
    data_run_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"{len(samples)} BAM file(s) processed -> run {run_id}")

    sample_metas = []
    for row in samples:
        sample_name = row["sample_name"]
        # A run directory can hold several libraries as subfolders (e.g.
        # Run1/alignements_libA/, Run1/alignements_libB/), each with samples
        # of the same name. Mirror that subfolder under the run so results
        # never collide, and keep the subfolder name as the sample's library.
        rel_dir = Path(row["library"])
        library = rel_dir.as_posix() if rel_dir != Path(".") else None
        sample_key = (rel_dir / sample_name).as_posix()
        bam_name = Path(row["bam_path"]).name
        data_library_dir = data_run_dir / rel_dir
        publish_raw_outputs(row, data_library_dir)
        logger.info(f"{bam_name} (sample: {sample_key}) -> {data_library_dir}")
        sample_metas.append(compute_sample_meta(data_library_dir, sample_name, Path(row["bam_path"]),
                                                args.genome, run_id, library, sample_key))

    if args.roi_bed:
        logger.info(f"ROI correspondence BED: {args.roi_bed}")
    rows_by_sample = load_sample_library_rows(data_run_dir, sample_metas, args.roi_bed)
    started_at = write_run_index(viz_run_dir, run_id, sample_metas, rows_by_sample)
    write_sample_group_indexes(viz_run_dir, run_id, sample_metas, rows_by_sample)

    # Only this run's pages: the dashboard (index_msi_*.html/.js, every run)
    # is built by the index subcommand once this run is added to the reports.
    logger.info(f"Report: {viz_run_dir / (run_id + '_index_run.html')}")


def main():
    parser = argparse.ArgumentParser(description="Python steps of the msi_detection-0.1.0.wdl workflow (inside its tasks).")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("list-to-bed", help="MSIsensor-pro .list -> BED6.")
    p.add_argument("list")
    p.add_argument("output")

    p = sub.add_parser("mark-list", help="Full .list, filter=PASS for sites in the panel BED, NO otherwise.")
    p.add_argument("list")
    p.add_argument("bed")
    p.add_argument("output")

    p = sub.add_parser("report", help="Publish raw outputs to data/ and generate this run's viz/runs/<run_id>/ pages.")
    p.add_argument("--run-id", required=True)
    p.add_argument("--samples", required=True, help="samples.tsv written by the workflow.")
    p.add_argument("--genome", required=True, help="Reference fasta (only recorded in the reports).")
    p.add_argument("--data-output", required=True)
    p.add_argument("--viz-output", required=True)
    p.add_argument("--log-output", required=True)
    p.add_argument("--roi-bed", default=None, help="ROI correspondence BED (from the roi bedtools_intersect task).")
    p.add_argument("--min-reads", type=int, required=True, help="A site counts only with >= this many reads.")
    p.add_argument("--pval-cutoff", type=float, required=True, help="A site is unstable when its p-value is above this.")
    p.add_argument("--sample-na-fraction", type=float, required=True,
                   help="Sample status NA when more than this fraction of its loci listed on its sample page (a p-value in at least one library, or an ROI marker) are NA or MSU.")
    p.add_argument("--consensus-roi-cutoff", type=int, required=True,
                   help="Sample status MSI when at least this many ROI markers are unstable in consensus.")

    args = parser.parse_args()

    if args.command == "report":
        report(args)
        return

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    if args.command == "list-to-bed":
        convert_list_to_bed(Path(args.list), Path(args.output))
    elif args.command == "mark-list":
        apply_panel_filter_to_list(Path(args.list), Path(args.bed), Path(args.output))


if __name__ == "__main__":
    main()
