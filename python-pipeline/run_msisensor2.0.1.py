#!/usr/bin/env python3
"""
run_msisensor2.py

Runs MSIsensor-pro on a BAM file via an Apptainer image. The reference
fasta must be given explicitly via -g/--genome (or the "genome" key in
--config) as a path -- it is not auto-detected from the BAM, and it is
not a build keyword ("hg38", ...).

Parameters can be passed on the command line, or via a JSON config file
(--config), or a mix of both: any value present in the JSON is used as a
default, and any matching CLI flag overrides it. JSON keys match the long
option names (with dashes replaced by underscores), except for these more
readable aliases:
    "pathToRun"    -> -b/--bam (a single BAM, or a directory of BAMs for one run)
    "msisensor-pro" -> -s/--sif (the MSIsensor-pro Apptainer image)
    "bedtools"     -> --bedtools-sif (the bedtools Apptainer image)
e.g.:

    {
        "pathToRun": "/data/bams/",
        "output": "/data/results/",
        "genome": "/data/refData/hg38/Homo_sapiens_assembly38.fasta",
        "msisensor-pro": "apptainer-recipes/msisensor-pro-1.3.0.sif",
        "bedtools": "apptainer-recipes/bedtools-2.31.1.sif",
        "threads": 8,
        "recursive": true,
        "panel_bed": "/data/panels/Panel_Juno_v3_Juno.bed",
        "roi_bed": "/data/panels/msi_markers.bed"
    }

Each run of the script creates a run folder under both <output>/viz/runs/<run_id>/
and <output>/data/runs/<run_id>/, named after the input directory (e.g. "Run1";
re-running it updates that same folder in place). viz/ holds every generated
.html/.js/.json report; data/ holds every raw msisensor-pro output file -- the
two trees mirror the same runs/<run_id>/<library>/ layout, just rooted
differently. --viz-output/--data-output/--log-output ("viz_output"/
"data_output"/"log_output" in the JSON) override where each actually lives,
e.g. to put data/ on a large/slow disk and viz/+log/ on a small/fast one;
--output is then only their default parent (<output>/viz, <output>/data,
<output>/log). If -b/pathToRun holds several library subfolders (e.g.
Run1/alignements_libA/, Run1/alignements_libB/) with the same sample names,
each library's files still sit in their own folder so they never collide
(flat inside it though: no per-sample subfolder). Console output is also
logged to a per-run file under log/, named <run_id>_MSI-<version>_<timestamp>.log.
HTML pages are kept
up to date at every level:

    <output>/viz/index_msi_run.html                              - list of all runs
    <output>/viz/runs/<run_id>/<run_id>_index_run.html            - one row per sample, one column-group per library
    <output>/viz/runs/<run_id>/<sample>_index_sample.html         - that sample's result pooled across its libraries

(a run with no library subfolders still gets a <sample>_index_sample.html,
with a single library entry -- every run has exactly these two levels below
index_msi_run.html.)

Usage:
    python3 run_msisensor2.py -b tumor.bam -o /data/results/ -g /data/ref/hg38.fasta
    python3 run_msisensor2.py -b tumor.bam -o /data/results/ -g /data/ref/hg38.fasta -t 4
    python3 run_msisensor2.py --config params.json

    # -b can also be a directory: every *.bam it contains is processed
    # (--recursive to descend into subdirectories).
    python3 run_msisensor2.py -b /data/bams/ -o /data/results/ -g /data/ref/hg38.fasta
    python3 run_msisensor2.py -b /data/bams/ -o /data/results/ -g /data/ref/hg38.fasta --recursive
"""

import argparse
import json
import logging
import shlex
import re
import subprocess
import sys
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from urllib.parse import quote

SCRIPT_DIR = Path(__file__).resolve().parent

RUNS_DIRNAME = "runs"
VIZ_DIRNAME = "viz"   # <output>/viz/    -- every generated .html/.js/.json (the reports)
DATA_DIRNAME = "data"  # <output>/data/  -- every raw msisensor-pro output file
LOG_DIRNAME = "log"   # <output>/log/    -- per-run log files

# Hardcoded (not queried from the container): update when apptainer-recipes/
# ships a different msisensor-pro build. Used in the log filename/content.
TOOL_VERSION = "1.3.0"

logger = logging.getLogger("run_msisensor2")

# Custom level for external commands actually executed (between INFO and
# WARNING, so it shows whenever INFO does): "[COMMAND] apptainer exec ...".
COMMAND = 25
logging.addLevelName(COMMAND, "COMMAND")


def log_command(cmd: list, stdout_path: Path = None) -> None:
    """Log an external command exactly as it is run, shell-quoted so it can be
    copy-pasted and re-run as-is (including its stdout redirection, if any)."""
    line = shlex.join(str(c) for c in cmd)
    if stdout_path is not None:
        line += f" > {shlex.quote(str(stdout_path))}"
    logger.log(COMMAND, line)


def setup_logging(run_id: str, log_dir: Path) -> Path:
    """Configure the module logger: console + a per-run file under log_dir.

    The tool version is baked into the log filename and logged as the first
    line, so a log file is self-describing.
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = log_dir / f"{run_id}_MSI-{TOOL_VERSION}_{timestamp}.log"

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")

    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    logger.info(f"msisensor-pro {TOOL_VERSION}")
    return log_path


def run_quiet(cmd: list) -> subprocess.CompletedProcess:
    """Run a subprocess without flooding the console: msisensor-pro's own
    (very verbose, per-chromosome/per-window) stdout is always discarded.
    stderr is captured and only logged (to console + log file) if the
    command actually fails, so real errors aren't lost."""
    result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0 and result.stderr:
        logger.error(f"Command failed (code {result.returncode}), stderr:\n{result.stderr}")
    return result


def fatal(msg: str) -> None:
    """Log an error and exit(1). Use instead of sys.exit() once logging is set up."""
    logger.error(msg)
    sys.exit(1)

REF_DATA_DIR = Path("/media/jp/SSD2/CHU/refData")
# static/  : user-provided, not touched by the pipeline (panel BEDs, custom ROI BEDs, ...).
# generated/: every intermediate reference file the pipeline derives on its own
# (scan output, converted BEDs, panel/ROI restrictions, ...).
STATIC_DIR = REF_DATA_DIR / "MSI" / "static"
GENERATED_DIR = REF_DATA_DIR / "MSI" / "generated"

# Every apptainer invocation binds this host disk at a *different* path
# inside the container (not an identity mount), so any path passed as a CLI
# argument to a containerized tool -- everything the tool reads or writes,
# e.g. -a/-b/-d/-o/-t -- must be translated with to_container_path() first.
# The .sif image path itself is the one exception: apptainer always resolves
# it on the host side, never through the bind, so it's passed unchanged.
HOST_MOUNT_ROOT = Path("/media/jp/SSD2")
CONTAINER_MOUNT_ROOT = Path("/InsideContainer")


def to_container_path(host_path) -> str:
    """Translate a host path under HOST_MOUNT_ROOT to its equivalent inside
    the apptainer container (bound at CONTAINER_MOUNT_ROOT)."""
    host_path = Path(host_path)
    try:
        return str(CONTAINER_MOUNT_ROOT / host_path.relative_to(HOST_MOUNT_ROOT))
    except ValueError:
        fatal(f"{host_path} is outside {HOST_MOUNT_ROOT}, the only disk bound into the container "
              f"(as {CONTAINER_MOUNT_ROOT}) -- move it under {HOST_MOUNT_ROOT} or extend the --bind.")

DEFAULT_SIF = SCRIPT_DIR.parent / "wdl" / "apptainer-recipes" / "msisensor-pro-1.3.0.sif"
DEFAULT_BEDTOOLS_SIF = SCRIPT_DIR.parent / "wdl" / "apptainer-recipes" / "bedtools-2.31.1.sif"



@contextmanager
def atomic_output(path: Path):
    """Yield a temporary ".<name>.part" path next to path, renamed to path only
    once the with-block completes. An interrupted write (Ctrl-C, crash, fatal())
    thus never leaves a truncated file under the final name, which
    is_up_to_date() would otherwise happily reuse on every later run."""
    tmp = path.with_name(f".{path.name}.part")
    try:
        yield tmp
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def is_up_to_date(output_path: Path, *input_paths: Path) -> bool:
    """True if output_path exists and is at least as recent as every input:
    a derived reference file is then reused as-is instead of being rebuilt,
    and editing an input (e.g. the panel/ROI BED) is enough to trigger a rebuild."""
    if not output_path.exists():
        return False
    out_mtime = output_path.stat().st_mtime
    return all(out_mtime >= p.stat().st_mtime for p in input_paths if p.exists())


def convert_list_to_bed(list_path: Path) -> Path:
    """Convert a MSIsensor-pro *.msisensor_pro.list reference file into a BED6
    file (chrom, start, end, name, score=0, strand=+), generated once next to
    the list.

    "location" in the .list is the 1-based start of the repeat; BED is
    0-based half-open, so start = location - 1 and end = start +
    repeat_unit_length * repeat_times.
    """
    bed_path = list_path.with_suffix(".bed")
    if bed_path.exists():
        return bed_path

    logger.info(f"Generating BED file from microsatellite list: {bed_path}")
    count = 0
    with atomic_output(bed_path) as tmp_path, open(list_path) as fh_in, open(tmp_path, "w") as fh_out:
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


def convert_panel_to_probe_bed(panel_bed_path: Path) -> Path:
    """Convert a Juno-style panel BED into a plain BED6: strip the UCSC
    "track ..." header line, and use each region's probe/target span
    (thickStart/thickEnd, columns 7-8) as start/end instead of the full
    amplicon (columns 2-3). Name is kept as-is; score=0, strand=+.

    Treated as static, not generated: once prepared for a panel it's a stable
    input in its own right (like the raw panel BED it comes from), reused
    as-is by every run rather than rebuilt from it -- so the output always
    goes to STATIC_DIR regardless of where the input panel BED lives.
    """
    STATIC_DIR.mkdir(parents=True, exist_ok=True)
    output_path = STATIC_DIR / f"{panel_bed_path.stem}.amplicons.bed"
    count = 0
    with open(panel_bed_path) as src, open(output_path, "w") as dst:
        for line in src:
            if line.startswith("track") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            chrom = cols[0]
            name = cols[3]
            start, end = cols[6], cols[7]
            dst.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t+\n")
            count += 1
    logger.info(f"Panel probe BED written: {output_path} ({count:,} regions)")
    return output_path


def _msi_bed_uses_chr_prefix(msi_bed_path: Path) -> bool:
    """Peek at the genome-wide microsatellite BED's first chromosome field to
    tell whether this reference build names contigs "chr1" (UCSC-style, e.g.
    hg19/hg38) or "1" (Ensembl-style, e.g. GRCh37.75) -- so the panel/ROI BED
    can be aligned to it before bedtools intersect, which matches chrom names
    literally and silently returns zero overlaps on a naming mismatch."""
    with open(msi_bed_path) as fh:
        for line in fh:
            if line.strip():
                return line.split("\t", 1)[0].startswith("chr")
    return True


def generate_panel_correspondence_bed(msi_bed_path: Path, panel_bed_path: Path, bedtools_sif: Path,
                                       output_suffix: str = "panel", append_position: bool = False) -> Path:
    """Intersect the genome-wide microsatellite BED with a panel/target or ROI BED
    via bedtools (in its own container), then collapse the overlaps with awk into
    <msi_bed_path stem>.<output_suffix>.bed (BED6): chrom, start, end,
    "<microsatellite name>#<region name>", score=0, strand=+.

    Used for both --panel-bed (sequencing panel amplicons) and --roi-bed
    (e.g. the classic BAT26/NR21/... MSI marker regions) -- same mechanism,
    different input BED and output_suffix.

    append_position appends the matched panel_bed_path region's own
    "#<chrom>:<start>-<end>" (chrom without a "chr" prefix, regardless of the
    BEDs' own naming convention) after the region name -- set for the
    --panel-bed step (where panel_bed_path is the probe/amplicon BED) so the
    probe's genomic position travels with its name. The chained --roi-bed step
    then just appends the marker name on top (name is already
    "<microsatellite>#<probe>#<position>"), giving a final
    "<microsatellite>#<probe>#<position>#<marker>" -- position before the
    marker name, not after, so the marker (the clinically meaningful id) stays
    the last, most-visible segment.

    Generated once and reused: only rebuilt when missing or older than
    msi_bed_path or panel_bed_path (see is_up_to_date).
    """
    if not bedtools_sif.exists():
        fatal(f"Apptainer image not found: {bedtools_sif}")
    if not panel_bed_path.exists():
        fatal(f"BED not found: {panel_bed_path}")

    output_path = msi_bed_path.parent / f"{msi_bed_path.stem}.{output_suffix}.bed"
    if is_up_to_date(output_path, msi_bed_path, panel_bed_path):
        return output_path
    normalized_panel_path = msi_bed_path.parent / f".{panel_bed_path.stem}.chrnorm.tmp"
    intersect_path = msi_bed_path.parent / f".{panel_bed_path.stem}.intersect.tmp"

    # panel_bed_path is always an already-cleaned BED here (e.g. the static
    # *.amplicons.bed convert_panel_to_probe_bed() produces, or a plain ROI
    # BED) -- no UCSC "track ..." header to strip. Just align chrom naming
    # ("chr1" vs "1") to the microsatellite BED's own convention, otherwise
    # bedtools intersect silently matches nothing.
    use_chr_prefix = _msi_bed_uses_chr_prefix(msi_bed_path)
    renamed = 0
    with open(panel_bed_path) as src, open(normalized_panel_path, "w") as dst:
        for line in src:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            has_prefix = cols[0].startswith("chr")
            if use_chr_prefix and not has_prefix:
                cols[0] = "chr" + cols[0]
                renamed += 1
            elif not use_chr_prefix and has_prefix:
                cols[0] = cols[0][3:]
                renamed += 1
            dst.write("\t".join(cols) + "\n")
    if renamed:
        logger.info(
            f"Normalized chromosome naming ({'added' if use_chr_prefix else 'stripped'} 'chr' "
            f"prefix) on {renamed:,} line(s) of {panel_bed_path.name} to match {msi_bed_path.name}."
        )

    cmd = [
        "apptainer", "exec",
        "--bind", f"{HOST_MOUNT_ROOT}:{CONTAINER_MOUNT_ROOT}",
        str(bedtools_sif),
        "bedtools", "intersect",
        "-a", to_container_path(msi_bed_path),
        "-b", to_container_path(normalized_panel_path),
        "-wa", "-wb",
        "-f", "1.0",  # the microsatellite site must be fully contained in the panel region
    ]
    log_command(cmd, intersect_path)
    with open(intersect_path, "w") as out:
        result = subprocess.run(cmd, stdout=out)
    if result.returncode != 0:
        normalized_panel_path.unlink(missing_ok=True)
        intersect_path.unlink(missing_ok=True)
        fatal(f"bedtools intersect failed (code {result.returncode}).")

    logger.info(f"Post-processing overlaps with awk: {output_path}")
    # One row per overlap: name becomes "<microsatellite>#<panel region>", with a
    # fixed score/strand (0, +) rather than a per-site overlap count. $7-$9 are
    # the matched panel_bed_path region's own chrom/start/end from this same
    # intersect record, so pairing them with $10 (its name) needs no extra
    # lookup or filtering even when a site overlaps more than one region --
    # each overlap is already its own output row.
    if append_position:
        awk_script = (
            'BEGIN{OFS="\\t"} {chrom=$7; sub(/^chr/,"",chrom); '
            'print $1,$2,$3,$4"#"$10"#"chrom":"$8"-"$9,0,"+"}'
        )
    else:
        awk_script = 'BEGIN{OFS="\\t"} {print $1,$2,$3,$4"#"$10,0,"+"}'
    awk_cmd = ["awk", awk_script, str(intersect_path)]
    with atomic_output(output_path) as tmp_path:
        log_command(awk_cmd, tmp_path)
        with open(tmp_path, "w") as out:
            result = subprocess.run(awk_cmd, stdout=out)
        normalized_panel_path.unlink(missing_ok=True)
        intersect_path.unlink(missing_ok=True)
        if result.returncode != 0:
            fatal(f"awk post-processing failed (code {result.returncode}).")

    logger.info(f"Panel correspondence BED written: {output_path}")
    return output_path


def filter_list_to_panel(list_path: Path, panel_bed_path: Path, output_suffix: str = "panel") -> Path:
    """Restrict the full genome-wide MSIsensor-pro .list to just the sites
    selected in <output_suffix>.bed (built by generate_panel_correspondence_bed
    from this same list), keeping the original .list column format untouched.

    Matches rows by (chromosome, start), where start = location - 1 -- the same
    site is exactly the BED interval that generate_panel_correspondence_bed
    computed for it, so no fuzzy matching is needed.

    Used only for the ROI/marker correspondence list (a small, deliberately
    restricted subset kept for reporting) -- see apply_panel_filter_to_list
    for the panel-wide list, which keeps every row instead of dropping any.
    """
    output_path = list_path.parent / f"{list_path.stem}.{output_suffix}.list"
    if is_up_to_date(output_path, list_path, panel_bed_path):
        return output_path

    keys = set()
    with open(panel_bed_path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2:
                keys.add((cols[0], cols[1]))

    count = 0
    with atomic_output(output_path) as tmp_path, open(list_path) as fh_in, open(tmp_path, "w") as fh_out:
        header = fh_in.readline()
        fh_out.write(header)
        idx = {name: i for i, name in enumerate(header.rstrip("\n").split("\t"))}
        for line in fh_in:
            cols = line.rstrip("\n").split("\t")
            if not cols or not cols[0]:
                continue
            start = str(int(cols[idx["location"]]) - 1)
            if (cols[idx["chromosome"]], start) in keys:
                fh_out.write(line if line.endswith("\n") else line + "\n")
                count += 1

    logger.info(f"Panel-restricted list written: {output_path} ({count:,} sites)")
    return output_path


def apply_panel_filter_to_list(list_path: Path, panel_bed_path: Path, output_suffix: str = "panel") -> Path:
    """Write <list_path stem>.<output_suffix>.list with every row of the full
    genome-wide list kept (same row count/order as list_path), but the
    "filter" column set to "PASS" for sites selected in <output_suffix>.bed
    (built by generate_panel_correspondence_bed from this same list) and
    "NO" for every other site.

    Unlike filter_list_to_panel (which drops non-matching rows entirely),
    this keeps the panel restriction visible/inspectable directly in the
    list -- every genome-wide site is still there, just tagged.

    Matches rows by (chromosome, start), where start = location - 1 -- the
    same site is exactly the BED interval that generate_panel_correspondence_bed
    computed for it, so no fuzzy matching is needed.
    """
    output_path = list_path.parent / f"{list_path.stem}.{output_suffix}.list"
    if is_up_to_date(output_path, list_path, panel_bed_path):
        return output_path

    keys = set()
    with open(panel_bed_path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2:
                keys.add((cols[0], cols[1]))

    pass_count = 0
    total = 0
    with atomic_output(output_path) as tmp_path, open(list_path) as fh_in, open(tmp_path, "w") as fh_out:
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


def ensure_msisensor_pro_list(genome: str, sif: Path, panel_bed: str = None,
                               roi_bed: str = None, bedtools_sif: Path = None) -> tuple:
    """Return (analysis_list_path, roi_correspondence_bed_path) for MSIsensor-pro,
    generating the -d list if missing.

    genome is the path to the reference fasta (not a build keyword). If
    panel_bed is given, the panel-restricted list is returned instead of the
    full genome-wide one -- analysis itself is scoped to the panel, not just
    the reporting; no separate "which list to analyze" parameter is needed,
    it's derived from the general list automatically.

    If roi_bed is also given, the ROI is cross-referenced against the
    panel-restricted sites (not the whole genome) -- the classic markers are
    expected to already sit inside the panel's own amplicons -- yielding
    <list>.panel.roi.bed/.list. Without a panel, it falls back to the full
    genome-wide BED, yielding plain <list>.roi.bed/.list. roi_correspondence_bed_path
    is None when roi_bed isn't given.
    """
    fasta_path = Path(genome)

    if not fasta_path.exists():
        fatal(f"Reference fasta not found: {fasta_path}")

    msi_dir = GENERATED_DIR
    msi_dir.mkdir(parents=True, exist_ok=True)
    list_path = msi_dir / f"{fasta_path.name}.msisensor_pro.list"
    analysis_list_path = list_path
    roi_correspondence_path = None

    def build_correspondence_beds(bed_path: Path) -> None:
        nonlocal analysis_list_path, roi_correspondence_path
        panel_correspondence = None
        if panel_bed:
            panel_correspondence = generate_panel_correspondence_bed(
                bed_path, Path(panel_bed), bedtools_sif, output_suffix="panel", append_position=True)
            analysis_list_path = apply_panel_filter_to_list(list_path, panel_correspondence, output_suffix="panel")
        if roi_bed:
            roi_source_bed = panel_correspondence or bed_path
            # generate_panel_correspondence_bed names its output after
            # roi_source_bed's own stem, which already ends in ".panel" when
            # chained from it -- so "roi" alone is enough there ("panel" +
            # "roi" -> "panel.roi.bed"). filter_list_to_panel always starts
            # from the one master list_path, which never carries ".panel" in
            # its name, so it needs the combined suffix spelled out.
            list_suffix = "panel.roi" if panel_correspondence else "roi"
            roi_correspondence_path = generate_panel_correspondence_bed(
                roi_source_bed, Path(roi_bed), bedtools_sif, output_suffix="roi")
            filter_list_to_panel(list_path, roi_correspondence_path, output_suffix=list_suffix)

    if list_path.exists():
        bed_path = convert_list_to_bed(list_path)
        build_correspondence_beds(bed_path)
        return analysis_list_path, roi_correspondence_path

    logger.info(f"Generating MSIsensor-pro reference file (scan): {list_path}")
    with atomic_output(list_path) as tmp_path:
        cmd = [
            "apptainer", "exec",
            "--bind", f"{HOST_MOUNT_ROOT}:{CONTAINER_MOUNT_ROOT}",
            str(sif),
            "msisensor-pro", "scan",
            "-d", to_container_path(fasta_path),
            "-o", to_container_path(tmp_path),
        ]
        logger.info("msisensor-pro scan output is very verbose; silencing it (not printed or logged).")
        log_command(cmd)
        result = run_quiet(cmd)
        if result.returncode != 0:
            fatal(f"msisensor-pro scan failed (code {result.returncode}).")
    bed_path = convert_list_to_bed(list_path)
    build_correspondence_beds(bed_path)
    return analysis_list_path, roi_correspondence_path


def run_msisensor_pro(bam_path: Path, output_prefix: str, sif: Path, list_path: Path, threads: int,
                       coverage_threshold: int) -> None:
    cmd = [
        "apptainer", "exec",
        "--bind", f"{HOST_MOUNT_ROOT}:{CONTAINER_MOUNT_ROOT}",
        str(sif),
        "msisensor-pro", "pro",
        "-d", to_container_path(list_path),
        "-t", to_container_path(bam_path),
        "-o", to_container_path(output_prefix),
        "-b", str(threads),
        "-c", str(coverage_threshold),
    ]
    logger.info("msisensor-pro pro output is very verbose; silencing it (not printed or logged).")
    log_command(cmd)
    result = run_quiet(cmd)
    if result.returncode != 0:
        fatal(f"msisensor-pro failed (code {result.returncode}).")
    logger.info(f"MSIsensor-pro finished. Results: {output_prefix}*")


def find_bam_files(bam_arg: Path, recursive: bool) -> list:
    """Resolve the -b argument: a single BAM, or every *.bam in a directory (glob)."""
    if bam_arg.is_dir():
        pattern = "**/*.bam" if recursive else "*.bam"
        bam_files = sorted(bam_arg.glob(pattern))
        if not bam_files:
            fatal(f"No .bam file found in {bam_arg} (pattern: {pattern}).")
        return bam_files
    if not bam_arg.exists():
        fatal(f"BAM file not found: {bam_arg}")
    return [bam_arg]


def load_config(config_path: Path) -> dict:
    """Load parameter defaults from a JSON file (keys match the CLI dest names)."""
    if not config_path.exists():
        sys.exit(f"Error: config file not found: {config_path}")
    with open(config_path) as fh:
        try:
            return json.load(fh)
        except json.JSONDecodeError as e:
            sys.exit(f"Error: invalid JSON in {config_path}: {e}")


def prepare_reference(args) -> tuple:
    """Resolve (sif_path, analysis_list_path, roi_correspondence_bed_path) once
    per run, before any BAM is processed. Every derived reference file (scan
    .list, .bed, panel/ROI correspondence .bed/.list) is generated only if
    missing or outdated, then reused by every sample and every later run."""
    sif_path = Path(args.sif)
    if not sif_path.exists():
        fatal(f"Apptainer image not found: {sif_path}")

    logger.info(f"Reference genome: {args.genome}")
    roi_correspondence_path = None
    if args.pro_list:
        list_path = Path(args.pro_list)
        if not list_path.exists():
            fatal(f"MSIsensor-pro -d file not found: {list_path}")
    else:
        list_path, roi_correspondence_path = ensure_msisensor_pro_list(
            args.genome, sif_path, args.panel_bed, args.roi_bed, Path(args.bedtools_sif))
    logger.info(f"MSIsensor-pro reference list: {list_path}")
    if roi_correspondence_path:
        logger.info(f"ROI correspondence BED: {roi_correspondence_path}")
    return sif_path, list_path, roi_correspondence_path


def process_bam(bam_path: Path, output_dir: Path, sample_name: str, sif_path: Path, list_path: Path, args) -> None:
    """Run msisensor-pro on a single BAM against the already-prepared reference list."""
    output_prefix = str(output_dir / sample_name)
    run_msisensor_pro(bam_path, output_prefix, sif_path, list_path, args.threads, args.coverage_threshold)


# ---------------------------------------------------------------------------
# HTML reporting: RUN index -> SAMPLE index -> per-sample RESULT report
# ---------------------------------------------------------------------------

UNSTABLE_COLUMNS = [
    "chromosome", "location", "left_flank_bases", "repeat_times",
    "repeat_unit_bases", "right_flank_bases", "pro_p", "pro_q",
    "CovReads", "threshold",
]


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
                "threshold": float(get(cols, "threshold") or 0),
            })
    return rows


def row_is_unstable(row: dict) -> bool:
    """MSIsensor-pro's own unstable/stable call for one *_all row: a site
    goes into *_unstable exactly when pro_p > threshold (verified against real
    output -- the two files never disagree). Used to recompute that same call
    ourselves once rows come from *_all rather than the pre-filtered
    *_unstable file."""
    return row["pro_p"] > row["threshold"]


def count_data_lines(path: Path) -> int:
    """Count non-empty data lines with CovReads > 0 (excluding the header) in a
    MSIsensor-pro TSV output, e.g. *_all -- used to recompute the % of unstable
    sites ourselves, since the official summary file is empty in tumor-only mode.
    Lines with CovReads=0 are excluded from the count."""
    if not path.exists():
        return 0
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        cov_idx = idx.get("CovReads")
        if cov_idx is None:
            # Fallback: if CovReads column doesn't exist, count all lines
            return sum(1 for line in fh if line.strip())
        count = 0
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            try:
                cov_reads = int(cols[cov_idx]) if cov_idx < len(cols) else 0
                if cov_reads > 0:
                    count += 1
            except (ValueError, IndexError):
                # If parsing fails, count the line anyway
                count += 1
        return count


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
.index-table td{padding:11px 12px;border-bottom:1px solid var(--border);}
.index-table tr:last-child td{border-bottom:none;}
.index-table tr:hover td{background:var(--surface-alt);}
footer{color:var(--text-muted);font-size:.78rem;text-align:center;padding-top:6px;}
"""

FAVICON_PATH = Path(__file__).parent / "favicon-msi.svg"


def favicon_link() -> str:
    """<link rel="icon"> embedding favicon-msi.svg (next to this script) as a
    data URI, so every page carries its icon wherever viz/ ends up and even
    opened on its own. Empty if the file is missing."""
    if not FAVICON_PATH.exists():
        return ""
    return f'<link rel="icon" type="image/svg+xml" href="data:image/svg+xml,{quote(FAVICON_PATH.read_text().strip())}">'


# Shared <head> links of every page: web fonts + favicon.
HEAD_LINKS = ('<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;500;600&display=swap">'
              + favicon_link())


def extract_library_suffix(lib: str) -> str:
    """"alignements_libA" -> "libA": strips everything before the trailing
    "lib<x>" segment (e.g. the "alignements_" folder-naming prefix), falling
    back to the original string if it doesn't end in "lib<x>"."""
    m = re.search(r"(lib[A-Za-z0-9]+)$", lib, re.IGNORECASE)
    return m.group(1) if m else lib


MSI_CONSENSUS_ROI_CUTOFF = 3  # consensus ROI unstable sites >= this -> MSI


def build_run_table(samples: list, rows_by_sample: dict = None) -> dict:
    """Precompute (Phase 2) everything <run_id>_index_run.html displays: one row
    per sample, one column-group per library (e.g. alignements_libA,
    alignements_libB) side by side instead of a row per (sample, library) pair.

    Per library: unstable-site count over the whole analysed region
    (unstable_region_count), and how many of those land on a classic MSI marker
    (unstable_roi_count); a failed library is just {"failed": true}. "% unstable"
    (pct_region_unstable) is the simple average of each ok library's own one. When a sample has 2+ ok libraries, consensus
    counts give how many unstable sites (total, and ROI-only) every one of its
    libraries shares, and "call" is MSI when the consensus ROI count reaches
    MSI_CONSENSUS_ROI_CUTOFF."""
    rows_by_sample = rows_by_sample or {}
    libraries = sorted({s.get("library") or "-" for s in samples})

    by_sample = {}
    for s in samples:
        by_sample.setdefault(s["sample_name"], {})[s.get("library") or "-"] = s

    # rows come from *_all (every evaluated site, not just the unstable ones),
    # so these need their own unstable filter.
    def site_keys(rows):
        return {(r["chromosome"], r["location"]) for r in rows if row_is_unstable(r)}

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
            "pct_region_unstable": round(sum(ok_pcts) / len(ok_pcts), 1) if ok_pcts else None,
            "libraries": libs,
            "consensus_total": None,
            "consensus_roi": None,
            "call": None,
        }
        if len(ok_libs) >= 2:
            row["consensus_total"] = len(set.intersection(*(site_keys(lib_rows.get(lib, [])) for lib in ok_libs)))
            row["consensus_roi"] = len(set.intersection(*(roi_site_keys(lib_rows.get(lib, [])) for lib in ok_libs)))
            if not any(roi_covered(lib_rows.get(lib, [])) for lib in ok_libs):
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
  var sub = '';
  for (var i = 0; i <= D.libraries.length; i++) sub += '<th class="num">Total</th><th class="num">ROI</th>';
  var thead = '<tr><th rowspan="2">Sample</th><th rowspan="2">% unstable</th>' + groups
        + '<th class="num" colspan="2">UNSTABLE_CONSENSUS</th><th rowspan="2">Status</th></tr><tr>' + sub + '</tr>';
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
      + num(pct(r.pct_region_unstable)) + cells + num(n(r.consensus_total)) + num(n(r.consensus_roi))
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

    library_chips = []
    for library in sorted(entries_by_library):
        entry = entries_by_library[library]
        lib_label = extract_library_suffix(library)
        if not entry.get("failed"):
            lib_pct = f'{entry["pct_region_unstable"]}%' if entry.get("pct_region_unstable") is not None else "-"
            library_chips.append(
                f'<span class="meta-chip">{lib_label}: <strong>{lib_pct} unstable</strong> '
                f'({entry["unstable_region_count"]}/{entry.get("region_count", 0)})</span> '
            )
        else:
            library_chips.append(
                f'<span class="meta-chip">{lib_label}: <span class="status-pill failed">failed</span></span>'
            )

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
          roi_marker:r.roi_marker, libs:[]
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

  // MSI: every library clears p > 0.1 (all agree it's unstable). MSS: no
  // library clears it (all agree it's stable -- only reachable here if the
  // adjustable display filter below let such a row through, e.g. after
  // lowering Min p-value; applyFilters normally hides it). MSU: libraries
  // disagree -- some clear it, some don't -- an unresolved call, like the
  // chr4:1806012 case (libA 0.0965, libB 0.1037: one either side of 0.1).
  // Fixed at 0.1 regardless of state.minPval, which only controls which
  // rows are shown, not this call. Only libraries that actually covered the
  // site count (CovReads > 0; msisensor-pro reports pro_p as NaN at 0 reads);
  // NA when none did.
  var MSI_PVAL_CUTOFF = 0.1;
  function covered(l){ return l.CovReads > 0 && !isNaN(l.pro_p); }
  function siteStatus(g){
    var libs = g.libs.filter(covered);
    if (!libs.length) return "NA";
    var passing = libs.filter(function(l){ return l.pro_p > MSI_PVAL_CUTOFF; }).length;
    if (passing === libs.length) return "MSI";
    if (passing === 0) return "MSS";
    return "MSU";
  }
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
    if (key==="status") return STATUS_RANK[siteStatus(g)];
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
        var statusCell = STATUS_CHIP[siteStatus(g)];
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
            .replace("__SAMPLE__", sample_name)
            .replace("__GENOME__", genome)
            .replace("__GENERATED__", generated_at)
            .replace("__COUNT__", str(site_count))
            .replace("__LIBCOUNT__", str(lib_count))
            .replace("__LIBPLURAL__", "y" if lib_count == 1 else "ies")
            .replace("__DATA_JS_FILE__", f"{sample_name}_index_sample.js")
            .replace("__LIB_COLOR_JSON__", lib_color_json))


RUNS_DATA_FILENAME = "index_msi_run.js"
OUTPUT_INDEX_FILENAME = "index_msi_run.html"
SAMPLES_DATA_FILENAME = "index_msi_sample.js"


def write_data_js(js_path: Path, var_name: str, data) -> None:
    """Write data as a .js file assigning it to window.<var_name>. The HTML pages
    load it via <script src>: fetch() of a local .json is blocked by browsers when
    the page is opened from file://, a <script src> is not."""
    js_path.write_text(f"window.{var_name} = {json.dumps(data, indent=2)};\n")


def read_data_js(js_path: Path):
    """Inverse of write_data_js: the data a "window.<var> = <json>;" file holds."""
    text = js_path.read_text()
    return json.loads(text.split("=", 1)[1].strip().rstrip(";"))


def render_output_index_shell() -> str:
    """Static page shell. Loads run summaries from index_msi_run.js (the
    window.RUNS_DATA) via <script src>."""
    html = """<title>MSI Dashboard - Run View</title>
__HEAD_LINKS__
<style>__CSS__
.nav-links{{display:flex;gap:12px;margin-bottom:20px;}}
.nav-links a{{display:inline-flex;align-items:center;gap:6px;padding:8px 16px;background:var(--accent-soft);
  color:var(--accent);border-radius:6px;font-weight:500;transition:all 0.2s;text-decoration:none;}}
.nav-links a:hover{{background:var(--accent);color:white;}}
</style>
<div class="wrap">
  <header>
    <h1>MSI Dashboard - Run View</h1>
    <p>Every run of the MSI-Sensor-Pro pipeline, most recent first.</p>
    <div class="nav-links">
      <a href="./index_msi_sample.html">All Samples →</a>
    </div>
  </header>
  <div class="panel" style="padding:10px 14px;">
    <input type="text" id="searchInput" placeholder="Filter runs..."
      style="width:100%;border:1px solid var(--border);border-radius:8px;background:var(--bg);color:var(--text);padding:8px 10px;font:inherit;">
  </div>
  <div class="table-panel">
    <div class="table-scroll">
      <table class="index-table">
        <thead><tr><th>Run</th><th>Started</th><th class="num"># Samples</th></tr></thead>
        <tbody id="tbody"></tbody>
      </table>
    </div>
  </div>
  <footer>MSIsensor-pro pipeline report</footer>
</div>
<script src="./__RUNS_DATA_FILE__"></script>
<script>
(function(){
  "use strict";
  var runs = window.RUNS_DATA || [];
  function render(filterText){
        var f = (filterText||"").trim().toLowerCase();
        var rows = runs.filter(function(r){ return !f || r.run_id.toLowerCase().indexOf(f)!==-1; });
        var body = document.getElementById("tbody");
        if (!rows.length){
          body.innerHTML = '<tr class="empty-row"><td colspan="3">No run yet.</td></tr>';
          return;
        }
        body.innerHTML = rows.map(function(r){
          return '<tr><td><a href="./__RUNS_DIR__/'+r.run_id+'/'+r.run_id+'_index_run.html">'+r.run_id+'</a></td>'
            +'<td>'+r.started_at+'</td>'
            +'<td class="num tabular">'+r.samples_count+'</td></tr>';
        }).join("");
  }
  document.getElementById("searchInput").addEventListener("input", function(e){ render(e.target.value); });
  render("");
})();
</script>
"""
    return (html
            .replace("__HEAD_LINKS__", HEAD_LINKS)
            .replace("__CSS__", REPORT_CSS)
            .replace("__RUNS_DATA_FILE__", RUNS_DATA_FILENAME)
            .replace("__RUNS_DIR__", RUNS_DIRNAME))


def compute_sample_meta(data_library_dir: Path, sample_name: str, bam_path: Path,
                        genome: str, run_id: str, library: str = None, sample_key: str = None) -> dict:
    """Compute this (sample, library)'s stats from its raw msisensor-pro output
    (read from data/, where every sample of this library sits directly, named
    with the sample as a prefix). Nothing is written: the returned dict is one
    entry of the run's "samples" list, serialized once into <run_id>_index_run.js.
    """
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    unstable_path = data_library_dir / f"{sample_name}_unstable"
    unstable_region_count = len(parse_msisensor_tsv(unstable_path))

    # msisensor-pro's own summary output is empty in tumor-only mode, so the
    # "% unstable" MSI score is recomputed here from the raw line counts
    # instead: being in *_unstable at all is msisensor-pro's own confidence
    # call, no q-value bucketing of our own needed.
    region_count = count_data_lines(data_library_dir / f"{sample_name}_all")
    pct_region_unstable = round(100 * unstable_region_count / region_count, 1) if region_count else None

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
    run_data = {"run_id": run_id, "started_at": generated_at, "samples": samples}
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

        # Phase 2 -> Phase 3: serialize rows to the .js the HTML loads
        serialize_sample_rows_to_js(run_dir, sample_name, combined_rows)

        html = render_sample_group_report(run_id, sample_name, generated_at, combined_rows, entries)
        (run_dir / f"{sample_name}_index_sample.html").write_text(html)


def _run_summary(run_id: str, started_at: str, samples: list) -> dict:
    """One index_msi_run.js entry for a run, from its sample_metas (or the
    "samples" list stored in that run's own <run_id>_index_run.js). Counts distinct
    samples, not (sample, library) entries: a sample with 2 libraries must
    still count as 1 sample here, not 2."""
    return {
        "run_id": run_id,
        "started_at": started_at,
        "samples_count": len({s["sample_name"] for s in samples}),
    }


def rebuild_runs_manifest(viz_dir: Path) -> list:
    """Rebuild index_msi_run.js from scratch by scanning every run's own
    viz/runs/<run_id>/<run_id>_index_run.js, instead of relying on the incremental
    history append_run_to_manifest() normally keeps by editing
    index_msi_run.js in place.

    Used to repair/regenerate index_msi_run.js (e.g. after deleting it to test)
    without re-running every past run: index_msi_run.html itself never needs
    touching, it only ever reads this file via window.RUNS_DATA."""
    runs = []
    runs_dir = viz_dir / RUNS_DIRNAME
    if runs_dir.is_dir():
        for run_dir in sorted(d for d in runs_dir.iterdir() if d.is_dir()):
            js_path = run_data_js_path(run_dir, run_dir.name)
            if not js_path.exists():
                continue
            try:
                run_data = read_data_js(js_path)
            except (json.JSONDecodeError, IndexError):
                continue
            runs.append(_run_summary(run_data["run_id"], run_data.get("started_at", ""),
                                     run_data.get("samples", [])))
    runs.sort(key=lambda r: r["started_at"], reverse=True)
    write_data_js(viz_dir / RUNS_DATA_FILENAME, "RUNS_DATA", runs)
    return runs


def append_run_to_manifest(viz_dir: Path, run_id: str, started_at: str, samples: list) -> None:
    """Append (or replace) this run's summary in index_msi_run.js.

    This is the only file that changes on every run: index_msi_run.html itself
    is written once and just reads window.RUNS_DATA from here, so re-running
    the pipeline never needs to touch (or regenerate) the index page.
    """
    manifest_path = viz_dir / RUNS_DATA_FILENAME
    runs = []
    if manifest_path.exists():
        text = manifest_path.read_text()
        start, end = text.find("["), text.rfind("]")
        if start != -1 and end != -1:
            try:
                runs = json.loads(text[start:end + 1])
            except json.JSONDecodeError:
                runs = []

    runs = [r for r in runs if r.get("run_id") != run_id]
    runs.append(_run_summary(run_id, started_at, samples))
    runs.sort(key=lambda r: r["started_at"], reverse=True)
    write_data_js(manifest_path, "RUNS_DATA", runs)


def ensure_output_index_page(viz_dir: Path) -> None:
    """Write viz/index_msi_run.html. Only rewritten when the shell template itself
    changed: run data lives in index_msi_run.js, which append_run_to_manifest()
    keeps up to date."""
    index_path = viz_dir / OUTPUT_INDEX_FILENAME
    html = render_output_index_shell()
    if index_path.exists() and index_path.read_text() == html:
        return
    index_path.write_text(html)


def ensure_consolidated_sample_page(viz_dir: Path) -> None:
    """Write viz/index_msi_sample.html + its one data file index_msi_sample.js
    (window.SAMPLES_DATA): every sample of every run, one row each, gathered
    from each run's own <run_id>_index_run.js table (same columns as the run
    pages, plus a leading Run column)."""
    runs_dir = viz_dir / RUNS_DIRNAME
    if not runs_dir.is_dir():
        return

    labels = set()
    samples = []
    for run_dir in sorted(d for d in runs_dir.iterdir() if d.is_dir()):
        js_path = run_data_js_path(run_dir, run_dir.name)
        if not js_path.exists():
            continue
        try:
            run_data = read_data_js(js_path)
        except (json.JSONDecodeError, IndexError) as e:
            logger.warning(f"Error reading {js_path}: {e}")
            continue
        lib_labels = run_data.get("library_labels", {})
        for row in run_data.get("table", []):
            libs = {lib_labels.get(lib, lib): v for lib, v in row["libraries"].items()}
            labels.update(libs)
            samples.append({
                "run_id": run_data["run_id"],
                "sample_name": row["sample_name"],
                "pct_region_unstable": row["pct_region_unstable"],
                "libraries": libs,
                "consensus_total": row["consensus_total"],
                "consensus_roi": row["consensus_roi"],
                "call": row["call"],
            })

    write_data_js(viz_dir / SAMPLES_DATA_FILENAME, "SAMPLES_DATA",
                  {"library_labels": sorted(labels), "samples": samples})

    html = """<title>MSI Dashboard - Sample View</title>
__HEAD_LINKS__
<style>__CSS__
.nav-links{display:flex;gap:12px;margin-bottom:20px;}
.nav-links a{display:inline-flex;align-items:center;gap:6px;padding:8px 16px;background:var(--accent-soft);
  color:var(--accent);border-radius:6px;font-weight:500;transition:all 0.2s;text-decoration:none;}
.nav-links a:hover{background:var(--accent);color:white;}
</style>
<div class="wrap">
  <header>
    <h1>MSI Dashboard - Sample View</h1>
    <p><span id="sampleCount">0</span> samples from <span id="runCount">0</span> run(s).</p>
    <div class="nav-links">
      <a href="./__OUTPUT_INDEX_FILE__">&larr; Back to Runs</a>
    </div>
  </header>
  <div class="table-panel">
    <div class="table-scroll">
      <table class="index-table">
        <thead id="thead"></thead>
        <tbody id="tbody"></tbody>
      </table>
    </div>
  </div>
  <footer>MSIsensor-pro pipeline - Consolidated report</footer>
  <div class="scroll-top" id="scrollTop">↑</div>
</div>
<script src="./__DATA_JS_FILE__"></script>
<script>
(function(){
  "use strict";
  var D = window.SAMPLES_DATA || {library_labels: [], samples: []};
  function n(v){ return (v === null || v === undefined) ? "-" : Number(v).toLocaleString("en-US"); }
  function num(v){ return '<td class="num tabular">' + v + '</td>'; }
  var runs = new Set(D.samples.map(function(s){ return s.run_id; }));
  document.getElementById("sampleCount").textContent = D.samples.length;
  document.getElementById("runCount").textContent = runs.size;

  var groups = D.library_labels.map(function(l){
    return '<th class="num" colspan="2">UNSTABLE_SITE_' + l.toUpperCase() + '</th>';
  }).join("");
  var sub = '';
  for (var i = 0; i <= D.library_labels.length; i++) sub += '<th class="num">Total</th><th class="num">ROI</th>';
  document.getElementById("thead").innerHTML =
    '<tr><th rowspan="2">Run</th><th rowspan="2">Sample</th><th rowspan="2">% unstable</th>' + groups
    + '<th class="num" colspan="2">UNSTABLE_CONSENSUS</th><th rowspan="2">Status</th></tr><tr>' + sub + '</tr>';

  document.getElementById("tbody").innerHTML = D.samples.map(function(s){
    var cells = D.library_labels.map(function(l){
      var e = s.libraries[l];
      if (!e) return num("-") + num("-");
      if (e.failed) return '<td colspan="2"><span class="status-pill failed">failed</span></td>';
      return num(n(e.unstable_region_count)) + num(e.unstable_roi_count);
    }).join("");
    var status = (s.call && s.call !== "NA") ? '<span class="chip ' + (s.call === "MSI" ? "mid" : "low") + '">' + s.call + '</span>'
               : s.call === "NA" ? '<span class="lib-empty">NA</span>'
                        : '<span class="lib-empty">—</span>';
    return '<tr><td>' + s.run_id + '</td>'
      + '<td><a href="./__RUNS_DIR__/' + s.run_id + '/' + s.sample_name + '_index_sample.html">' + s.sample_name + '</a></td>'
      + num(s.pct_region_unstable === null ? "-" : s.pct_region_unstable + "%") + cells
      + num(n(s.consensus_total)) + num(n(s.consensus_roi)) + '<td>' + status + '</td></tr>';
  }).join("");

  var scrollTop = document.getElementById("scrollTop");
  window.addEventListener("scroll", function(){ scrollTop.classList.toggle("show", window.pageYOffset > 300); });
  scrollTop.addEventListener("click", function(){ window.scrollTo({top: 0, behavior: "smooth"}); });
})();
</script>
"""
    (viz_dir / "index_msi_sample.html").write_text(html
        .replace("__HEAD_LINKS__", HEAD_LINKS)
        .replace("__CSS__", REPORT_CSS)
        .replace("__OUTPUT_INDEX_FILE__", OUTPUT_INDEX_FILENAME)
        .replace("__RUNS_DIR__", RUNS_DIRNAME)
        .replace("__DATA_JS_FILE__", SAMPLES_DATA_FILENAME))


def main():
    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument("--config", default=None)
    pre_args, _ = pre_parser.parse_known_args()

    parser = argparse.ArgumentParser(description="Run MSIsensor-pro on a BAM via Apptainer.")
    parser.add_argument("--config", default=None,
                         help="JSON file providing parameter defaults (any matching CLI flag overrides it).")
    parser.add_argument("-b", "--bam", default=None,
                         help="Input BAM file (indexed, .bai expected), or a directory containing BAMs "
                              "(glob *.bam, all processed).")
    parser.add_argument("-o", "--output", default=None,
                         help="Output directory. Each run creates data/runs/<run_id>/<library>/ "
                              "with the raw results; HTML/JS/JSON reports are (re)generated under "
                              "viz/ (index_msi_run.html), at the output, run and sample level. "
                              "--viz-output/--data-output override where viz/ and data/ actually "
                              "live (e.g. on separate disks); this path is only their default parent.")
    parser.add_argument("--viz-output", default=None,
                         help="Directory for every generated .html/.js/.json report (defaults to "
                              "<output>/viz). Lets viz/ live on a different disk/path than data/.")
    parser.add_argument("--data-output", default=None,
                         help="Directory for every raw msisensor-pro output file (defaults to "
                              "<output>/data). Lets data/ live on a different disk/path than viz/.")
    parser.add_argument("--log-output", default=None,
                         help="Directory for per-run log files (defaults to <output>/log). "
                              "Lets log/ live on a different disk/path than viz/ and data/.")
    parser.add_argument("-g", "--genome", default=None,
                         help="Path to the reference fasta used to generate the MSIsensor-pro -d file "
                              "(via 'msisensor-pro scan'). Required (via CLI flag or --config), "
                              "unless --pro-list is given.")
    parser.add_argument("-s", "--sif", default=str(DEFAULT_SIF),
                         help="Path to the MSIsensor-pro Apptainer image (.sif).")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads.")
    parser.add_argument("--coverage-threshold", type=int, default=1,
                         help="Passed to msisensor-pro pro's -c: coverage threshold for MSI analysis "
                              "(msisensor-pro's own default is 15 for WGS, 20 for WXS). This pipeline "
                              "defaults to 1.")
    parser.add_argument("--pro-list", default=None,
                         help="Already generated MSIsensor-pro -d (microsatellites) file (skips the scan).")
    parser.add_argument("--panel-bed", default=None,
                         help="Panel/target BED file. If given, intersects it with the genome-wide "
                              "microsatellite BED (via bedtools) into <...>.panel.bed, "
                              "next to the reference list/BED.")
    parser.add_argument("--roi-bed", default=None,
                         help="Region-of-interest BED (e.g. the classic BAT26/NR21/... MSI marker "
                              "loci). If given, intersects it with the genome-wide microsatellite "
                              "BED (via bedtools) into <...>.roi.bed, same mechanism as --panel-bed.")
    parser.add_argument("--bedtools-sif", default=str(DEFAULT_BEDTOOLS_SIF),
                         help="Path to the bedtools Apptainer image (.sif), used with --panel-bed.")
    parser.add_argument("--recursive", action="store_true",
                         help="If -b is a directory, look for *.bam recursively in subdirectories.")

    CONFIG_KEY_ALIASES = {
        "pathToRun": "bam",
        "msisensor-pro": "sif",
        "bedtools": "bedtools_sif",
    }

    if pre_args.config:
        config = load_config(Path(pre_args.config))
        for alias, dest in CONFIG_KEY_ALIASES.items():
            if alias in config:
                config[dest] = config.pop(alias)
        unknown = set(config) - {a.dest for a in parser._actions}
        if unknown:
            sys.exit(f"Error: unknown key(s) in config file: {', '.join(sorted(unknown))}")
        parser.set_defaults(**config)

    args = parser.parse_args()

    if not args.bam or not args.output:
        sys.exit("Error: --bam/-b and --output/-o are required (via CLI flags or --config).")
    if not args.genome and not args.pro_list:
        sys.exit("Error: --genome/-g is required (via CLI flag or --config), unless --pro-list is given.")

    bam_arg = Path(args.bam)

    # Name the run after the input directory (e.g. "Run1"), so re-running it
    # updates that same run in place instead of piling up timestamped folders.
    # Falls back to the BAM's filename when -b/pathToRun points at a single file.
    run_id = bam_arg.name if bam_arg.is_dir() else bam_arg.stem

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    # viz/, data/ and log/ default to subfolders of --output, but
    # --viz-output/--data-output/--log-output (CLI flag or JSON key) let each
    # live at its own path instead -- e.g. data/ on a large/slow disk, viz/
    # and log/ on a small/fast one.
    viz_dir = Path(args.viz_output) if args.viz_output else output_dir / VIZ_DIRNAME
    data_dir = Path(args.data_output) if args.data_output else output_dir / DATA_DIRNAME
    log_dir = Path(args.log_output) if args.log_output else output_dir / LOG_DIRNAME
    viz_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    setup_logging(run_id, log_dir)

    bam_files = find_bam_files(bam_arg, args.recursive)

    # Every generated .html/.js/.json report lives under viz/; every raw
    # msisensor-pro output file lives under data/ -- both mirror the same
    # runs/<run_id>/<library>/ layout, just rooted differently.
    viz_run_dir = viz_dir / RUNS_DIRNAME / run_id
    data_run_dir = data_dir / RUNS_DIRNAME / run_id
    viz_run_dir.mkdir(parents=True, exist_ok=True)
    data_run_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"{len(bam_files)} BAM file(s) to process -> run {run_id}")

    sif_path, list_path, roi_correspondence_path = prepare_reference(args)
    genome = args.genome

    failures = []
    sample_metas = []
    for bam_path in bam_files:
        sample_name = bam_path.stem
        # A directory input can hold several libraries as subfolders (e.g.
        # Run1/alignements_libA/, Run1/alignements_libB/), each with samples
        # of the same name. Mirror that subfolder under the run so results
        # never collide, and keep the subfolder name as the sample's library.
        # No per-sample subfolder inside it, though: every sample's files sit
        # directly in the library folder, named with the sample as a prefix
        # (exactly how msisensor-pro's own -o prefix already names them).
        rel_dir = bam_path.parent.relative_to(bam_arg) if bam_arg.is_dir() else Path(".")
        library = rel_dir.as_posix() if rel_dir != Path(".") else None
        sample_key = (rel_dir / sample_name).as_posix()
        data_library_dir = data_run_dir / rel_dir
        data_library_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"=== Processing {bam_path.name} (sample: {sample_key}) ===")
        try:
            process_bam(bam_path, data_library_dir, sample_name, sif_path, list_path, args)
            meta = compute_sample_meta(data_library_dir, sample_name, bam_path, genome,
                                       run_id, library, sample_key)
            sample_metas.append(meta)
        except SystemExit:
            logger.error(f"Failed on {bam_path.name} (see above for details)")
            failures.append(bam_path.name)
            sample_metas.append({"sample_name": sample_name, "library": library, "sample_key": sample_key, "failed": True})

    rows_by_sample = load_sample_library_rows(data_run_dir, sample_metas, roi_correspondence_path)
    started_at = write_run_index(viz_run_dir, run_id, sample_metas, rows_by_sample)
    write_sample_group_indexes(viz_run_dir, run_id, sample_metas, rows_by_sample)

    # index_msi_run.html is written once and never rewritten (see
    # ensure_output_index_page), so if it's already there but index_msi_run.js got
    # deleted by accident, the normal incremental append below would only
    # restore *this* run -- every other past run would vanish from the list.
    # Rebuild the full history from viz/runs/*/*_index_run.js first in that case.
    if (viz_dir / OUTPUT_INDEX_FILENAME).exists() and not (viz_dir / RUNS_DATA_FILENAME).exists():
        logger.info(f"{RUNS_DATA_FILENAME} missing but {OUTPUT_INDEX_FILENAME} exists -- rebuilding it from viz/runs/*/*_index_run.js")
        rebuild_runs_manifest(viz_dir)

    append_run_to_manifest(viz_dir, run_id, started_at, sample_metas)
    ensure_output_index_page(viz_dir)
    ensure_consolidated_sample_page(viz_dir)
    logger.info(f"Report: {viz_dir / OUTPUT_INDEX_FILENAME}")

    if failures:
        fatal(f"{len(failures)}/{len(bam_files)} BAM(s) failed: {', '.join(failures)}")


if __name__ == "__main__":
    main()
