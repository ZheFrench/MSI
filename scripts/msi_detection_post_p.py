#!/usr/bin/env python3
"""
msi_detection_post_p.py

Run after the msi_detection-0.1.0.wdl workflow, outside it (it writes into
the reports folder). The steps inside the workflow are in
msi_detection_main_p.py (reference lists) and msi_detection_report_p.py (this
run's report), whose page code it reuses (import: both scripts sit in
scripts/).

    index        --viz DIR [--run-outputs outputs.json]
        Copies the run's runs/<run_id>/ into DIR/runs/, then rebuilds the
        dashboard from every DIR/runs/*/:

    DIR/index_msi_run.js, index_msi_sample.js     - rewritten (all runs, all samples)
    DIR/index_msi_run.html, index_msi_sample.html - written once (only if the template changed)

Usage (BASH.sprocket.cron.sh / BASH.sprocket.batch.sh, after each run):
    msi_detection_post_p.py index --viz VIZ --run-outputs <run>/outputs.json
"""

import argparse
import json
import logging
import shutil
import sys
from pathlib import Path

# msi_detection_report_p.py sits next to this script (scripts/): importable
# whatever the current directory or python's -I/-P flags.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from msi_detection_report_p import (HEAD_LINKS, OUTPUT_INDEX_FILENAME, REPORT_CSS,  # noqa: E402
                                  RUNS_DIRNAME, logger, run_data_js_path, write_data_js)

# Dashboard data files, in the reports folder
RUNS_DATA_FILENAME = "index_msi_run.js"
SAMPLES_DATA_FILENAME = "index_msi_sample.js"


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
    viz/runs/<run_id>/<run_id>_index_run.js: viz/runs/ is the only source of
    truth, so the list always matches the run folders present.
    index_msi_run.html itself never needs touching, it only ever reads this
    file via window.RUNS_DATA."""
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


def ensure_output_index_page(viz_dir: Path) -> None:
    """Write viz/index_msi_run.html. Only rewritten when the shell template itself
    changed: run data lives in index_msi_run.js (rebuild_runs_manifest)."""
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
    rules = []  # each run's own sample status rule (its p_process), shown once each
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
        rule = run_data.get("status_rule")
        if rule and rule not in rules:
            rules.append(rule)
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
    <p>__SAMPLE_STATUS_RULE__</p>
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
    <div class="table-foot">
      <span id="rangeLabel">-</span>
      <div class="pager">
        <button id="prevBtn" type="button">&larr; Prev</button>
        <span id="pageLabel">Page 1</span>
        <button id="nextBtn" type="button">Next &rarr;</button>
      </div>
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
  // UNSTABLE_CONSENSUS (headers in black, values in bold), as on the run pages.
  var sub = '';
  for (var i = 0; i < D.library_labels.length; i++) sub += '<th class="num">Total</th><th class="num">ROI</th>';
  sub += '<th class="num consensus">Total</th><th class="num consensus">ROI</th>';
  document.getElementById("thead").innerHTML =
    '<tr><th rowspan="2">Run</th><th rowspan="2">Sample</th><th rowspan="2">% unstable</th>' + groups
    + '<th class="num consensus" colspan="2">UNSTABLE_CONSENSUS</th><th rowspan="2">Status</th></tr><tr>' + sub + '</tr>';

  // PAGE_SIZE samples per page, Prev/Next to browse the others.
  var PAGE_SIZE = 20, page = 0, last = Math.max(0, Math.ceil(D.samples.length / PAGE_SIZE) - 1);
  function renderPage(){
    var start = page * PAGE_SIZE;
    document.getElementById("tbody").innerHTML = D.samples.slice(start, start + PAGE_SIZE).map(function(s){
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
        + '<td class="num tabular consensus">' + n(s.consensus_total) + '</td><td class="num tabular consensus">' + n(s.consensus_roi) + '</td>'
        + '<td>' + status + '</td></tr>';
    }).join("");
    var total = D.samples.length;
    document.getElementById("rangeLabel").textContent = total
      ? "Samples " + (start + 1) + "\u2013" + Math.min(start + PAGE_SIZE, total) + " of " + total.toLocaleString("en-US") : "No sample";
    document.getElementById("pageLabel").textContent = "Page " + (page + 1) + " / " + (last + 1);
    document.getElementById("prevBtn").disabled = page <= 0;
    document.getElementById("nextBtn").disabled = page >= last;
  }
  document.getElementById("prevBtn").addEventListener("click", function(){ if (page > 0){ page--; renderPage(); } });
  document.getElementById("nextBtn").addEventListener("click", function(){ if (page < last){ page++; renderPage(); } });
  renderPage();

  var scrollTop = document.getElementById("scrollTop");
  window.addEventListener("scroll", function(){ scrollTop.classList.toggle("show", window.pageYOffset > 300); });
  scrollTop.addEventListener("click", function(){ window.scrollTo({top: 0, behavior: "smooth"}); });
})();
</script>
"""
    html = (html
        .replace("__HEAD_LINKS__", HEAD_LINKS)
        .replace("__CSS__", REPORT_CSS)
        .replace("__OUTPUT_INDEX_FILE__", OUTPUT_INDEX_FILENAME)
        .replace("__RUNS_DIR__", RUNS_DIRNAME)
        .replace("__SAMPLE_STATUS_RULE__", "<br>".join(rules))
        .replace("__DATA_JS_FILE__", SAMPLES_DATA_FILENAME))
    # Written once: only rewritten when the template itself changed.
    html_path = viz_dir / "index_msi_sample.html"
    if not (html_path.exists() and html_path.read_text() == html):
        html_path.write_text(html)


# ---------------------------------------------------------------------------
# index: add a run's pages to the reports folder, rebuild the dashboard
# ---------------------------------------------------------------------------

def index(args) -> None:
    """Add this run's viz/runs/<run_id>/ (copied from the ta03_report call
    folder named by --run-outputs, so the sprocket run stays complete) to the
    reports folder, then rebuild the dashboard from every viz/runs/*/: index_msi_run.js and
    index_msi_sample.js are rewritten, index_msi_run.html and
    index_msi_sample.html only when their template changed. Re-runnable at
    any time (e.g. after deleting a run folder)."""
    viz_dir = Path(args.viz)
    (viz_dir / RUNS_DIRNAME).mkdir(parents=True, exist_ok=True)

    if args.run_outputs:
        outputs = json.loads(Path(args.run_outputs).read_text())
        new_runs = Path(outputs["msi_detection.report_viz"]) / RUNS_DIRNAME
        for run_dir in sorted(d for d in new_runs.iterdir() if d.is_dir()):
            dest = viz_dir / RUNS_DIRNAME / run_dir.name
            if dest.exists():
                shutil.rmtree(dest)
                logger.info(f"Replacing {dest}")
            shutil.copytree(run_dir, dest)
            logger.info(f"Added {dest}")

    runs = rebuild_runs_manifest(viz_dir)
    ensure_output_index_page(viz_dir)
    ensure_consolidated_sample_page(viz_dir)
    logger.info(f"Dashboard: {viz_dir / OUTPUT_INDEX_FILENAME} ({len(runs)} run(s))")


def main():
    parser = argparse.ArgumentParser(description="Dashboard of the msi_detection-0.1.0.wdl runs (after the workflow).")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("index", help="Add a run's pages to the reports folder and rebuild the dashboard.")
    p.add_argument("--viz", required=True, help="Reports folder (VIZ of BASH.sprocket.cron.sh).")
    p.add_argument("--run-outputs", default=None,
                   help="outputs.json of a run: its report_viz runs/<run_id>/ is copied into --viz/runs/.")

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    if args.command == "index":
        index(args)


if __name__ == "__main__":
    main()
