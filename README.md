# msi_detection-0.1.0.wdl

Everything about the WDL workflow lives in this folder (`pipeline-jp/`).


WDL workflow run with [sprocket](https://sprocket.bio). It runs MSIsensor-pro (tumor-only) on every BAM of a run directory,
restricts the analysis to a sequencing panel and/or the classic MSI markers (BAT26, NR21, ...), and
updates the HTML reports (all runs → one run → one sample).

[![msi_detection metro map](docs/diagram/msi_detection_metromap.preview.png)](docs/diagram/msi_detection_metromap.preview.png)

(Full-size image: [`docs/diagram/msi_detection_metromap.preview.png`](docs/diagram/msi_detection_metromap.preview.png), copied there by `BASH.doc.sh`; editable diagram (local, not in git): `diagram/msi_detection_metromap.drawio`, PNG re-exported with `diagram/BASH.render_png.sh`.)

## Contents

- [Configuration (to set at install)](#configuration-to-set-at-install)
- [Quick start](#quick-start)
  - [Batch mode: params JSON and TSV export](#batch-mode-params-json-and-tsv-export)
- [Tasks](#tasks)
  - [Resources](#resources)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Containers: development vs production](#containers-development-vs-production)
- [Files](#files)
- [HTML documentation](#html-documentation)

## Configuration (to set at install)

Nothing is read from the environment: every path and Slurm setting is hardcoded in the files below
(values here: this workstation). Every path must be the same on all the Slurm nodes (shared file
system): the containers mount them at the same path.

**`BASH.sprocket.cron.sh`** (cron mode)

| Variable / line | Value | Role |
|---|---|---|
| `PATH` (line 6, commented) | `/home/jp/.cargo/bin:$PATH` | cron's PATH is minimal: uncomment it if `sprocket` is not found |
| `WDL_DIR` | `/media/jp/SSD2/CHU/code/pipeline-jp` | this folder |
| `VIZ` | `/media/jp/SSD2/CHU/test/viz` | reports and dashboard (written by `index`) |
| `WORK` | `/media/jp/SSD2/CHU/test/work` | sprocket's runs (`sprocket/`) and the params JSON (`config/`) |
| `RUNS_DIR` | `/media/jp/SSD2/CHU/input/fake-data/runs` | one folder per sequencing run, with its `SampleSheet.csv` |
| `BAMS_DIR` | `/media/jp/SSD2/CHU/input/fake-data/fake_bam` | one BAM folder per Experiment Name |
| `LOCKER_ADIVAR` | `/media/jp/SSD2/CHU/config/locker/adivar-testNS.json` | ADIVAR's locker (read only) |
| `LOCKER_MSI` | `/media/jp/SSD2/CHU/config/locker/msi-sensor.locker.txt` | MSI locker (runs already taken), created if missing |
| `srun` in `python_slurm_job` | `--partition=prod --cpus-per-task=1 --mem=1G` | Slurm job of the python steps (`pick_run`, `experiment`, `pre_process`, `index`) |

**`companion-scripts/BASH.sprocket.batch.sh`** (batch mode): `WDL_DIR` and the same `srun` line (its output folder
`OUT_DIR`, which gets `viz/` and `work/`, is a required argument).

**`backends.conf/sprocket.prod.toml`**

| Key | Value | Role |
|---|---|---|
| `[run.backends.slurm.sbatch] args` | `["--partition=prod"]` | partition (and any other `sbatch` option) of the tasks |
| `[run.task] cache_dir` | `/media/jp/SSD2/CHU/test/work/sprocket_cache` | sprocket's call cache (shared by all the runs) |
| `[run.backends.slurm] interval` | `5` | seconds between two `sacct` polls: short for local tests, raise it on the cluster |
| `max_concurrency`, `[run.workflow.scatter] concurrency` | `50`, `10` | jobs submitted at once, scatter elements in parallel |

**`params/template.json`** (copied by `pre_process` for each run; `run` and `samples` stay empty, filled per run)

| Key | Value |
|---|---|
| `msi_detection.references.genome` | `/media/jp/SSD2/CHU/refData/hg19/Homo_sapiens.GRCh37.75.wout_contigs.fa` |
| `msi_detection.references.bed_panel` | `/media/jp/SSD2/CHU/refData/MSI/static/Panel_Juno_v3_Juno.amplicons.bed` |
| `msi_detection.references.bed_roi` | `/media/jp/SSD2/CHU/refData/MSI/static/ROI.hg19.bed` |
| `msi_detection.params` | `threads` (8, CPUs of each `msisensor-pro pro` only), `coverage_threshold`, `p_process` thresholds |

**`params/containers.json`**: the workflow's images and python helpers, the only place they are
set (absolute paths; images as `file://` + path). These five WDL inputs have no default: both
scripts pass this file before each run's params JSON (`sprocket run ... @containers.json @<params>.json`;
for a key in both, the last file wins, so a params JSON can still override an image for its run).
Forgetting it stops sprocket before anything runs (`missing required input`).

| Key | Value |
|---|---|
| `msi_detection.msisensorpro_container` | `file:///media/jp/SSD2/CHU/code/pipeline-jp/apptainer-recipes/msisensor-pro-1.3.0.sif` |
| `msi_detection.bedtools_container` | `file:///media/jp/SSD2/CHU/code/pipeline-jp/apptainer-recipes/bedtools-2.31.1.sif` |
| `msi_detection.helper_container` | `file:///media/jp/SSD2/CHU/code/pipeline-jp/apptainer-recipes/python-3.12.sif` |
| `msi_detection.msi_main_helper` | `/media/jp/SSD2/CHU/code/pipeline-jp/main-scripts/msi_detection_main_p.py` |
| `msi_detection.msi_report_helper` | `/media/jp/SSD2/CHU/code/pipeline-jp/main-scripts/msi_detection_report_p.py` |

The images of the scripts' own Slurm jobs follow `WDL_DIR`: `$WDL_DIR/apptainer-recipes/python-3.12.sif`
(`python_slurm_job` of `BASH.sprocket.cron.sh` / `BASH.sprocket.batch.sh`: pick_run, experiment,
pre_process, index) and `$WDL_DIR/apptainer-recipes/sprocket-doc-0.31.0.sif` (`BASH.doc.sh`).

**crontab** (not in this folder), for example every 30 min with a log:

```sh
*/30 * * * * /media/jp/SSD2/CHU/code/pipeline-jp/BASH.sprocket.cron.sh >> /media/jp/SSD2/CHU/test/cron.log 2>&1
```

## Quick start

```sh
cd /media/jp/SSD2/CHU/code/pipeline-jp
./BASH.sprocket.cron.sh                       # cron mode: pre_process picks the run left to analyse, then runs it
./companion-scripts/BASH.sprocket.batch.sh CONFIG_DIR OUT_DIR   # batch mode: every params JSON of CONFIG_DIR, outputs in OUT_DIR
```

The batch mode runs `sprocket` in the current shell (hours for a whole group): launch it in
`screen` / `tmux`, or a lost terminal / SSH session stops it. Not needed for the cron mode (no terminal).

`BASH.sprocket.cron.sh` (cron) first runs `msi_detection_pre_p.py` three times, each container
mounting only what that step reads:

1. `--print-run` (job `msi-pick_run`; the two lockers): prints the run finished by ADIVAR and not yet
   in the MSI locker;
2. `--print-experiment` (job `msi-experiment`; + that run's `$RUNS_DIR/<run>/SampleSheet.csv`): checks
   its `[Header]` and prints its Experiment Name, the name of its BAM folder;
3. without option (job `msi-pre_process`; + that run's `$BAMS_DIR/<Experiment Name>/`, the MSI
   locker and `$WORK/config/` writable): writes its params JSON in `$WORK/config/` and adds the run
   to the locker.

Neither the whole `$RUNS_DIR` nor the whole `$BAMS_DIR` is mounted. A rejected sample sheet or a
missing BAM folder stops the script before the run is locked (retried at the next launch).
Then both scripts, for each params JSON:

1. sets the destinations, hardcoded in the script (not in the params): `VIZ` (reports folder) and
   `WORK` (sprocket's runs), created if needed;
2. runs `sprocket run -s -c $SPROCKET_CONFIG` (`backends.conf/sprocket.prod.toml`) in
   `$WORK/sprocket/runs/msi_detection/<timestamp>/`, one folder per task;
3. runs `msi_detection_post_p.py index` as a Slurm job (`srun`, in `python-3.12.sif`): copies the run's
   `viz/runs/<run_id>/` into `$VIZ/runs/` and rebuilds the dashboard from every
   `$VIZ/runs/*/` (`index_msi_run.js`, `index_msi_sample.js`; the two `.html` are written once,
   only rewritten if their template changes).

The workflow itself never writes outside its tasks: the dashboard covers every run, a state that
changes with each run, so it is built after the run (`index`), not by a task (a task reading the
current reports would not be reproducible nor cacheable). `index` can be re-run at any time, e.g.
after deleting a run folder. Everything else stays in sprocket's run folder (paths in its
`outputs.json`), every task keeping its `work/` folder (~300 MB per run).

### Batch mode: params JSON and TSV export

The params JSON of a batch come from `companion-scripts/make_configs.py`: given a group folder of
BAM folders (one per run: `input/colon-msi`, `input/foie-mss`...) and run names, on the command line
and/or in a text file, it writes `<run>_params.json` (copy of `params/template.json`,
`run` = `<group folder>/<run>`) in `--output-dir` (required: the launcher's `OUTPUT_DIR`). The
group folder picks the run, the same run name being in several groups.

The runs file has one line per run, or per run + sample ID to analyse only some samples
(columns separated by tabs or spaces, `#` comments allowed):

```
Run154_Juno                  # every BAM of the run (but the controls)
Run171_Juno  C22D03556T      # only these samples of Run171_Juno
Run171_Juno  C22D03587T
```

The sample IDs of a run go to `"msi_detection.inputs".samples`: `ta01_find_bams` then only keeps the
BAMs whose name contains one of them (every library). A run without BAM, or with a sample ID
matching no BAM of the run (checked before the workflow), gets no params JSON: its lines are written
in `<output-dir>/failed-runs.txt` (same format as the runs file, for a rerun; removed when no run
fails) and the script exits 1, the other runs being written.

After the batch, `companion-scripts/export_tsv.py` writes the results of its reports folder as a
TSV, one line per sample: `results_msi_<date>.tsv` (date `YYYY-MM-DD_HH-MM-SS`: each export is a new file), in the launcher's `OUTPUT_DIR` (required), in the reports folder with the python script alone (`--output-dir`). It recomputes nothing: it reads
the dashboard data (`index_msi_sample.js`, then each sample's `runs/<run_id>/<sample>_index_sample.js`),
so the table matches the dashboard. Columns: `run_id`, `sample_name`, `call`, `consensus_roi`,
`consensus_total`, `pct_region_unstable`; per library `<lib>_unstable_region_count`,
`<lib>_unstable_roi_count`; per ROI marker `<marker>` (its site status, the same for every library),
then `<lib>_<marker>_pro_p` and `<lib>_<marker>_cov` for each library (no mean). Empty cell: no value.

Both run through a launcher, like the other python steps: a Slurm job (`srun`, jobs
`msi-make_configs` / `msi-export_tsv`) in `python-3.12.sif` (`python3 -I`), mounting the code
read-only, what the script reads (read-only) and the folder it writes:

```sh
S=/media/jp/SSD2/CHU
./companion-scripts/BASH.make_configs.sh $S/input/foie-mss $S/config/runs-test.txt $S/config/config-test
./companion-scripts/BASH.sprocket.batch.sh $S/config/config-test $S/output/config-test
./companion-scripts/BASH.export_tsv.sh $S/output/config-test/viz $S/output/config-test      # -> $S/output/config-test/results_msi_<date>.tsv
```

The python scripts can also be run directly on the host (only the standard library):
`make_configs.py BAMS_DIR [RUN ...] [--runs-file TXT] --output-dir DIR`,
`export_tsv.py --viz VIZ_DIR [--output-dir DIR | --output TSV]`.

Run by hand, with a params JSON written by `pre_process` (or a filled copy of `params/template.json`):

```sh
sprocket run -s -c backends.conf/sprocket.prod.toml msi_detection-0.1.0.wdl @params/containers.json @$WORK/config/<date>-<run>-<Experiment Name>.json
```

## Tasks

In execution order. Tools are wrapped in **modules**, imported by the workflow:

```wdl
import "module/msisensorpro-1.3.0.wdl" as msisensorpro   # tasks scan, pro
import "module/bedtools-2.31.1.wdl" as bedtools           # task intersect
```

Module tasks are generic (the tool's own options only), declare their image in
`requirements.container` and can be reused by other workflows; the pipeline-specific steps are
tasks of `msi_detection-0.1.0.wdl`, in a generic python3 + awk image (`helper_container`,
`apptainer-recipes/python-3.12.sif`); the python script is given to them as a `File`, not built into the image.

Calls are named after their line on the diagram, `t<line><rank>_<step>`: `ta` BAMs (green),
`tb` msisensor-pro list (blue, panel steps included), `tc` ROI markers (red). A WDL name can't
start with a digit, hence the letter; sprocket's `calls/<call>/` folders group the same way.

| Call (task) | Runs | What it does | Main output |
|---|---|---|---|
| `ta01_find_bams` (`find_bams`) | python image, bash | Lists every `*.bam` of the run directory (`inputs.run`) (sub-folders with `params.recursive`), with its index (`.bam.bai` or `.bai`), its library (folder relative to `input_dir`, `.` if none) and its sample name (BAM basename). `inputs.run` is a `Directory` (mounted by the engine); paths are returned relative to it. With `inputs.samples`, keeps only the BAMs whose name contains one of these IDs (every library). Fails if a BAM has no index, no BAM is found, or a sample ID matches no BAM. | `bams`: one row per BAM |
| `tb01_sensorpro_scan` (`msisensorpro.scan`) | msisensor-pro | `msisensor-pro scan` of the genome: every microsatellite site. Skipped when `references.pro_list` is given. Not redone from one run to the next: sprocket's call cache (`[run.task] cache = "on"`) reuses the previous scan (same genome and image), as long as that run's folder is kept. | `<genome>.msisensor_pro.list` |
| `tb02_list_to_bed` (`list_to_bed`) | python image | Converts the list to BED6 (start = location − 1, end = start + unit length × repeat times). | `<list>.bed` |
| `tb03_panel_chrom` (`align_chrom_naming`) | python image, awk | Only with `bed_panel`. Copy of the panel BED with its chromosome names (`chr1` / `1`) aligned on the list's: bedtools matches them literally. | `<panel>.chrnorm.bed` |
| `tb04_panel_intersect` (`bedtools.intersect`) | bedtools | `bedtools intersect -wa -wb -f 1.0`: the sites fully inside a panel region, one row per overlap. | `<list>.panel.intersect.tsv` |
| `tb05_panel_correspondence` (`name_correspondence`) | python image, awk | Turns the intersection into a BED6 named `<site>#<probe>#<chrom:start-end>`. | `<list>.panel.bed` |
| `tb06_mark_list` (`mark_list`) | python image | Only with `bed_panel`. Full list, same rows, `filter` = `PASS` for the panel sites, `NO` elsewhere: **the list analysed by msisensor-pro**. | `<list>.panel.list` |
| `tc01_roi_chrom`, `tc02_roi_intersect`, `tc03_roi_correspondence` | python image / bedtools / python image | Only with `bed_roi`. Same three steps with the ROI markers, on the panel sites when there is a panel (genome-wide otherwise). Name ends with `#<marker>`. | `<list>.panel.roi.bed` |
| `ta02_sensorpro_pro` × N (`msisensorpro.pro`) | msisensor-pro | Scatter, one per BAM: `msisensor-pro pro -d <analysis list> -t <bam> -b threads -c coverage_threshold`. A failure on any BAM fails the workflow. | `<sample>`, `_all`, `_dis`, `_unstable` |
| `ta03_report` (`report`) | python image | Gets every BAM's results as `File`s, copies them to `data/<run_id>/<library>/`, builds this run's reports in `viz/runs/<run_id>/` (run page, sample pages; the all-runs dashboard is built by `index` after the run) and writes the log; writes nothing outside its task. | `viz/`, `data/`, `log/` (`Directory`) |

The python steps of the workflow are subcommands of two self-contained scripts (a WDL task gets
its script alone): [`main-scripts/msi_detection_main_p.py`](main-scripts/msi_detection_main_p.py)
(`msi_main_helper`: `list-to-bed`, `mark-list`, the reference list) and
[`main-scripts/msi_detection_report_p.py`](main-scripts/msi_detection_report_p.py) (`msi_report_helper`:
`report`). They are kept apart because sprocket's call cache keys on each task's input files:
editing the report must not invalidate the list steps, whose outputs are whole-genome lists
(~100 MB each). [`main-scripts/msi_detection_post_p.py`](main-scripts/msi_detection_post_p.py) only holds
`index`, run after the workflow by the BASH scripts; it imports the page code (CSS, `<head>`,
`.js` writer) of `msi_detection_report_p.py`, next to it in `main-scripts/`.

### Resources

What each job reserves from Slurm (`requirements` of the task: `cpu`, `memory`; sprocket passes them
to `sbatch`). Only `ta02_sensorpro_pro` takes `params.threads`: every other task keeps its default
`threads = 1`. No time limit is set: the partition's applies.

| Job | Task / script | CPU | Memory | Set in |
|---|---|---|---|---|
| `msi-pick_run`, `msi-experiment`, `msi-pre_process` (cron only) | `msi_detection_pre_p.py` | 1 | 1 GiB | `srun` line of `BASH.sprocket.cron.sh` |
| `ta01_find_bams` | `find_bams` | 1 | 512 MiB | `msi_detection-0.1.0.wdl` |
| `tb01_sensorpro_scan` | `msisensorpro.scan` | 1 | 4 GiB | `module/msisensorpro-1.3.0.wdl` |
| `tb02_list_to_bed` | `list_to_bed` | 1 | 2 GiB | `msi_detection-0.1.0.wdl` |
| `tb03_panel_chrom`, `tc01_roi_chrom` | `align_chrom_naming` | 1 | 512 MiB | `msi_detection-0.1.0.wdl` |
| `tb04_panel_intersect`, `tc02_roi_intersect` | `bedtools.intersect` | 1 | 2 GiB | `module/bedtools-2.31.1.wdl` |
| `tb05_panel_correspondence`, `tc03_roi_correspondence` | `name_correspondence` | 1 | 512 MiB | `msi_detection-0.1.0.wdl` |
| `tb06_mark_list` | `mark_list` | 1 | 2 GiB | `msi_detection-0.1.0.wdl` |
| `ta02_sensorpro_pro` × N (one per BAM) | `msisensorpro.pro` | `params.threads` (8: `pro -b`) | 4 GiB | CPU: params JSON; memory: `module/msisensorpro-1.3.0.wdl` |
| `ta03_report` | `report` | 1 | 2 GiB | `msi_detection-0.1.0.wdl` |
| `msi-index` | `msi_detection_post_p.py index` | 1 | 1 GiB | `srun` line of `BASH.sprocket.cron.sh` / `BASH.sprocket.batch.sh` |
| `msi-make_configs`, `msi-export_tsv` (batch only, run by hand) | `make_configs.py`, `export_tsv.py` | 1 | 1 GiB | `srun` line of `BASH.make_configs.sh` / `BASH.export_tsv.sh` |

At most `max_concurrency` (50) jobs at once and `[run.workflow.scatter] concurrency` (10) BAMs in
parallel (`sprocket.prod.toml`): up to 10 × 8 = 80 CPUs and 40 GiB during `ta02_sensorpro_pro`.

## Inputs

Keys of the params JSON, all prefixed `msi_detection.` (see [`params/template.json`](params/template.json)).
`inputs`, `references` and `params` are objects (WDL structs); where results go is not a
params key (see `BASH.sprocket.cron.sh`):

```json
"msi_detection.inputs":     {"run": "...", "samples": ["E25D04852T", "..."]},
"msi_detection.references": {"genome": "...", "bed_panel": "...", "bed_roi": "...", "pro_list": "..."},
"msi_detection.params":     {"threads": 8, "coverage_threshold": 1, "recursive": true, "p_process": {"min_reads": 200, "pval_cutoff": 0.1, "sample_na_fraction": 0.5, "consensus_roi_cutoff": 3}},
```

Optional members can be left out of their object.

| Object / input | Member | Required | Description |
|---|---|---|---|
| `inputs` (`InputDirs`) | `run` | yes | Run directory; its name is the `run_id` |
| | `samples` | no | Sample IDs: only the BAMs whose name contains one of them are analysed (an ID with no BAM fails the run). Without it: every BAM but the controls |
| `references` (`References`) | `genome` | yes | Fasta the BAMs are aligned on |
| | `pro_list` | no | Existing scan list: skips the scan |
| | `bed_panel` | no | Plain panel BED (no `track` line, e.g. `*.amplicons.bed`) |
| | `bed_roi` | no | MSI markers BED (e.g. `ROI.hg19.bed`) |
| `params` (`RunParams`) | `threads` | yes | `pro -b` and the CPUs reserved by each `ta02_sensorpro_pro` (only task using it) |
| | `coverage_threshold` | yes | `pro -c` |
| | `recursive` | yes | Also look for BAMs in sub-folders (one per library) |
| | `p_process.min_reads` | yes | Reports: a site counts (% unstable, site status) only with ≥ this many reads (e.g. 200) |
| | `p_process.pval_cutoff` | yes | Reports: a site is unstable when its p-value is above this (e.g. 0.1) |
| | `p_process.sample_na_fraction` | yes | Reports: sample status NA when more than this fraction of the loci listed on its sample page are NA or MSU (e.g. 0.5): ROI markers, and sites with a p-value in at least one library (loci without any, e.g. 0 reads everywhere, are left out) |
| | `p_process.consensus_roi_cutoff` | yes | Reports: sample status MSI when at least this many ROI markers are unstable in consensus, MSS otherwise (e.g. 3) |
| `msisensorpro_container`, `bedtools_container` | | yes (`params/containers.json`) | Image of each tool (see below) |
| `helper_container` | | yes (`params/containers.json`) | Generic python3 + awk image of the python/bash tasks |
| `msi_main_helper` | | yes (`params/containers.json`) | Python reference list steps, a `File` |
| `msi_report_helper` | | yes (`params/containers.json`) | Python report step, a `File` |

The genome, panel and ROI must use the same build as the BAMs (e.g. GRCh37 + `ROI.hg19.bed`).

## Outputs

Returned by `ta03_report` as `Directory` outputs, in its `work/` folder of the run directory;
`index` then adds `runs/<run_id>/` to `$VIZ` and builds the dashboard there:

```
<VIZ>/index_msi_run.html                       all runs (written once, reads index_msi_run.js)
<VIZ>/index_msi_sample.html                    all samples of all runs (written once, reads index_msi_sample.js)
<VIZ>/runs/<run_id>/<run_id>_index_run.html     one row per sample, one column group per library
<VIZ>/runs/<run_id>/<sample>_index_sample.html  one sample, pooled across its libraries
data/<run_id>/<library>/<sample>{,_all,_dis,_unstable}
log/<run_id>_MSI_<timestamp>.log
```

sprocket's own run directory (with `BASH.sprocket.cron.sh`: `$WORK/sprocket/runs/msi_detection/<timestamp>/`, `_latest` → last run):

```
output.log        sprocket's log
inputs.json       resolved inputs
outputs.json      workflow outputs (also printed by sprocket run)
calls/<call>/attempts/0/{command,stdout,stderr,work/}   each call (ta01_find_bams, ...)
```

Workflow outputs (`outputs.json`): `msi_reference_list`, `msi_analysis_list`,
`panel_correspondence_bed`, `roi_correspondence_bed`, `msisensorpro_all`,
`msisensorpro_unstable`, `report_viz`, `report_data`, `report_log`
(Directory).

## Containers: development vs production

Tool tasks only name their image and call the tool directly:

```wdl
command <<<
    msisensor-pro pro -d "~{msi_list}" -t "~{bam}" -o "~{sample_name}" ...
>>>
requirements {
    container: container
}
```

Development and production use the same backend, `slurm_apptainer`: sprocket submits every task
with `sbatch` and runs it with `apptainer exec` in its image, mounting its `File`/`Directory`
inputs and its work directory (nothing else of the host is visible).

| | Development | Production |
|---|---|---|
| Config | `backends.conf/sprocket.prod.toml` | `backends.conf/sprocket.prod.toml` |
| Slurm | single node on this workstation (`backends.conf/slurm-local/install.sh`) | the cluster, partition `prod`, 600 min |
| `*_container` (`params/containers.json`) | `file://.../apptainer-recipes/<tool>.sif` | `file:///bioinfo/apptainer_sif/<tool>_<version>.sif` |

The workflow's python scripts are deliberately **not** in the image: they are `File` inputs (`msi_main_helper`, `msi_report_helper`),
versioned with the WDL, so editing them needs no image rebuild (a frozen image with the scripts inside
only if the workflow is ever distributed without this folder).

Every task only reads its inputs and writes in its own directory (the run directory is a
`Directory` input, the reports are `Directory` outputs). The reports folder (`$VIZ`) is only
written by `index`, after the run.

## Files

| File | Role |
|---|---|
| `msi_detection-0.1.0.wdl` | The workflow |
| `module/msisensorpro-1.3.0.wdl`, `module/bedtools-2.31.1.wdl` | Tool modules (tasks `scan`, `pro`; `intersect`) |
| `main-scripts/msi_detection_pre_p.py` | Before the workflow: params JSON of the runs left to analyse (run by `BASH.sprocket.cron.sh`) |
| `main-scripts/msi_detection_main_p.py` | Reference list steps of the workflow (`msi_main_helper`) |
| `main-scripts/msi_detection_report_p.py` | Report step of the workflow (`msi_report_helper`): this run's pages (page icon embedded) |
| `main-scripts/msi_detection_post_p.py` | After the workflow: `index`, the all-runs dashboard (run by `BASH.sprocket.cron.sh` / `BASH.sprocket.batch.sh`) |
| `apptainer-recipes/python-3.12.def` | Generic python3 + awk image of the python/bash tasks |
| `params/containers.json` | Images and python helpers of the workflow (passed to every `sprocket run`) |
| `params/template.json` | Params JSON template (GRCh37, Juno v3 panel, ROI hg19): `pre_process` fills `run` per run |
| `backends.conf/sprocket.prod.toml` | sprocket config: Slurm + Apptainer |
| `backends.conf/slurm-local/` | Single-node Slurm of this workstation (`sudo ./install.sh`) |
| `BASH.sprocket.cron.sh`, `companion-scripts/BASH.sprocket.batch.sh` | `sprocket run` with this folder's config: cron mode (run found by pre_process) / batch mode (folder of params JSON) |
| `companion-scripts/make_configs.py`, `BASH.make_configs.sh` | Batch mode, before: `<run>_params.json` of the given runs of a BAM group folder (`failed-runs.txt`: runs without BAM); its Slurm/Apptainer launcher |
| `companion-scripts/export_tsv.py`, `BASH.export_tsv.sh` | Batch mode, after: results of a reports folder as a TSV, one line per sample; its Slurm/Apptainer launcher |
| `BASH.doc.sh`, `apptainer-recipes/sprocket-doc-0.31.0.def` | HTML documentation (`docs/`, generated) and its image |
| `diagram/msi_detection_metromap.drawio` / `.preview.png`, `diagram/BASH.render_png.sh` | Workflow diagram, its PNG export |

## HTML documentation

Generated by sprocket from the WDL's `meta` / `parameter_meta`, with this README as home page:

```sh
./BASH.doc.sh        # sprocket dev doc . --output docs --index-page README.md ... + diagram copy
```

Open [`docs/index.html`](docs/index.html). The search bar needs the pages to be served over http (browsers block its
index files from `file://`):

```sh
python3 -m http.server -d docs 8000      # then http://localhost:8000
```

`BASH.doc.sh` runs `sprocket dev doc` in `apptainer-recipes/sprocket-doc-0.31.0.sif` (sprocket +
Node.js + `pagefind`, pre-cached for the search index): nothing to install on the host, no download.
Rebuild the image with:

```sh
cd apptainer-recipes && apptainer build --fakeroot sprocket-doc-0.31.0.sif sprocket-doc-0.31.0.def
```


