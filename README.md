#  Microsatellite instability (MSI) detection Workflow

## Introduction

This workflow detects microsatellite instability (MSI). It is written in [WDL](https://openwdl.org/) and run with
[sprocket](https://sprocket.bio) on a Slurm + Apptainer backend.
It uses two tool modules, MSIsensor-pro and bedtools, and three Apptainer images: one per tool plus
a generic Python image for the other steps.

It runs MSIsensor-pro (tumor-only) on every BAM of a run directory, once the upstream pipeline
(ADIVAR) has finished that run. Controls are skipped: BAMs whose name starts with `TN` or `TPH`.

Before that, `msisensor-pro scan` lists every microsatellite of the genome (once: the result is
cached, or given as `pro_list`). bedtools then intersects this list with the panel BED, so that
only the microsatellites inside the panel regions are analysed, and with the ROI BED, which names
the classic MSI markers (BAT25, BAT26, NR21...) among them.


In cron mode, the pipeline compares two lockers: ADIVAR's (JSON, `status` 3 = finished) and the MSI
locker (text, one run per line = already analysed). A run finished by ADIVAR and not yet in the MSI
locker is analysed, then added to it.

As a safeguard, the cron only launches when exactly one run
differs between the two lockers: with more, it stops with an error and analyses nothing.

For the picked run, the pipeline reads a single field of its `SampleSheet.csv`: the Experiment Name, which gives the BAM folder to analyse.

<table>
<tr>
<th>ADIVAR locker</th>
<th>MSI locker (built from it)</th>
</tr>
<tr>
<td><pre><code>
{
    "261001_XV501631_0734_ATESTR141": {
        "status": 3,
        "last-status": 3
    },
    "261002_XV501631_0542_ATESTR145": {
        "status": 2,
        "last-status": 2
    }
}
</code></pre></td>
<td><pre><code>
261001_XV501631_0734_ATESTR141
</code></pre>
(ATESTR145, still running, is left out: the cron analyses it once ADIVAR finishes it.)</td>
</tr>
</table>

Cron mode, e.g. every 4 hours on Thursdays and Fridays:

```sh
0 */4 * * 4,5 /path/to/pipeline-jp/BASH.sprocket.cron.sh >> /path/to/BASH.sprocket.cron.log 2>&1
```


[![msi_detection metro map](docs/diagram/msi_detection_metromap.preview.png)](docs/diagram/msi_detection_metromap.preview.png)


## Contents
- [Introduction](#introduction)
- [Before the first cron launch](#before-the-first-cron-launch)
- [Configuration](#configuration)
- [Quick start](#quick-start)
  - [Batch mode: params JSON and TSV export](#batch-mode-params-json-and-tsv-export)
- [Tasks](#tasks)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Files](#files)
- [Companion scripts](#companion-scripts)

## Before the first cron launch

So that the cron does not analyse every past run, fill the MSI locker
with the runs ADIVAR has already finished.

```sh
# same lockers as LOCKER_ADIVAR / LOCKER_MSI of BASH.sprocket.cron.sh
L=/path/to/locker
./companion-scripts/make_msi_locker.py --locker-adivar $L/adivar-testNS.json --locker-msi $L/msi-sensor.locker.txt --dry-run   # prints the runs to add
./companion-scripts/make_msi_locker.py --locker-adivar $L/adivar-testNS.json --locker-msi $L/msi-sensor.locker.txt
```

## Configuration

Edit these files by hand at install:

- `params/template.json`: references and thresholds of each run
- `params/install.json`: Apptainer images and Python scripts
- `backends.conf/sprocket.prod.toml`: sprocket and Slurm settings
- `BASH.sprocket.cron.sh`: folders and lockers of the cron mode
- finally, edit the crontab, to launch `BASH.sprocket.cron.sh` 

`params/template.json`: (copied by `pre_process` for each run; `run` and `samples` stay empty, filled per run). Define a FASTA genome (`genome`), a BED of your panel (`bed_panel`) and a
BED of the ROI markers (`bed_roi`), the marker name in the 4th column (custom naming), e.g.
`ROI.hg19.bed`:

```text
chr2	47641530	47641644	BAT26	0	+
chr2	95849332	95849414	NR24	0	+
chr4	55598161	55598266	BAT25	0	+
chr11	102193479	102193564	NR27	0	+
chr11	125490736	125490816	NR22	0	+
chr13	31722591	31722667	HT17	0	+
chr14	23652317	23652397	NR21	0	+
```

Set the workflow's thresholds in the same file, under `msi_detection.params` (current values in
[`params/template.json`](params/template.json), each key described in [Inputs](#inputs)).

[`params/install.json`](params/install.json): paths of the Apptainer images and Python scripts the WDL workflow uses
(details below).

**`BASH.sprocket.cron.sh`** (cron mode)

| Variable | Role |
|---|---|
| `WDL_DIR` | this folder |
| `VIZ` | reports and dashboard (written by `index`) |
| `WORK` | sprocket's runs (`sprocket/`) and the params JSON (`config/`) |
| `RUNS_DIR` | parent folder of the sequencing run folders: the cron reads `$RUNS_DIR/<run>/SampleSheet.csv` (`<run>`: name in ADIVAR's locker) |
| `BAMS_DIR` | parent folder of the BAM folders, each named after an Experiment Name: the run analysed is `$BAMS_DIR/<Experiment Name>/` |
| `LOCKER_ADIVAR` | ADIVAR's locker (read only) |
| `LOCKER_MSI` | MSI locker (runs already taken), created if missing |


**`companion-scripts/BASH.sprocket.batch.sh`** (batch mode): `WDL_DIR` and the same `srun` line; `VIZ`
and `WORK`, the same two folders as above, are its required arguments.

**`backends.conf/sprocket.prod.toml`**

| Key | Value | Role |
|---|---|---|
| `[run.backends.slurm.sbatch] args` | `["--partition=prod"]` | partition (and any other `sbatch` option) of the tasks |
| `[run.task] cache_dir` | a folder | sprocket's call cache (shared by all the runs) |
| `[run.backends.slurm] interval` | `5` | seconds between two `sacct` polls: short for local tests, raise it on the cluster |
| `max_concurrency`, `[run.workflow.scatter] concurrency` | `50`, `10` | jobs submitted at once, scatter elements in parallel |

**`params/template.json`** (copied by `pre_process` for each run; `run` and `samples` stay empty, filled per run)

| Key | Role |
|---|---|
| `msi_detection.references.genome` | reference FASTA (e.g. GRCh37) |
| `msi_detection.references.bed_panel` | BED of the sequencing panel |
| `msi_detection.references.bed_roi` | BED of the ROI markers |
| `msi_detection.params` | `threads` (CPUs of each `msisensor-pro pro` only), `coverage_threshold`, `p_process` thresholds |

**`params/install.json`**: what depends on the installation, not on the run: the workflow's
images and python helpers (scripts given to the tasks as `File`, not built into the image), the
only place they are set (absolute paths; images as `file://` + path). These five WDL inputs have no default: both
scripts pass this file before each run's params JSON (`sprocket run ... @install.json @<params>.json`;
for a key in both, the last file wins, so a params JSON can still override an image for its run).
Forgetting it stops sprocket before anything runs (`missing required input`).

| Key | Points to |
|---|---|
| `msi_detection.msisensorpro_container` | `apptainer-recipes/msisensor-pro-1.3.0.sif` |
| `msi_detection.bedtools_container` | `apptainer-recipes/bedtools-2.31.1.sif` |
| `msi_detection.python_container` | `apptainer-recipes/python-3.12.sif` |
| `msi_detection.msi_main_helper` | `main-scripts/msi_detection_main_p.py` |
| `msi_detection.msi_report_helper` | `main-scripts/msi_detection_report_p.py` |

The images of the scripts' own Slurm jobs follow `WDL_DIR`: `$WDL_DIR/apptainer-recipes/python-3.12.sif`
(`python_slurm_job` of `BASH.sprocket.cron.sh` / `BASH.sprocket.batch.sh`: pick_run, experiment,
pre_process, index) and `$WDL_DIR/apptainer-recipes/sprocket-doc-0.31.0.sif` (`BASH.doc.sh`).

**crontab** (not in this folder), for example every 30 min with a log:

```sh
*/30 * * * * /path/to/pipeline-jp/BASH.sprocket.cron.sh >> /path/to/cron.log 2>&1
```


## Quick start

```sh
cd /path/to/pipeline-jp
./BASH.sprocket.cron.sh                       # cron mode: pre_process picks the run left to analyse, then runs it
./companion-scripts/BASH.sprocket.batch.sh CONFIG_DIR VIZ WORK   # batch mode: every params JSON of CONFIG_DIR
```

The batch mode runs `sprocket` in the current shell (hours for a whole group): launch it in
`screen` / `tmux`, or a lost terminal / SSH session stops it. Not needed for the cron mode (no terminal).

`BASH.sprocket.cron.sh` (cron) first runs `msi_detection_pre_p.py` three times, each container
mounting only what that step reads:

1. `--print-run` (job `msi-pick_run`; the two lockers): prints the run finished by ADIVAR and not yet
   in the MSI locker;
2. `--print-experiment` (job `msi-experiment`; + that run's `$RUNS_DIR/<run>/SampleSheet.csv`): prints
   the Experiment Name of its `[Header]`, the name of its BAM folder;
3. without option (job `msi-pre_process`; + that run's `$BAMS_DIR/<Experiment Name>/`, the MSI
   locker and `$WORK/config/` writable): writes its params JSON in `$WORK/config/` and adds the run
   to the locker.

Neither the whole `$RUNS_DIR` nor the whole `$BAMS_DIR` is mounted. A sample sheet without
Experiment Name or a missing BAM folder stops the script before the run is locked (retried at the next launch).
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
S=/path/to/data
./companion-scripts/BASH.make_configs.sh $S/input/foie-mss $S/config/runs-test.txt $S/config/config-test
./companion-scripts/BASH.sprocket.batch.sh $S/config/config-test $S/output/config-test/viz $S/output/config-test/work
./companion-scripts/BASH.export_tsv.sh $S/output/config-test/viz $S/output/config-test      # -> $S/output/config-test/results_msi_<date>.tsv
```

The python scripts can also be run directly on the host (only the standard library):
`make_configs.py BAMS_DIR [RUN ...] [--runs-file TXT] --output-dir DIR`,
`export_tsv.py --viz VIZ_DIR [--output-dir DIR | --output TSV]`.

Run by hand, with a params JSON written by `pre_process` (or a filled copy of `params/template.json`):

```sh
sprocket run -s -c backends.conf/sprocket.prod.toml msi_detection-0.1.0.wdl @params/install.json @$WORK/config/<date>-<run>-<Experiment Name>.json
```

## Tasks

Calls are named after their line on the diagram, `t<line><rank>_<step>`: `ta` BAMs (green),
`tb` msisensor-pro list (blue, panel steps included), `tc` ROI markers (red).
CPU and Memory are what each job reserves from Slurm (`requirements` of the task, passed to
`sbatch`); only `ta02_sensorpro_pro` takes `params.threads`, every other task keeps its default
`threads = 1`. No time limit is set: the partition's applies.

| Call (task) | Runs | What it does | Main output | CPU | Memory | Set in |
|---|---|---|---|---|---|---|
| `ta01_find_bams` (`find_bams`) | python image, bash | Lists every `*.bam` of the run directory (`inputs.run`) (sub-folders with `params.recursive`), with its index (`.bam.bai` or `.bai`), its library (folder relative to `input_dir`, `.` if none) and its sample name (BAM basename). `inputs.run` is a `Directory` (mounted by the engine); paths are returned relative to it. With `inputs.samples`, keeps only the BAMs whose name contains one of these IDs (every library). Fails if a BAM has no index, no BAM is found, or a sample ID matches no BAM. | `bams`: one row per BAM | 1 | 512 MiB | `msi_detection-0.1.0.wdl` |
| `tb01_sensorpro_scan` (`sensorpro.scan`) | msisensor-pro | `msisensor-pro scan` of the genome: every microsatellite site. Skipped when `references.pro_list` is given. Not redone from one run to the next: sprocket's call cache (`[run.task] cache = "on"`) reuses the previous scan (same genome and image), as long as that run's folder is kept. | `<genome>.msisensor_pro.list` | 1 | 4 GiB | `module/msisensorpro-1.3.0.wdl` |
| `tb02_list_to_bed` (`list_to_bed`) | python image | Converts the list to BED6 (start = location − 1, end = start + unit length × repeat times). | `<list>.bed` | 1 | 2 GiB | `msi_detection-0.1.0.wdl` |
| `tb03_panel_chrom` (`align_chrom_naming`) | python image, awk | Only with `bed_panel`. Copy of the panel BED with its chromosome names (`chr1` / `1`) aligned on the list's: bedtools matches them literally. | `<panel>.chrnorm.bed` | 1 | 512 MiB | `msi_detection-0.1.0.wdl` |
| `tb04_panel_intersect` (`bedtools.intersect`) | bedtools | `bedtools intersect -wa -wb -f 1.0`: the sites fully inside a panel region, one row per overlap. | `<list>.panel.intersect.tsv` | 1 | 2 GiB | `module/bedtools-2.31.1.wdl` |
| `tb05_panel_correspondence` (`name_correspondence`) | python image, awk | Turns the intersection into a BED6 named `<site>#<probe>#<chrom:start-end>`. | `<list>.panel.bed` | 1 | 512 MiB | `msi_detection-0.1.0.wdl` |
| `tb06_mark_list` (`mark_list`) | python image | Only with `bed_panel`. Full list, same rows, `filter` = `PASS` for the panel sites, `NO` elsewhere: **the list analysed by msisensor-pro**. | `<list>.panel.list` | 1 | 2 GiB | `msi_detection-0.1.0.wdl` |
| `tc01_roi_chrom` (`align_chrom_naming`) | python image, awk | Only with `bed_roi`. Same as `tb03` with the ROI BED. | `<roi>.chrnorm.bed` | 1 | 512 MiB | `msi_detection-0.1.0.wdl` |
| `tc02_roi_intersect` (`bedtools.intersect`) | bedtools | Same as `tb04` with the ROI markers, on the panel sites when there is a panel (genome-wide otherwise). | `<list>.panel.roi.intersect.tsv` | 1 | 2 GiB | `module/bedtools-2.31.1.wdl` |
| `tc03_roi_correspondence` (`name_correspondence`) | python image, awk | Same as `tb05`; name ends with `#<marker>`. | `<list>.panel.roi.bed` | 1 | 512 MiB | `msi_detection-0.1.0.wdl` |
| `ta02_sensorpro_pro` × N (`sensorpro.pro`) | msisensor-pro | Scatter, one per BAM: `msisensor-pro pro -d <analysis list> -t <bam> -b threads -c coverage_threshold`. A failure on any BAM fails the workflow. | `<sample>`, `_all`, `_dis`, `_unstable` | `params.threads` (8: `pro -b`) | 4 GiB | CPU: params JSON; memory: `module/msisensorpro-1.3.0.wdl` |
| `ta03_report` (`report`) | python image | Reads every BAM's `_all` (`File`s, nothing copied) and builds this run's reports in `viz/runs/<run_id>/` (run page, sample pages; the all-runs dashboard is built by `index` after the run); log on its stderr; writes nothing outside its task. | `viz/` (`Directory`) | 1 | 2 GiB | `msi_detection-0.1.0.wdl` |

The python steps of the workflow are subcommands of two self-contained scripts (a WDL task gets
its script alone): [`main-scripts/msi_detection_main_p.py`](main-scripts/msi_detection_main_p.py)
(`msi_main_helper`: `list-to-bed`, `mark-list`, the reference list) and
[`main-scripts/msi_detection_report_p.py`](main-scripts/msi_detection_report_p.py) (`msi_report_helper`:
`report`). They are kept apart because sprocket's call cache keys on each task's input files:
editing the report must not invalidate the list steps, whose outputs are whole-genome lists
(~100 MB each). [`main-scripts/msi_detection_post_p.py`](main-scripts/msi_detection_post_p.py) only holds
`index`, run after the workflow by the BASH scripts; it imports the page code (CSS, `<head>`,
`.js` writer) of `msi_detection_report_p.py`, next to it in `main-scripts/`.

## Inputs

Keys of the params JSON, all prefixed `msi_detection.` (see [`params/template.json`](params/template.json)).
`inputs`, `references` and `params` are objects (WDL structs); 

```json
"msi_detection.inputs":     {"run": "...", "samples": ["E25D04852T", "..."]},
"msi_detection.references": {"genome": "...", "bed_panel": "...", "bed_roi": "...", "pro_list": "..."},
"msi_detection.params":     {"threads": 8, "coverage_threshold": 1, "recursive": true, "p_process": {"min_reads": 200, "pval_cutoff": 0.1, "sample_na_fraction": 0.5, "consensus_roi_cutoff": 3}},
```


| Key | Required | Description |
|---|---|---|
| `inputs.run` | yes | Run directory; its name is the `run_id` |
| `inputs.samples` | no | Sample IDs: only the BAMs whose name contains one of them are analysed (an ID with no BAM fails the run). Without it: every BAM but the controls |
| `references.genome` | yes | Fasta the BAMs are aligned on |
| `references.pro_list` | no | Existing scan list: skips the scan |
| `references.bed_panel` | no | Plain panel BED (no `track` line, e.g. `*.amplicons.bed`) |
| `references.bed_roi` | no | MSI markers BED (e.g. `ROI.hg19.bed`) |
| `params.threads` | yes | `pro -b` and the CPUs reserved by each `ta02_sensorpro_pro` (only task using it) |
| `params.coverage_threshold` | yes | `pro -c` |
| `params.recursive` | yes | Also look for BAMs in sub-folders (one per library) |
| `params.p_process.min_reads` | yes | Reports: a site counts (% unstable, site status) only with ≥ this many reads (e.g. 200) |
| `params.p_process.pval_cutoff` | yes | Reports: a site is unstable when its p-value is above this (e.g. 0.1) |
| `params.p_process.sample_na_fraction` | yes | Reports: sample status NA when more than this fraction of the loci listed on its sample page are NA or MSU (e.g. 0.5): ROI markers, and sites with a p-value in at least one library (loci without any, e.g. 0 reads everywhere, are left out) |
| `params.p_process.consensus_roi_cutoff` | yes | Reports: sample status MSI when at least this many ROI markers are unstable in consensus, MSS otherwise (e.g. 3) |


The genome, panel and ROI must use the same build as the BAMs (e.g. GRCh37 + `ROI.hg19.bed`).

## Outputs

Two main folders, `VIZ` and `WORK`:

- `VIZ`: the HTML reports and the all-runs dashboard, the folder to open (described first below);
- `WORK`: sprocket's runs, with every task's files and logs (tree at the end of this section).
  The call cache is elsewhere, in `cache_dir` of `backends.conf/sprocket.prod.toml`.

`ta03_report` returns `viz/` as a `Directory` output, in its folder of the sprocket run:

```text
$WORK/sprocket/runs/msi_detection/<timestamp>/calls/ta03_report/attempts/0/work/viz/
```

`index` then copies its `runs/<run_id>/` to `$VIZ` and
builds the dashboard there:

```text
<VIZ>/
├── index_msi_run.html              all runs (written once)
├── index_msi_run.js                its data, rewritten by each index
├── index_msi_sample.html           all samples of all runs (written once)
├── index_msi_sample.js             its data, rewritten by each index
└── runs/
    ├── Run154_Juno/                one folder per run (<run_id>)
    │   ├── Run154_Juno_index_run.html     
    │   ├── Run154_Juno_index_run.js
    │   ├── NICK_index_sample.html    
    │   └── NICK_index_sample.js
    └── Run171_Juno/
        ├── Run171_Juno_index_run.html
        ├── Run171_Juno_index_run.js
        ├── JIM_index_sample.html
        ├── JIM_index_sample.js
        ├── BOB_index_sample.html
        └── BOB_index_sample.js
```

The raw msisensor-pro outputs are not copied: they stay in each `ta02_sensorpro_pro` folder of the
sprocket run (`_all` and `_unstable` listed in `outputs.json`); 

When every call came from the cache, `outputs.json` points to
the run that computed them, possibly in another `WORK`.

sprocket's work folder, `WORK`:

```text
$WORK/
└── sprocket/
    ├── .lock                                   sprocket's lock
    ├── sprocket.db                             sprocket's database (runs, index)
    ├── index/
    │   └── msi_detection/
    │       └── Run154_Juno_params/             one folder per params JSON (--index-on)
    │           ├── outputs.json  -> ../../../runs/msi_detection/<timestamp>/outputs.json
    │           ├── viz           -> .../calls/ta03_report/attempts/0/work/viz
    │           ├── NICK_all, NICK_unstable     -> .../calls/ta02_sensorpro_pro-<N>/attempts/0/work/
    │           └── <genome>.msisensor_pro.panel{.bed,.list,.roi.bed}  -> .../calls/tb05|tb06|tc03_*/attempts/0/work/
    └── runs/
        └── msi_detection/
            ├── _latest -> <last timestamp>
            └── 2026-10-27_081009457281388/     one folder per sprocket run
                ├── inputs.json                 resolved inputs
                ├── outputs.json                workflow outputs
                ├── output.log                  sprocket's log
                └── calls/
                    ├── ta01_find_bams/
                    ├── tb01_sensorpro_scan/
                    ├── tb02_list_to_bed/ … tb06_mark_list/
                    ├── tc01_roi_chrom/ … tc03_roi_correspondence/
                    ├── ta02_sensorpro_pro-0/ … -<N>/      one per BAM (scatter); <sample>{,_all,_dis,_unstable} in attempts/0/work/
                    └── ta03_report/
                        ├── inputs.json, outputs.json
                        ├── tmp/
                        └── attempts/0/
                            ├── command                 script run
                            ├── stdout, stderr          task output (report log: stderr)
                            ├── sbatch_command, job_id, slurm.std{out,err}
                            ├── apptainer_command, apptainer.std{out,err}
                            └── work/                   working folder: viz/ (+ samples.tsv)
```

A call taken from the cache has no `attempts/`, only `inputs.json` and `outputs.json`: its results
are in the run that computed them.

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
| `params/install.json` | Images and python helpers of the workflow (passed to every `sprocket run`) |
| `params/template.json` | Params JSON template (GRCh37, Juno v3 panel, ROI hg19): `pre_process` fills `run` per run |
| `backends.conf/sprocket.prod.toml` | sprocket config: Slurm + Apptainer |
| `backends.conf/slurm-local/` | Single-node Slurm of this workstation (`sudo ./install.sh`) |
| `BASH.sprocket.cron.sh`, `companion-scripts/BASH.sprocket.batch.sh` | `sprocket run` with this folder's config: cron mode (run found by pre_process) / batch mode (folder of params JSON) |
| `companion-scripts/make_configs.py`, `BASH.make_configs.sh` | Batch mode, before: `<run>_params.json` of the given runs of a BAM group folder (`failed-runs.txt`: runs without BAM); its Slurm/Apptainer launcher |
| `companion-scripts/export_tsv.py`, `BASH.export_tsv.sh` | Batch mode, after: results of a reports folder as a TSV, one line per sample; its Slurm/Apptainer launcher |
| `companion-scripts/make_msi_locker.py` | Cron mode, before the first launch: adds every run finished by ADIVAR to the MSI locker |
| `BASH.doc.sh`, `apptainer-recipes/sprocket-doc-0.31.0.def` | HTML documentation (`docs/`, generated) and its image |
| `diagram/msi_detection_metromap.drawio` / `.preview.png`, `diagram/BASH.render_png.sh` | Workflow diagram, its PNG export |

## Companion scripts

Helper scripts around the workflow, in `companion-scripts/`. The `BASH.*.sh` launchers run their
python script as a Slurm job (`srun`) in `python-3.12.sif`; the python scripts only need the
standard library and can also be run directly on the host.

| Script | Mode | What it does |
|---|---|---|
| [`BASH.make_configs.sh`](companion-scripts/BASH.make_configs.sh) → [`make_configs.py`](companion-scripts/make_configs.py) | batch, before | Writes one `<run>_params.json` per run of a runs file (optionally restricted to some samples); runs without BAM go to `failed-runs.txt` |
| [`BASH.sprocket.batch.sh`](companion-scripts/BASH.sprocket.batch.sh) | batch | Runs the workflow on every params JSON of a folder, then rebuilds the dashboard |
| [`BASH.export_tsv.sh`](companion-scripts/BASH.export_tsv.sh) → [`export_tsv.py`](companion-scripts/export_tsv.py) | batch, after | Exports a reports folder as a TSV, one line per sample (`results_msi_<date>.tsv`) |
| [`make_msi_locker.py`](companion-scripts/make_msi_locker.py) | cron, before the first launch | Adds every run already finished by ADIVAR to the MSI locker |

```sh
./companion-scripts/BASH.make_configs.sh BAMS_DIR RUNS_FILE OUTPUT_DIR
./companion-scripts/BASH.sprocket.batch.sh CONFIG_DIR VIZ WORK         # in screen / tmux
./companion-scripts/BASH.export_tsv.sh VIZ OUTPUT_DIR
./companion-scripts/make_msi_locker.py --locker-adivar ADIVAR_JSON --locker-msi MSI_TXT [--dry-run]
```

Details: [Batch mode](#batch-mode-params-json-and-tsv-export) (runs file format, TSV columns) and
[Before the first cron launch](#before-the-first-cron-launch).
