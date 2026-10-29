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
0 */4 * * 4,5 /path/to/pipeline-jp/BASH.msi_detection.cron.sh >> /path/to/BASH.msi_detection.cron.log 2>&1
```


[![msi_detection metro map](docs/diagram/msi_detection_metromap.preview.png)](docs/diagram/msi_detection_metromap.preview.png)


## Contents
- [Introduction](#introduction)
- [Configuration](#configuration)
- [The cron steps](#the-cron-steps)
- [Tasks](#tasks)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Companion scripts](#companion-scripts)


## Configuration


**Fill the MSI locker with the runs ADIVAR has already finished.** The cron then only analyses the
runs that are not in the MSI locker.

```sh
# same lockers as LOCKER_ADIVAR / LOCKER_MSI of BASH.msi_detection.cron.sh
L=/path/to/locker
./companion-scripts/make_msi_locker.py --locker-adivar $L/adivar-testNS.json --locker-msi $L/msi-sensor.locker.txt --dry-run   # prints the runs to add
./companion-scripts/make_msi_locker.py --locker-adivar $L/adivar-testNS.json --locker-msi $L/msi-sensor.locker.txt
```

**Edit these files by hand at install:**

- [`params/template.json`](params/template.json): references and thresholds of each run
- [`params/install.json`](params/install.json): Apptainer images (`/media/jp/SSD2/CHU/containers/*.sif`, the `CONTAINERS_DIR` of the shell scripts; their recipes stay in `apptainer-recipes/*.def`) and Python scripts
- [`backends.conf/sprocket.cron.toml`](backends.conf/sprocket.cron.toml), [`backends.conf/sprocket.batch.toml`](backends.conf/sprocket.batch.toml): sprocket and Slurm settings of each mode (same settings, one call cache each)
- [`BASH.msi_detection.cron.sh`](BASH.msi_detection.cron.sh): folders and lockers of the cron mode
- finally, edit the crontab, to launch `BASH.msi_detection.cron.sh` 

`params/template.json`: Define a FASTA genome (`genome`), a BED of your panel (`bed_panel`) and a
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

**`BASH.msi_detection.cron.sh`** (cron mode)

| Variable | Role |
|---|---|
| `PIPELINE_DIR` | this folder |
| `CONTAINERS_DIR` | Apptainer images (`*.sif`) |
| `SCRATCH_DIR` | scratch folder, where the script runs and everything but the reports is computed: sprocket's runs (`Sprocket/`), the params JSON (`Config/`), the copies of ADIVAR's locker and of the run's sample sheet (`Temp/`) and the call caches (`cache/cron`, `cache/batch`) |
| `RESULTS_DIR` | final reports folder: `index` copies each run's reports into it and rebuilds the dashboard |
| `RUNS_DIR` | parent folder of the sequencing run folders: the cron copies `$RUNS_DIR/<run>/bcl2fastq/SampleSheet.csv` (`<run>`: name in ADIVAR's locker) into `$SCRATCH_DIR/Temp/<run>-SampleSheet.csv` and reads only that copy |
| `BAMS_DIR` | parent folder of the BAM folders, each named after an Experiment Name: the run analysed is `$BAMS_DIR/<Experiment Name>/` |
| `LOCKER_ADIVAR` | ADIVAR's locker (read only): copied at each launch into `$SCRATCH_DIR/Temp/`, the only file the containers read, so that the three jobs see the same state |
| `LOCKER_MSI` | MSI locker (runs already taken), created if missing |



## The cron steps


`BASH.msi_detection.cron.sh` (cron) first runs `msi_detection_pre_p.py` three times, each container
mounting only what that step reads:

1. `--print-run` (job `msi-pick_run`; the copy of ADIVAR's locker and the MSI locker): prints the run finished by ADIVAR and not yet
   in the MSI locker;
2. `--print-experiment` (job `msi-experiment`; + the copy of that run's sample sheet, `$SCRATCH_DIR/Temp/<run>-SampleSheet.csv`, given as `--samplesheet`): prints
   the Experiment Name of its `[Header]`, the name of its BAM folder;
3. without option (job `msi-pre_process`; + that run's `$BAMS_DIR/<Experiment Name>/`, the MSI
   locker and `$SCRATCH_DIR/Config/` writable): writes its params JSON in `$SCRATCH_DIR/Config/` and adds the run
   to the locker.

`$RUNS_DIR` is never mounted (the sample sheet is copied by the script, on the host), nor the whole `$BAMS_DIR`. A sample sheet without
Experiment Name or a missing BAM folder stops the script before the run is locked (retried at the next launch).
Then both scripts, for each params JSON:

1. sets the destinations `RESULTS_DIR` (final reports folder) and `SCRATCH_DIR` (sprocket's runs), created if needed:
   hardcoded in `BASH.msi_detection.cron.sh`.
2. runs `sprocket run -s -c $SPROCKET_CONFIG` (`backends.conf/sprocket.cron.toml` / `sprocket.batch.toml`) in
   `$SCRATCH_DIR/Sprocket/runs/msi_detection/<timestamp>/`, one folder per task;
3. runs `msi_detection_post_p.py index` as a Slurm job (`srun`, in `python-MSI-3.12.sif`): copies the run's
   `viz/runs/<run_id>/` into `$RESULTS_DIR/runs/` and rebuilds the dashboard from every
   `$RESULTS_DIR/runs/*/` (`index_msi_run.js`, `index_msi_sample.js`; the two `.html` are written once,
   only rewritten if their template changes).

The workflow itself never writes outside its tasks: the dashboard covers every run, a state that
changes with each run, so it is built after the run (`index`), not by a task (a task reading the
current reports would not be reproducible nor cacheable). `index` can be re-run at any time, e.g.
after deleting a run folder. Everything else stays in sprocket's run folder (paths in its
`outputs.json`), every task keeping its `work/` folder (~300 MB per run).


## Tasks

Calls are named after their line on the diagram, `t<line><rank>_<step>`: `ta` BAMs (green),
`tb` msisensor-pro list (blue, panel steps included), `tc` ROI markers (red).
CPU and Memory are what each job reserves from Slurm (`requirements` of the task, passed to
`sbatch`); only `ta02_sensorpro_pro` takes `params.threads`, every other task keeps its default
`threads = 1`. No time limit is set: the partition's applies.

| Call (task) | Runs | What it does | Main output | CPU | Memory | Set in |
|---|---|---|---|---|---|---|
| `ta01_find_bams` (`find_bams`) | python image, bash | Lists every `*.bam` of the run directory (`inputs.run`) (sub-folders with `params.recursive`), with its index (`.bam.bai` or `.bai`), its library (folder relative to `input_dir`, `.` if none) and its sample name (BAM basename). `inputs.run` is a `Directory` (mounted by the engine); paths are returned relative to it. With `inputs.samples`, keeps only the BAMs whose name contains one of these IDs (every library). Fails if a BAM has no index, no BAM is found, or a sample ID matches no BAM. | `bams`: one row per BAM | 1 | 512 MiB | `msi_detection-1.0.0.wdl` |
| `tb01_sensorpro_scan` (`sensorpro.scan`) | msisensor-pro | `msisensor-pro scan` of the genome: every microsatellite site. Skipped when `references.pro_list` is given. Not redone from one run to the next: sprocket's call cache (`[run.task] cache = "on"`) reuses the previous scan (same genome and image), as long as that run's folder is kept. | `<genome>.msisensor_pro.list` | 1 | 4 GiB | `module/msisensorpro-1.3.0.wdl` |
| `tb02_list_to_bed` (`list_to_bed`) | python image | Converts the list to BED6 (start = location − 1, end = start + unit length × repeat times). | `<list>.bed` | 1 | 2 GiB | `msi_detection-1.0.0.wdl` |
| `tb03_panel_chrom` (`align_chrom_naming`) | python image, awk | Only with `bed_panel`. Copy of the panel BED with its chromosome names (`chr1` / `1`) aligned on the list's: bedtools matches them literally. | `<panel>.chrnorm.bed` | 1 | 512 MiB | `msi_detection-1.0.0.wdl` |
| `tb04_panel_intersect` (`bedtools.intersect`) | bedtools | `bedtools intersect -wa -wb -f 1.0`: the sites fully inside a panel region, one row per overlap. | `<list>.panel.intersect.tsv` | 1 | 2 GiB | `module/bedtools-2.31.1.wdl` |
| `tb05_panel_correspondence` (`name_correspondence`) | python image, awk | Turns the intersection into a BED6 named `<site>#<probe>#<chrom:start-end>`. | `<list>.panel.bed` | 1 | 512 MiB | `msi_detection-1.0.0.wdl` |
| `tb06_mark_list` (`mark_list`) | python image | Only with `bed_panel`. Full list, same rows, `filter` = `PASS` for the panel sites, `NO` elsewhere: **the list analysed by msisensor-pro**. | `<list>.panel.list` | 1 | 2 GiB | `msi_detection-1.0.0.wdl` |
| `tc01_roi_chrom` (`align_chrom_naming`) | python image, awk | Only with `bed_roi`. Same as `tb03` with the ROI BED. | `<roi>.chrnorm.bed` | 1 | 512 MiB | `msi_detection-1.0.0.wdl` |
| `tc02_roi_intersect` (`bedtools.intersect`) | bedtools | Same as `tb04` with the ROI markers, on the panel sites when there is a panel (genome-wide otherwise). | `<list>.panel.roi.intersect.tsv` | 1 | 2 GiB | `module/bedtools-2.31.1.wdl` |
| `tc03_roi_correspondence` (`name_correspondence`) | python image, awk | Same as `tb05`; name ends with `#<marker>`. | `<list>.panel.roi.bed` | 1 | 512 MiB | `msi_detection-1.0.0.wdl` |
| `ta02_sensorpro_pro` × N (`sensorpro.pro`) | msisensor-pro | Scatter, one per BAM: `msisensor-pro pro -d <analysis list> -t <bam> -b threads -c coverage_threshold`. A failure on any BAM fails the workflow. | `<sample>`, `_all`, `_dis`, `_unstable` | `params.threads` (8: `pro -b`) | 4 GiB | CPU: params JSON; memory: `module/msisensorpro-1.3.0.wdl` |
| `ta03_report` (`report`) | python image | Reads every BAM's `_all` (`File`s, nothing copied) and builds this run's reports in `viz/runs/<run_id>/` (run page, sample pages; the all-runs dashboard is built by `index` after the run); log on its stderr; writes nothing outside its task. | `viz/` (`Directory`) | 1 | 2 GiB | `msi_detection-1.0.0.wdl` |


## Inputs

Keys of the params JSON, all prefixed `msi_detection.` (see [`params/template.json`](params/template.json)).
`inputs`, `references` and `params` are objects (WDL structs); 

| Key | Required | Description |
|---|---|---|
| `inputs.run` | yes | Run directory; its name is the `run_id` |
| `inputs.samples` | no | Sample IDs: only the BAMs named `<ID>.bam` (exact name, every library) are analysed (an ID with no BAM fails the run). Without it: every BAM but the controls |
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

Two main folders, `RESULTS_DIR` and `SCRATCH_DIR`:

- `RESULTS_DIR`: the HTML reports and the all-runs dashboard, the folder to open (described first below);
- `SCRATCH_DIR`: sprocket's runs, with every task's files and logs (tree at the end of this section),
  and the call caches, `cache_dir` of `backends.conf/sprocket.cron.toml` (cron, `cache/cron`) and
  `sprocket.batch.toml` (batch, `cache/batch`): one cache per mode.

`ta03_report` returns `viz/` as a `Directory` output, in its folder of the sprocket run:

```text
$SCRATCH_DIR/Sprocket/runs/msi_detection/<timestamp>/calls/ta03_report/attempts/0/work/viz/
```

`index` then copies its `runs/<run_id>/` to `$RESULTS_DIR`, adds the run's msisensor-pro `_all` files to it, and
builds the dashboard there:

```text
<RESULTS_DIR>/
├── index_msi_run.html              all runs (written once)
├── index_msi_run.js                its data, rewritten by each index
├── index_msi_sample.html           all samples of all runs (written once)
├── index_msi_sample.js             its data, rewritten by each index
└── runs/
    ├── Run154_Juno/                one folder per run (<run_id>)
    │   ├── Run154_Juno_index_run.html     
    │   ├── Run154_Juno_index_run.js
    │   ├── NICK_index_sample.html    
    │   ├── NICK_index_sample.js
    │   └── Raw/                    every BAM's msisensor-pro _all, <library>_<sample>_all
    │       ├── libA_NICK_all
    │       └── libB_NICK_all
    └── Run171_Juno/
        ├── Run171_Juno_index_run.html
        ├── Run171_Juno_index_run.js
        ├── JIM_index_sample.html
        ├── JIM_index_sample.js
        ├── BOB_index_sample.html
        └── BOB_index_sample.js
```

Of the raw msisensor-pro outputs, only the `_all` files are copied (in `runs/<run_id>/Raw/`, as
`<library>_<sample>_all`, the same sample name being in every library): the others (`_unstable`, `_dis`) stay in each
`ta02_sensorpro_pro` folder of the sprocket run (`_all` and `_unstable` listed in `outputs.json`).

When every call came from the cache, `outputs.json` points to
the run that computed them.

sprocket's work folder, in `SCRATCH_DIR`:

```text
$SCRATCH_DIR/
└── Sprocket/
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

## Companion scripts

Helper scripts around the workflow, in `companion-scripts/`. The `BASH.*.sh` launchers run their
python script as a Slurm job (`srun`) in `python-MSI-3.12.sif`;

| Script | Mode | Inputs | What it does |
|---|---|---|---|
| [`make_msi_locker.py`](companion-scripts/make_msi_locker.py) | cron, before the first launch | `--locker-adivar` ADIVAR's locker JSON, `--locker-msi` MSI locker (text, created if missing), `--dry-run` | Adds every run already finished by ADIVAR to the MSI locker |
| [`BASH.make_configs.sh`](companion-scripts/BASH.make_configs.sh) → [`make_configs.py`](companion-scripts/make_configs.py) | batch, before | `BAMS_DIR` group folder of BAM folders, `RUNS_FILE` runs file (run, or run + sample ID per line), `OUTPUT_DIR` | Writes one `<run>_params.json` per run of a runs file (optionally restricted to some samples); runs without BAM go to `failed-runs.txt` |
| [`BASH.msi_detection.batch.sh`](companion-scripts/BASH.msi_detection.batch.sh) | batch | `CONFIGS_DIR` folder of params JSON, `RESULTS_DIR` final reports folder (sprocket's runs: `$SCRATCH_DIR/Sprocket/`) | Runs the workflow on every params JSON of a folder, then rebuilds the dashboard |
| [`BASH.export_tsv.sh`](companion-scripts/BASH.export_tsv.sh) → [`export_tsv.py`](companion-scripts/export_tsv.py) | batch, after | `RESULTS_DIR` reports folder (read only), `OUTPUT_DIR` where the TSV goes | Exports a reports folder as a TSV, one line per sample (`results_msi_<date>.tsv`) |

```sh
./companion-scripts/make_msi_locker.py --locker-adivar ADIVAR_JSON --locker-msi MSI_TXT [--dry-run]
./companion-scripts/BASH.make_configs.sh BAMS_DIR RUNS_FILE OUTPUT_DIR
./companion-scripts/BASH.msi_detection.batch.sh CONFIGS_DIR RESULTS_DIR  # in screen / tmux
./companion-scripts/BASH.export_tsv.sh RESULTS_DIR OUTPUT_DIR
```

**Note:** `BASH.msi_detection.batch.sh` runs sprocket with its own config,
[`backends.conf/sprocket.batch.toml`](backends.conf/sprocket.batch.toml), not the cron's
`sprocket.cron.toml`: same settings, but a different `cache_dir`, so that the batch and cron call
caches stay apart (a batch never reuses a cron result, nor the reverse).

