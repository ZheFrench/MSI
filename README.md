#  Microsatellite instability (MSI) detection Workflow

## Introduction

This workflow detects microsatellite instability (MSI). It is written in [WDL](https://openwdl.org/) 1.1 and run with
[Cromwell](https://github.com/broadinstitute/cromwell) 92 on a Slurm + Apptainer backend.
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

**Install Cromwell.** It needs Java 17 or later (`JAVA` of the shell scripts,
`/usr/lib/jvm/java-17-openjdk-amd64/bin/java`); the jars are not in git:

```sh
mkdir -p soft/cromwell && cd soft/cromwell
curl -LO https://github.com/broadinstitute/cromwell/releases/download/92/cromwell-92.jar
curl -LO https://github.com/broadinstitute/cromwell/releases/download/92/womtool-92.jar
```

Nothing else is needed on the host but Slurm and Apptainer: every Python step, the launchers' own
included, runs in `python-MSI-3.12.sif`.

**Edit these files by hand at install:**

- [`params/template.json`](params/template.json): references and thresholds of each run
- [`params/install.json`](params/install.json): Apptainer images, named Docker-style `<name>:<tag>` (`msisensor-pro:1.3.0`, `bedtools:2.31.1`, `python-msi:3.12`): the backend runs `$CONTAINERS_DIR/<name>-<tag>.sif` (case ignored; recipes in `apptainer-recipes/*.def`), and Python scripts
- [`backends.conf/cromwell.conf`](backends.conf/cromwell.conf): Cromwell settings, the same for both modes: Slurm partition (`prod`), concurrent jobs, call cache
- [`BASH.msi_detection.cron.sh`](BASH.msi_detection.cron.sh): folders, lockers, `JAVA` and `CROMWELL_JAR` of the cron mode
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
| `REFDATA_DIR` | references (genome, panel and ROI BEDs, `pro_list`): bound read-only in the task containers with `PIPELINE_DIR`. On another filesystem than `SCRATCH_DIR`, Cromwell soft-links the task inputs instead of hard-linking them, and the link targets must be visible in the containers. Also set in the batch and test launchers |
| `SCRATCH_DIR` | scratch folder, where the script runs and everything but the reports is computed: Cromwell's runs and call cache (`Cromwell/`), the params JSON (`Config/`) and the copies of ADIVAR's locker and of the run's sample sheet (`Temp/`) |
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

1. sets the destinations `RESULTS_DIR` (final reports folder) and `SCRATCH_DIR` (Cromwell's runs), created if needed:
   hardcoded in `BASH.msi_detection.cron.sh`.
2. merges `params/install.json` and the params JSON into `$SCRATCH_DIR/Cromwell/inputs/<name>.json`
   (Cromwell takes a single inputs JSON; the params JSON wins);
3. runs `cromwell run` (`backends.conf/cromwell.conf`) from `$SCRATCH_DIR/Cromwell/`, the run directory
   bound read-only in the task containers (`APPTAINER_BIND`): one Slurm job per task, in
   `cromwell-executions/msi_detection/<workflow id>/`; Cromwell's log goes to `logs/<name>.log`, its
   metadata (status, outputs, every task) to `metadata/<name>.json`;
4. runs `msi_detection_post_p.py index` as a Slurm job (`srun`, in `python-MSI-3.12.sif`): extracts the
   run's `viz.tar.gz` and copies its `runs/<run_id>/` into `$RESULTS_DIR/runs/`, then rebuilds the dashboard from every
   `$RESULTS_DIR/runs/*/` (`index_msi_run.js`, `index_msi_sample.js`; the two `.html` are written once,
   only rewritten if their template changes).

The workflow itself never writes outside its tasks: the dashboard covers every run, a state that
changes with each run, so it is built after the run (`index`), not by a task (a task reading the
current reports would not be reproducible nor cacheable). `index` can be re-run at any time, e.g.
after deleting a run folder. Everything else stays in Cromwell's run folder (paths in its
metadata JSON), every task keeping its `execution/` folder (~300 MB per run).


## Tasks

Calls are named after their line on the diagram, `t<line><rank>_<step>`: `ta` BAMs (green),
`tb` msisensor-pro list (blue, panel steps included), `tc` ROI markers (red).
CPU and Memory are what each job reserves from Slurm (`runtime` of the task, passed to
`sbatch -c` / `--mem`; the defaults of `backends.conf/cromwell.conf` are never used); only `ta02_sensorpro_pro` takes `params.threads`, every other task keeps its default
`threads = 1`. No time limit is set: the partition's applies. Each task names its image with the
`apptainer` runtime attribute, not `container`: Cromwell 92 fails to hand `container`/`docker` to a
Slurm backend under WDL 1.1 ([issue #7846](https://github.com/broadinstitute/cromwell/issues/7846),
fixed after 92: back to `container` then).

| Call (task) | Runs | What it does | Main output | CPU | Memory | Set in |
|---|---|---|---|---|---|---|
| `ta01_find_bams` (`find_bams`) | python image, bash | Lists every `*.bam` of the run directory (`inputs.run`) (sub-folders with `params.recursive`), with its index (`.bam.bai` or `.bai`), its library (folder relative to `input_dir`, `.` if none) and its sample name (BAM basename). `inputs.run` is a `String` (WDL 1.1 has no `Directory`; the launchers bind it in the containers); paths are returned relative to it. With `inputs.samples`, keeps only the BAMs whose name contains one of these IDs (every library). Fails if a BAM has no index, no BAM is found, or a sample ID matches no BAM. | `bams`: one row per BAM | 1 | 512 MiB | `msi_detection-1.0.0.wdl` |
| `tb01_sensorpro_scan` (`sensorpro.scan`) | msisensor-pro | `msisensor-pro scan` of the genome: every microsatellite site. Skipped when `references.pro_list` is given. Not redone from one run to the next: Cromwell's call cache reuses the previous scan (same genome and image), as long as that run's folder is kept. | `<genome>.msisensor_pro.list` | 1 | 4 GiB | `module/msisensorpro-1.3.0.wdl` |
| `tb02_list_to_bed` (`list_to_bed`) | python image | Converts the list to BED6 (start = location − 1, end = start + unit length × repeat times). | `<list>.bed` | 1 | 2 GiB | `msi_detection-1.0.0.wdl` |
| `tb03_panel_chrom` (`align_chrom_naming`) | python image, awk | Only with `bed_panel`. Copy of the panel BED with its chromosome names (`chr1` / `1`) aligned on the list's: bedtools matches them literally. | `<panel>.chrnorm.bed` | 1 | 512 MiB | `msi_detection-1.0.0.wdl` |
| `tb04_panel_intersect` (`bedtools.intersect`) | bedtools | `bedtools intersect -wa -wb -f 1.0`: the sites fully inside a panel region, one row per overlap. | `<list>.panel.intersect.tsv` | 1 | 2 GiB | `module/bedtools-2.31.1.wdl` |
| `tb05_panel_correspondence` (`name_correspondence`) | python image, awk | Turns the intersection into a BED6 named `<site>#<probe>#<chrom:start-end>`. | `<list>.panel.bed` | 1 | 512 MiB | `msi_detection-1.0.0.wdl` |
| `tb06_mark_list` (`mark_list`) | python image | Only with `bed_panel`. Full list, same rows, `filter` = `PASS` for the panel sites, `NO` elsewhere: **the list analysed by msisensor-pro**. | `<list>.panel.list` | 1 | 2 GiB | `msi_detection-1.0.0.wdl` |
| `tc01_roi_chrom` (`align_chrom_naming`) | python image, awk | Only with `bed_roi`. Same as `tb03` with the ROI BED. | `<roi>.chrnorm.bed` | 1 | 512 MiB | `msi_detection-1.0.0.wdl` |
| `tc02_roi_intersect` (`bedtools.intersect`) | bedtools | Same as `tb04` with the ROI markers, on the panel sites when there is a panel (genome-wide otherwise). | `<list>.panel.roi.intersect.tsv` | 1 | 2 GiB | `module/bedtools-2.31.1.wdl` |
| `tc03_roi_correspondence` (`name_correspondence`) | python image, awk | Same as `tb05`; name ends with `#<marker>`. | `<list>.panel.roi.bed` | 1 | 512 MiB | `msi_detection-1.0.0.wdl` |
| `ta02_sensorpro_pro` × N (`sensorpro.pro`) | msisensor-pro | Scatter, one per BAM: `msisensor-pro pro -d <analysis list> -t <bam> -b threads -c coverage_threshold`. A failure on any BAM fails the workflow. | `<sample>`, `_all`, `_dis`, `_unstable` | `params.threads` (8: `pro -b`) | 4 GiB | CPU: params JSON; memory: `module/msisensorpro-1.3.0.wdl` |
| `ta03_report` (`report`) | python image | Reads every BAM's `_all` (`File`s, nothing copied) and builds this run's reports in `viz/runs/<run_id>/` (run page, sample pages; the all-runs dashboard is built by `index` after the run), archived with its `samples.tsv` (WDL 1.1 has no `Directory` output); log on its stderr; writes nothing outside its task. | `viz.tar.gz` | 1 | 2 GiB | `msi_detection-1.0.0.wdl` |


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
- `SCRATCH_DIR`: Cromwell's runs, with every task's files and logs, and its call cache (tree at
  the end of this section): `Cron/Cromwell/` for the cron, `<batch scratch>/Cromwell/` for the batch,
  so one call cache per mode.

`ta03_report` returns `viz.tar.gz` (this run's `runs/<run_id>/` and its `samples.tsv`), in its folder of the Cromwell run:

```text
$SCRATCH_DIR/Cromwell/cromwell-executions/msi_detection/<workflow id>/call-ta03_report/execution/viz.tar.gz
```

`index` then extracts it, copies its `runs/<run_id>/` to `$RESULTS_DIR`, adds the run's msisensor-pro `_all` files to it, and
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
`ta02_sensorpro_pro` folder of the Cromwell run (`_all` and `_unstable` listed in its metadata `outputs`).

Cromwell's folder, in `SCRATCH_DIR` (Cromwell runs from it: every path below is relative to it):

```text
$SCRATCH_DIR/Cromwell/
├── cromwell-db/                            call cache database (HSQLDB files)
├── inputs/<name>.json                      install.json + params JSON, merged (one per params JSON)
├── metadata/<name>.json                    Cromwell's metadata: status, failures, outputs, every task
├── logs/<name>.log                         Cromwell's whole log (start-up, database included)
└── cromwell-executions/
    ├── workflow-logs/workflow.<id>.log     this run's own log, without Cromwell's start-up
    └── msi_detection/
        └── <workflow id>/                  one folder per run (id in its metadata)
            ├── call-ta01_find_bams/
            ├── call-tb01_sensorpro_scan/
            ├── call-tb02_list_to_bed/ … call-tb06_mark_list/
            ├── call-tc01_roi_chrom/ … call-tc03_roi_correspondence/
            ├── call-ta02_sensorpro_pro/
            │   └── shard-0/ … shard-<N>/   one per BAM (scatter); <sample>{,_all,_dis,_unstable} in execution/
            └── call-ta03_report/
                ├── inputs/                 its input files, hard-linked (soft-linked across filesystems)
                └── execution/
                    ├── script              script run
                    ├── script.submit       sbatch command (Slurm job msi-cromwell_<id>_<call>)
                    ├── stdout, stderr      task output (report log: stderr)
                    ├── rc                  return code
                    └── viz/, samples.tsv, viz.tar.gz
```

A call taken from the cache has a `cacheCopy/execution/` folder instead, with only its declared
outputs (hard links to the run that computed them) and `rc`, `stdout`, `stderr`. The call cache
keeps a call's results as long as the run folder that computed it is kept; it keys on the input
files (size and date), the command and the runtime attributes, so the image *name*, not the
`.sif` content: a rebuilt `.sif` needs a new tag.

To run one params JSON by hand, with a summary at the end (status, failures, each task, cached or
not, outputs), use [`BASH.msi_detection.cromwell_test.sh`](companion-scripts/BASH.msi_detection.cromwell_test.sh)
(see [Companion scripts](#companion-scripts)).

## Companion scripts

Helper scripts around the workflow, in `companion-scripts/`. The `BASH.*.sh` launchers run their
python script as a Slurm job (`srun`) in `python-MSI-3.12.sif`;

| Script | Mode | Inputs | What it does |
|---|---|---|---|
| [`make_msi_locker.py`](companion-scripts/make_msi_locker.py) | cron, before the first launch | `--locker-adivar` ADIVAR's locker JSON, `--locker-msi` MSI locker (text, created if missing), `--dry-run` | Adds every run already finished by ADIVAR to the MSI locker |
| [`BASH.make_configs.sh`](companion-scripts/BASH.make_configs.sh) → [`make_configs.py`](companion-scripts/make_configs.py) | batch, before | `BAMS_DIR` group folder of BAM folders, `RUNS_FILE` runs file (run, or run + sample ID per line), `OUTPUT_DIR` | Writes one `<run>_params.json` per run of a runs file (optionally restricted to some samples); runs without BAM go to `failed-runs.txt` |
| [`BASH.msi_detection.batch.sh`](companion-scripts/BASH.msi_detection.batch.sh) | batch | `CONFIGS_DIR` folder of params JSON, `RESULTS_DIR` final reports folder (Cromwell's runs: `$SCRATCH_DIR/Cromwell/`) | Runs the workflow on every params JSON of a folder, then rebuilds the dashboard |
| [`BASH.msi_detection.cromwell_test.sh`](companion-scripts/BASH.msi_detection.cromwell_test.sh) | test | `PARAMS_JSON` (default: Run141_Juno, fake data), `WORK_DIR` Cromwell's folder (default: `Scratch/MSI/Cromwell-jp`) | `womtool validate`, then the workflow on one params JSON (no locker, no reports folder), then a summary from its metadata; rerun in the same `WORK_DIR` to reuse its call cache |
| [`BASH.export_tsv.sh`](companion-scripts/BASH.export_tsv.sh) → [`export_tsv.py`](companion-scripts/export_tsv.py) | batch, after | `RESULTS_DIR` reports folder (read only), `OUTPUT_DIR` where the TSV goes | Exports a reports folder as a TSV, one line per sample (`results_msi_<date>.tsv`) |

```sh
./companion-scripts/make_msi_locker.py --locker-adivar ADIVAR_JSON --locker-msi MSI_TXT [--dry-run]
./companion-scripts/BASH.make_configs.sh BAMS_DIR RUNS_FILE OUTPUT_DIR
./companion-scripts/BASH.msi_detection.batch.sh CONFIGS_DIR RESULTS_DIR  # in screen / tmux
./companion-scripts/BASH.msi_detection.cromwell_test.sh [PARAMS_JSON] [WORK_DIR]
./companion-scripts/BASH.export_tsv.sh RESULTS_DIR OUTPUT_DIR
```

**Note:** `BASH.msi_detection.batch.sh` and the cron use the same
[`backends.conf/cromwell.conf`](backends.conf/cromwell.conf), but each runs Cromwell from its own
`$SCRATCH_DIR/Cromwell/`, where its call cache database (`cromwell-db/`) is: the batch and cron
call caches stay apart (a batch never reuses a cron result, nor the reverse).

