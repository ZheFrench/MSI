version 1.3

# MSI-Detection 0.1 - WDL pipelines.
# Copyright (C) 2026 MoBiDiC
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# WDL port of run_msisensor2.0.1.py. Every external tool is wrapped in a
# module (module/<tool>.wdl, imported below), whose tasks only name their
# Apptainer image (requirements.container) and call the tool directly: the
# backend runs it (sprocket slurm_apptainer, on this workstation's single-node
# Slurm in development, on the cluster in production).
# The pure-Python steps are subcommands of two scripts, run in the generic
# python3 image (helper_container): msi_main_helper (msi_detection_main_p.py,
# list/BED conversions) and msi_post_helper (msi_detection_post_p.py, HTML
# reports).
#
# Calls are named after their line on the diagram, t<line><rank>_<step>
# (a: BAMs, b: msisensor-pro list, c: ROI markers; WDL names can't start
# with a digit), so sprocket's calls/ folders group the same way:
#   ta01_find_bams ------------------------------------------------------------.
#   tb01_sensorpro_scan (or pro_list) -> tb02_list_to_bed                        |
#     -> [panel] tb03_panel_chrom -> tb04_panel_intersect -> tb05_panel_correspondence -> tb06_mark_list
#     -> [roi]   tc01_roi_chrom -> tc02_roi_intersect -> tc03_roi_correspondence
#   scatter: ta02_sensorpro_pro (one per BAM) <--------------------------------'
#   ta03_report: viz/, data/, log/ as Directory outputs (in the run folder)

import "module/msisensorpro-1.3.0.wdl" as msisensorpro
import "module/bedtools-2.31.1.wdl" as bedtools

# What to analyse, given in the params JSON as one object:
#   "msi_detection.inputs": {"run": "<run directory>"}
# (a Directory: the engine mounts it in the tasks that use it)
struct InputDirs {
	Directory run
}

# Reference files, given in the params JSON as one object:
#   "msi_detection.references": {"genome": "...", "bed_panel": "...", "bed_roi": "..."}
# (bed_panel, bed_roi, pro_list optional: leave them out to skip that step)
struct References {
	File genome
	File? pro_list
	File? bed_panel
	File? bed_roi
}

# Run settings, given in the params JSON as one object:
#   "msi_detection.params": {"threads": 8, "coverage_threshold": 1, "recursive": true,
#                            "p_process": {"min_reads": 200, "pval_cutoff": 0.1, "sample_na_fraction": 0.5}}
# p_process: report settings (ta03_report), msisensor-pro itself doesn't use them.
struct PProcess {
	Int min_reads
	Float pval_cutoff
	Float sample_na_fraction
}

struct RunParams {
	Int threads
	Int coverage_threshold
	Boolean recursive
	PProcess p_process
}

# Where results go is not a workflow input: every output stays in sprocket's
# run folder (the workflow writes nothing outside its tasks). BASH.sprocket.sh
# sets the destinations: -o for sprocket's runs, and the reports folder that
# `msi_detection_post_p.py index` fills after each run (adds runs/<run_id>/,
# rebuilds the dashboard of every run). The scan is not redone from one run
# to the next: sprocket's call cache reuses the previous result (same genome,
# image and threads).

workflow msi_detection {

	meta {
		author: "Jean-Philippe Villemin"
		email: "jeanphilipe.villemin(at)chu-montpellier.fr"
		version: "0.1.0"
		date: "2026-10-05"
		description: "MSI detection with msisensor-pro (tumor-only) on every BAM of a run directory, restricted to a sequencing panel and/or classic MSI markers, with HTML reports per run and per sample."
		outputs: {
			msi_reference_list: "Genome-wide msisensor-pro scan list (<genome>.msisensor_pro.list).",
			msi_analysis_list: "List actually analysed by msisensor-pro pro: the panel-marked list (filter PASS/NO) with a panel, the genome-wide one otherwise.",
			panel_correspondence_bed: "Microsatellite sites inside the panel regions, named <site>#<probe>#<chrom:start-end> (only with bed_panel).",
			roi_correspondence_bed: "Microsatellite sites inside the ROI markers, named ...#<marker> (only with bed_roi).",
			msisensorpro_all: "Per BAM: msisensor-pro <sample>_all (every evaluated site); empty for a failed BAM.",
			msisensorpro_unstable: "Per BAM: msisensor-pro <sample>_unstable (unstable sites only); empty for a failed BAM.",
			report_viz: "This run's HTML/JS reports, runs/<run_id>/ (added to the reports folder by msi_detection_post_p.py index, which also rebuilds the dashboard).",
			report_data: "This run's raw msisensor-pro outputs, <run_id>/<library>/<sample>{,_all,_dis,_unstable}.",
			report_log: "This run's report log, <run_id>_MSI_<timestamp>.log.",
			all_bams_ok: "false when msisensor-pro failed on some BAM (reported as failed in the reports).",
		}
	}

	input {
		InputDirs inputs

		References references

		# Tool images and python helper: set here, not in the params JSON
		# (still overridable there: production images in the prod params JSON).
		String msisensorpro_container = "file:///media/jp/SSD2/CHU/code/pipeline-jp/wdl/apptainer-recipes/msisensor-pro-1.3.0.sif"
		String bedtools_container = "file:///media/jp/SSD2/CHU/code/pipeline-jp/wdl/apptainer-recipes/bedtools-2.31.1.sif"
		String helper_container = "file:///media/jp/SSD2/CHU/code/pipeline-jp/wdl/apptainer-recipes/python-3.12.sif"
		# Controls, not analysed: BAMs whose name starts with one of these.
		Array[String] exclude_prefixes = ["TN", "TPH"]

		File msi_main_helper = "/media/jp/SSD2/CHU/code/pipeline-jp/wdl/scripts/msi_detection_main_p.py"
		File msi_post_helper = "/media/jp/SSD2/CHU/code/pipeline-jp/wdl/scripts/msi_detection_post_p.py"

		RunParams params
	}

	# Named after the input directory (e.g. "Run141_Juno"), so re-running it
	# updates that same run in place in the reports.
	String run_id = basename(sub("~{inputs.run}", "/+$", ""))

	call find_bams as ta01_find_bams {
		input:
			threads = params.threads,
			container = helper_container,
			input_dir = inputs.run,
			recursive = params.recursive,
			exclude_prefixes = exclude_prefixes,
	}

	# --- Reference list: genome-wide scan, then optional panel/ROI restriction ---

	if (!defined(references.pro_list)) {
		call msisensorpro.scan as tb01_sensorpro_scan {
			input:
				threads = params.threads,
				genome = references.genome,
				container = msisensorpro_container,
		}
	}
	File msi_list = select_first([references.pro_list, tb01_sensorpro_scan.msi_list])
	String list_stem = basename(msi_list, ".list")

	call list_to_bed as tb02_list_to_bed {
		input:
			threads = params.threads,
			container = helper_container,
			msi_list = msi_list,
			output_name = list_stem + ".bed",
			msi_helper = msi_main_helper,
	}

	if (defined(references.bed_panel)) {
		call align_chrom_naming as tb03_panel_chrom {
			input:
				threads = params.threads,
				container = helper_container,
				region_bed = select_first([references.bed_panel]),
				reference_bed = tb02_list_to_bed.bed,
		}

		# Sites fully inside a panel region (-f 1.0), one row per overlap.
		call bedtools.intersect as tb04_panel_intersect {
			input:
				threads = params.threads,
				a = tb02_list_to_bed.bed,
				b = tb03_panel_chrom.bed,
				output_name = list_stem + ".panel.intersect.tsv",
				wa = true,
				wb = true,
				min_overlap_a = 1.0,
				container = bedtools_container,
		}

		call name_correspondence as tb05_panel_correspondence {
			input:
				threads = params.threads,
				container = helper_container,
				intersection = tb04_panel_intersect.intersection,
				output_name = list_stem + ".panel.bed",
				append_position = true,
		}

		call mark_list as tb06_mark_list {
			input:
				threads = params.threads,
				container = helper_container,
				msi_list = msi_list,
				panel_correspondence_bed = tb05_panel_correspondence.bed,
				output_name = list_stem + ".panel.list",
				msi_helper = msi_main_helper,
		}
	}
	# With a panel, msisensor-pro pro analyses the panel-marked list (PASS/NO)
	# instead of the full genome-wide one.
	File analysis_list = select_first([tb06_mark_list.panel_list, msi_list])

	if (defined(references.bed_roi)) {
		# The ROI markers are expected to sit inside the panel's own amplicons:
		# cross-reference them against the panel sites when there is a panel,
		# against the whole genome otherwise.
		File roi_source_bed = select_first([tb05_panel_correspondence.bed, tb02_list_to_bed.bed])
		String roi_stem = basename(roi_source_bed, ".bed") + ".roi"

		call align_chrom_naming as tc01_roi_chrom {
			input:
				threads = params.threads,
				container = helper_container,
				region_bed = select_first([references.bed_roi]),
				reference_bed = roi_source_bed,
		}

		call bedtools.intersect as tc02_roi_intersect {
			input:
				threads = params.threads,
				a = roi_source_bed,
				b = tc01_roi_chrom.bed,
				output_name = roi_stem + ".intersect.tsv",
				wa = true,
				wb = true,
				min_overlap_a = 1.0,
				container = bedtools_container,
		}

		call name_correspondence as tc03_roi_correspondence {
			input:
				threads = params.threads,
				container = helper_container,
				intersection = tc02_roi_intersect.intersection,
				output_name = roi_stem + ".bed",
		}
	}

	# --- One msisensor-pro pro per BAM ---

	scatter (bam_row in ta01_find_bams.bams) {
		# bam_row: bam, bai (paths relative to inputs.run), library ("." if
		# none), sample_name
		String bam_path = "~{inputs.run}/~{bam_row[0]}"
		String sample_name = bam_row[3]
		String library = bam_row[2]

		call msisensorpro.pro as ta02_sensorpro_pro {
			input:
				bam = bam_path,
				bai = "~{inputs.run}/~{bam_row[1]}",
				prefix = bam_row[3],
				msi_list = analysis_list,
				threads = params.threads,
				coverage_threshold = params.coverage_threshold,
				# a failed BAM is reported (status = failed) by ta03_report
				tolerate_failure = true,
				container = msisensorpro_container,
		}
	}

	# --- Reports ---

	call report as ta03_report {
		input:
			threads = params.threads,
			container = helper_container,
			run_id = run_id,
			sample_names = sample_name,
			libraries = library,
			bam_paths = bam_path,
			statuses = ta02_sensorpro_pro.status,
			all_files = ta02_sensorpro_pro.all,
			unstable_files = ta02_sensorpro_pro.unstable,
			dis_files = ta02_sensorpro_pro.dis,
			summary_files = ta02_sensorpro_pro.summary,
			roi_correspondence_bed = tc03_roi_correspondence.bed,
			genome = "~{references.genome}",
			min_reads = params.p_process.min_reads,
			pval_cutoff = params.p_process.pval_cutoff,
			sample_na_fraction = params.p_process.sample_na_fraction,
			msi_helper = msi_post_helper,
	}

	output {
		File msi_reference_list = msi_list
		File msi_analysis_list = analysis_list
		File? panel_correspondence_bed = tb05_panel_correspondence.bed
		File? roi_correspondence_bed = tc03_roi_correspondence.bed
		Array[File] msisensorpro_all = ta02_sensorpro_pro.all
		Array[File] msisensorpro_unstable = ta02_sensorpro_pro.unstable
		Directory report_viz = ta03_report.viz
		Directory report_data = ta03_report.data
		Directory report_log = ta03_report.log
		Boolean all_bams_ok = ta03_report.all_bams_ok
	}

	parameter_meta {
		inputs: "run: run directory; every *.bam in it is processed (one sub-folder per library, e.g. alignments_libA/, with params.recursive = true). Its name is the run_id."
		references: "genome: reference fasta the BAMs are aligned on, scanned by msisensor-pro scan; pro_list (optional): already generated scan list of genome, skips the scan; bed_panel (optional): panel BED (plain BED, no track line, e.g. *.amplicons.bed), analysis restricted to its microsatellites; bed_roi (optional): classic MSI markers (BAT26, NR21, ...), highlighted in the reports."
		params: "threads: msisensor-pro pro -b (threads per BAM); coverage_threshold: msisensor-pro pro -c (msisensor-pro's own default: 15 WGS, 20 WXS); recursive: also look for BAMs in the run directory's sub-folders (one per library); p_process (reports only): min_reads, a site counts (% unstable, site status) only with >= this many reads, pval_cutoff, a site is unstable when its p-value is above this, and sample_na_fraction, sample status NA when more than this fraction of its loci are NA or MSU."
		msisensorpro_container: "msisensor-pro image [default: apptainer-recipes/msisensor-pro-1.3.0.sif]; file:///...sif on slurm_apptainer."
		bedtools_container: "bedtools image [default: apptainer-recipes/bedtools-2.31.1.sif]; file:///...sif on slurm_apptainer."
		helper_container: "Generic python3 + awk image of the python/bash tasks [default: apptainer-recipes/python-3.12.sif]; the script is not in it."
		exclude_prefixes: "BAMs whose name starts with one of these are not analysed (controls) [default: TN, TPH]."
		msi_main_helper: "msi_detection_main_p.py, given to list_to_bed and mark_list as a File [default: scripts/msi_detection_main_p.py]."
		msi_post_helper: "msi_detection_post_p.py, given to the report task as a File [default: scripts/msi_detection_post_p.py]."
	}
}

task find_bams {
	meta {
		description: "List every *.bam of input_dir (except those whose name starts with one of exclude_prefixes) with its index, library folder and sample name, as paths relative to input_dir (the caller rebuilds them from its own input_dir value)."
		outputs: {
			bams: "One row per BAM: bam path, index path (both relative to input_dir), library (folder relative to input_dir, \".\" if none), sample name (BAM basename).",
		}
	}

	input {
		Directory input_dir
		Boolean recursive
		Array[String] exclude_prefixes = []
		Int threads = 1

		String container
	}

	command <<<
		set -euo pipefail
		dir="~{input_dir}"
		dir="${dir%/}"

		# Paths relative to input_dir: inside a container, input_dir is mounted
		# at its own path, meaningless outside the task.
		(cd "$dir" && LC_ALL=C find -L . ~{if recursive then "" else "-maxdepth 1"} -name '*.bam' -type f) \
			| sed 's|^\./||' | LC_ALL=C sort > all_bams.txt
		# Controls (BAM name starting with one of exclude_prefixes) left out.
		exclude="~{sep("|", exclude_prefixes)}"
		if [ -n "$exclude" ]; then
			{ grep -Ev "(^|/)($exclude)[^/]*\.bam$" all_bams.txt || true; } > bam_list.txt
			{ grep -E "(^|/)($exclude)[^/]*\.bam$" all_bams.txt || true; } | sed 's/^/Excluded (control): /' >&2
		else
			cp all_bams.txt bam_list.txt
		fi
		while read -r bam; do
			bai="$bam.bai"
			[ -e "$dir/$bai" ] || bai="${bam%.bam}.bai"
			[ -e "$dir/$bai" ] || { echo "No index (.bam.bai or .bai) for $bam" >&2; exit 1; }
			printf '%s\t%s\t%s\t%s\n' "$bam" "$bai" "$(dirname "$bam")" "$(basename "$bam" .bam)"
		done < bam_list.txt > bams.tsv

		[ -s bams.tsv ] || { echo "No .bam file found in $dir (recursive: ~{recursive})" >&2; exit 1; }
	>>>

	output {
		Array[Array[String]] bams = read_tsv("bams.tsv")
	}

	parameter_meta {
		container: "python3 + awk image (helper_container)."
		threads: "CPUs reserved for the task."
		input_dir: "Run directory to scan for *.bam."
		recursive: "Descend into sub-folders (one per library)."
		exclude_prefixes: "BAMs whose name starts with one of these are left out (controls)."
	}

	requirements {
		container: container
		cpu: threads
		memory: "512 MiB"
	}
}

task list_to_bed {
	meta {
		description: "msisensor-pro .list -> BED6 (one interval per microsatellite site)."
		outputs: {
			bed: "BED6 of the sites: start = location - 1, end = start + unit length x repeat times, name = <times>-<unit>-<left flank>-<right flank>.",
		}
	}

	input {
		File msi_list
		String output_name
		File msi_helper
		Int threads = 1

		String container
	}

	command <<<
		set -euo pipefail
		python3 "~{msi_helper}" list-to-bed "~{msi_list}" "~{output_name}"
	>>>

	output {
		File bed = output_name
	}

	parameter_meta {
		container: "python3 + awk image (helper_container)."
		threads: "CPUs reserved for the task."
		msi_list: "msisensor-pro scan list."
		output_name: "Output BED file name."
		msi_helper: "Path to msi_detection_main_p.py."
	}

	requirements {
		container: container
		cpu: threads
		memory: "2 GiB"
	}
}

task align_chrom_naming {
	meta {
		description: "Copy of region_bed with its chromosome names aligned on reference_bed's convention (\"chr1\" UCSC vs \"1\" Ensembl, from its first line): bedtools matches names literally and silently finds no overlap on a mismatch."
		outputs: {
			bed: "region_bed with chromosome names aligned (empty lines dropped).",
		}
	}

	input {
		File region_bed
		File reference_bed
		String output_name = basename(region_bed, ".bed") + ".chrnorm.bed"
		Int threads = 1

		String container
	}

	command <<<
		set -euo pipefail
		if awk -F'\t' 'NF { found = 1; exit ($1 ~ /^chr/) ? 0 : 1 } END { if (!found) exit 0 }' "~{reference_bed}"; then
			use_chr=1
		else
			use_chr=0
		fi
		awk -F'\t' -v OFS='\t' -v use_chr="$use_chr" '
			!NF { next }
			use_chr && $1 !~ /^chr/ { $1 = "chr" $1 }
			!use_chr && $1 ~ /^chr/ { $1 = substr($1, 4) }
			{ print }
		' "~{region_bed}" > "~{output_name}"
	>>>

	output {
		File bed = output_name
	}

	parameter_meta {
		container: "python3 + awk image (helper_container)."
		threads: "CPUs reserved for the task."
		region_bed: "Panel or ROI BED."
		reference_bed: "Microsatellite BED whose naming convention is followed."
		output_name: "Output BED name [default: <region_bed basename>.chrnorm.bed]."
	}

	requirements {
		container: container
		cpu: threads
		memory: "512 MiB"
	}
}

task name_correspondence {
	meta {
		description: "Turn a `bedtools intersect -wa -wb` of the microsatellite BED6 (A) with a panel/ROI BED (B) into a correspondence BED6: one row per overlap, name <site>#<region>, score 0, strand +."
		outputs: {
			bed: "Correspondence BED6, name <site>#<region>[#<chrom>:<start>-<end>].",
		}
	}

	input {
		File intersection
		String output_name
		# Append the matched region's own "#<chrom>:<start>-<end>" (chrom without
		# "chr") after its name: set for the panel step, so the probe position
		# travels with its name; the chained ROI step then appends the marker
		# name last: "<site>#<probe>#<position>#<marker>".
		Boolean append_position = false
		Int threads = 1

		String container
	}

	command <<<
		set -euo pipefail
		# $7-$9 are the matched region's own chrom/start/end, $10 its name.
		if [ "~{append_position}" = true ]; then
			awk 'BEGIN{OFS="\t"} {chrom=$7; sub(/^chr/,"",chrom); print $1,$2,$3,$4"#"$10"#"chrom":"$8"-"$9,0,"+"}' \
				"~{intersection}" > "~{output_name}"
		else
			awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4"#"$10,0,"+"}' "~{intersection}" > "~{output_name}"
		fi
		echo "$(wc -l < "~{output_name}") site(s) in ~{output_name}" >&2
	>>>

	output {
		File bed = output_name
	}

	parameter_meta {
		container: "python3 + awk image (helper_container)."
		threads: "CPUs reserved for the task."
		intersection: "bedtools intersect -wa -wb output (A: microsatellite BED6, B: panel/ROI BED)."
		output_name: "Output BED name."
		append_position: "Append the region's #<chrom>:<start>-<end> to the name (panel step)."
	}

	requirements {
		container: container
		cpu: threads
		memory: "512 MiB"
	}
}

task mark_list {
	meta {
		description: "Full msisensor-pro .list with filter = PASS for the panel sites, NO otherwise."
		outputs: {
			panel_list: "Full scan list, same rows and order, filter column = PASS for the panel sites and NO elsewhere: the list msisensor-pro pro analyses.",
		}
	}

	input {
		File msi_list
		File panel_correspondence_bed
		String output_name
		File msi_helper
		Int threads = 1

		String container
	}

	command <<<
		set -euo pipefail
		python3 "~{msi_helper}" mark-list "~{msi_list}" "~{panel_correspondence_bed}" "~{output_name}"
	>>>

	output {
		File panel_list = output_name
	}

	parameter_meta {
		container: "python3 + awk image (helper_container)."
		threads: "CPUs reserved for the task."
		msi_list: "Genome-wide scan list."
		panel_correspondence_bed: "Panel correspondence BED (panel_intersect)."
		output_name: "Output list file name."
		msi_helper: "Path to msi_detection_main_p.py."
	}

	requirements {
		container: container
		cpu: threads
		memory: "2 GiB"
	}
}

task report {
	meta {
		description: "Gather every BAM's raw msisensor-pro outputs into data/<run_id>/<library>/ and build this run's HTML reports into viz/runs/<run_id>/ (run and sample pages; the all-runs dashboard is built by msi_detection_post_p.py index), all returned as Directory outputs."
		outputs: {
			viz: "This run's HTML/JS reports: runs/<run_id>/.",
			data: "This run's raw outputs: <run_id>/<library>/<sample>{,_all,_dis,_unstable}.",
			log: "This run's log: <run_id>_MSI_<timestamp>.log.",
			all_bams_ok: "false when some BAM failed (reported as failed in the reports, which are still all written).",
		}
	}

	input {
		String run_id
		# One entry per BAM, in the same order in every array
		Array[String] sample_names
		Array[String] libraries
		Array[String] bam_paths
		Array[String] statuses
		Array[File] all_files
		Array[File] unstable_files
		Array[File] dis_files
		Array[File] summary_files
		File? roi_correspondence_bed
		String genome
		Int min_reads
		Float pval_cutoff
		Float sample_na_fraction
		File msi_helper
		Int threads = 1

		String container
	}

	command <<<
		set -euo pipefail
		mkdir -p viz data log

		# One line per BAM: sample, library, bam, status, _all, _unstable, _dis, summary
		paste ~{write_lines(sample_names)} ~{write_lines(libraries)} ~{write_lines(bam_paths)} \
			~{write_lines(statuses)} ~{write_lines(all_files)} ~{write_lines(unstable_files)} \
			~{write_lines(dis_files)} ~{write_lines(summary_files)} > samples.tsv

		# Exit 3 = every report written but some BAM failed: kept as
		# all_bams_ok = false instead of failing the task (and losing its outputs).
		rc=0
		python3 "~{msi_helper}" report \
			--run-id "~{run_id}" \
			--samples samples.tsv \
			--genome "~{genome}" \
			--data-output data \
			--viz-output viz \
			--log-output log \
			--min-reads ~{min_reads} \
			--pval-cutoff ~{pval_cutoff} \
			--sample-na-fraction ~{sample_na_fraction} \
			~{if defined(roi_correspondence_bed) then "--roi-bed \"" + select_first([roi_correspondence_bed]) + "\"" else ""} \
			|| rc=$?
		[ "$rc" = 0 ] || [ "$rc" = 3 ] || exit "$rc"
		echo "$rc" > report_rc
	>>>

	output {
		Directory viz = "viz"
		Directory data = "data"
		Directory log = "log"
		Boolean all_bams_ok = read_int("report_rc") == 0
	}

	parameter_meta {
		container: "python3 + awk image (helper_container)."
		threads: "CPUs reserved for the task."
		run_id: "Run name (the run directory's name)."
		sample_names: "Per BAM: sample name (BAM basename)."
		libraries: "Per BAM: library folder (\".\" if none)."
		bam_paths: "Per BAM: BAM path (only recorded in the reports)."
		statuses: "Per BAM: msisensor-pro status, ok or failed."
		all_files: "Per BAM: msisensor-pro _all (empty if failed)."
		unstable_files: "Per BAM: msisensor-pro _unstable (empty if failed)."
		dis_files: "Per BAM: msisensor-pro _dis (empty if failed)."
		summary_files: "Per BAM: msisensor-pro summary (empty if failed)."
		roi_correspondence_bed: "ROI correspondence BED, to flag the classic markers in the reports."
		genome: "Reference fasta path (only recorded in the reports)."
		min_reads: "A site counts (% unstable, site status) only with >= this many reads."
		pval_cutoff: "A site is unstable when its p-value is above this."
		sample_na_fraction: "Sample status NA when more than this fraction of its loci are NA or MSU."
		msi_helper: "msi_detection_post_p.py."
	}

	requirements {
		container: container
		cpu: threads
		memory: "2 GiB"
	}
}
