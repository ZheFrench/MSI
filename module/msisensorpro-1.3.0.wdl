version 1.1

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

# msisensor-pro (https://github.com/xjtu-omics/msisensor-pro) wrapped as WDL
# tasks. Each task names its image (runtime.apptainer) and calls
# `msisensor-pro` directly.
#
#   import "module/msisensorpro-1.3.0.wdl" as sensorpro
#   call sensorpro.scan { ... }   call sensorpro.pro { ... }

task scan {
	meta {
		author: "Jean-Philippe Villemin"
		email: "jeanphilipe.villemin(at)chu-montpellier.fr"
		date: "2026-10-05"
		description: "msisensor-pro scan: every microsatellite site of a reference genome."
		outputs: {
			msi_list: "Microsatellite list (msisensor-pro -d file); a link to the cached list when reused (no copy).",
		}
	}

	input {
		File genome
		String output_name = basename(genome) + ".msisensor_pro.list"

		String container
		Int threads = 1
	}

	command <<<
		set -euo pipefail

		# scan's own stdout is very verbose (per chromosome/window): discarded.
		msisensor-pro scan -d "~{genome}" -o "~{output_name}" > /dev/null
	>>>

	output {
		File msi_list = output_name
	}

	parameter_meta {
		threads: "CPUs reserved for the task."
		genome: "Reference fasta."
		output_name: "Output list name [default: <genome basename>.msisensor_pro.list]."
		container: "msisensor-pro image."
	}

	runtime {
		# Not container: (Cromwell 92 bug on Slurm, issue #7846; back to container: with Cromwell 93)
		apptainer: container
		cpu: threads
		memory: "4 GiB"
	}
}

task pro {
	meta {
		author: "Jean-Philippe Villemin"
		email: "jeanphilipe.villemin(at)chu-montpellier.fr"
		date: "2026-10-05"
		description: "msisensor-pro pro: MSI scoring of a tumor-only BAM. A msisensor-pro failure fails the task."
		outputs: {
			summary: "<prefix>: summary file (empty in tumor-only mode).",
			all: "<prefix>_all: every evaluated site with pro_p, threshold and coverage.",
			dis: "<prefix>_dis: read-length distribution per site.",
			unstable: "<prefix>_unstable: sites called unstable (pro_p > threshold).",
		}
	}

	input {
		File bam
		File bai
		File msi_list
		String prefix = basename(bam, ".bam")
		Int threads = 1
		Int? coverage_threshold

		String container
	}

	command <<<
		set -euo pipefail

		# msisensor-pro reads the index next to the BAM: bai is an input only so
		# the engine localizes it there (~{basename(bai)}).
		# pro's own stdout is very verbose (per chromosome/window): discarded;
		# its stderr stays in this task's stderr.
		msisensor-pro pro -d "~{msi_list}" -t "~{bam}" -o "~{prefix}" -b ~{threads} \
			~{if defined(coverage_threshold) then "-c " + select_first([coverage_threshold]) else ""} \
			> /dev/null
	>>>

	output {
		File summary = prefix
		File all = prefix + "_all"
		File dis = prefix + "_dis"
		File unstable = prefix + "_unstable"
	}

	parameter_meta {
		bam: "Tumor BAM."
		bai: "Its index, only declared so the engine keeps it next to the BAM."
		msi_list: "Microsatellite list to analyse (msisensor-pro -d; sites with filter NO are skipped)."
		prefix: "Output prefix [default: BAM basename]."
		threads: "msisensor-pro pro -b."
		coverage_threshold: "msisensor-pro pro -c [default: msisensor-pro's own, 15 WGS / 20 WXS]."
		container: "msisensor-pro image."
	}

	runtime {
		# Not container: (Cromwell 92 bug on Slurm, issue #7846; back to container: with Cromwell 93)
		apptainer: container
		cpu: threads
		memory: "4 GiB"
	}
}
