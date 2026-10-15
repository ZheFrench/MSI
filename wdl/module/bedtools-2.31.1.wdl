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

# bedtools (https://github.com/arq5x/bedtools2) wrapped as WDL tasks. Each task
# names its image (requirements.container) and calls `bedtools` directly.
#
#   import "module/bedtools-2.31.1.wdl" as bedtools
#   call bedtools.intersect { ... }

task intersect {
	meta {
		author: "Jean-Philippe Villemin"
		email: "jeanphilipe.villemin(at)chu-montpellier.fr"
		version: "0.1.0"
		date: "2026-10-05"
		description: "bedtools intersect of a with b. Chromosome names are matched literally (chr1 vs 1 gives no overlap)."
		outputs: {
			intersection: "bedtools intersect output (stdout), named output_name.",
		}
	}

	input {
		File a
		File b
		String output_name
		Boolean wa = false
		Boolean wb = false
		Float? min_overlap_a

		String container
		Int threads = 1
	}

	command <<<
		set -euo pipefail
		bedtools intersect -a "~{a}" -b "~{b}" \
			~{if wa then "-wa" else ""} ~{if wb then "-wb" else ""} \
			~{if defined(min_overlap_a) then "-f " + select_first([min_overlap_a]) else ""} \
			> "~{output_name}"
	>>>

	output {
		File intersection = output_name
	}

	parameter_meta {
		threads: "CPUs reserved for the task."
		a: "BED/GFF/VCF/BAM file A (-a)."
		b: "BED/GFF/VCF/BAM file B (-b)."
		output_name: "Output file name."
		wa: "-wa: write the original entry of A for each overlap."
		wb: "-wb: write the original entry of B for each overlap."
		min_overlap_a: "-f: minimum overlap as a fraction of A (1.0: A fully inside B)."
		container: "bedtools image."
	}

	requirements {
		container: container
		cpu: threads
		memory: "2 GiB"
	}
}
