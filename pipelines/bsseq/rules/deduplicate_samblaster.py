# PiGx BSseq Pipeline.
#
# Copyright 2019 Alexander Blume <alexander.blume@mdc-berlin.de>
#
# This file is part of the PiGx BSseq Pipeline.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ==========================================================================================
# Deduplication by samblaster 
#
# https://github.com/GregoryFaust/samblaster
#
rule samblaster_markdup_sort:
    input:
        DIR_mapped+"{sample}.bwameth.bam"
    output:
        bam = DIR_sorted+"{sample}.bwameth.sorted.markdup.bam",
        index = DIR_sorted+"{sample}.bwameth.sorted.markdup.bam.bai"
    params:
        threads=config['execution']['rules']['samblaster_markdup_sort']['threads'],
        memory=config['execution']['rules']['samblaster_markdup_sort']['memory']
    log:
        DIR_sorted+"{sample}_markdups.log"
    message: fmt("Deduplicating reads with samblaster for sample {wildcards.sample}")
    shell:
        nice("samtools", 
        ["view -h {input}"," | ", 
        tool("samblaster"),toolArgs("samblaster"),"2>> {log}","|",
        tool("samtools"),"sort",
         "-o {output.bam}", "-@ {params.threads}", 
         "-m {params.memory}", "-l 9","2>> {log}",";",
         tool("samtools"),"index {output.bam}"],("{log}"))
