# PiGx BSseq Pipeline.
#
# Copyright © 2019 Alexander Blume <alexander.blume@mdc-berlin.de>
# Copyright © 2018 Katarzyna Wreczycka katarzyna.wreczycka@mdc-berlin.de
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

# # ==========================================================================================
# # Alignment stats extracted with samtools
#  
# are included into multiqc report

rule idxstats_samtools:
  input:
    DIR_mapped+"{prefix}.bam",
  output:
    DIR_mapped+"{prefix}.idxstats.txt"
  shell:
    nice("samtools" ,["idxstats {input} > {output}"])


rule stat_samtools:
  input:
    DIR_mapped+"{prefix}.bam",
  output:
    DIR_mapped+"{prefix}.stats.txt"
  shell:
    nice("samtools", ['stats {input} > {output}'])


rule flagstat_samtools:
  input:
    DIR_mapped+"{prefix}.bam",
  output:
    DIR_mapped+"{prefix}.flagstat.txt"
  shell:
    nice("samtools" ,["flagstat {input} > {output}"])