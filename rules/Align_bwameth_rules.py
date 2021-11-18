# PiGx BSseq Pipeline.
#
# Copyright 2019 Alexander Blume <alexander.blume@mdc-berlin.de>
# Copyright © 2017 Katarzyna Wreczycka <katwre@gmail.com>
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

### This file was originally taken and modified from  
### https://github.com/katwre/makeWGBSnake/blob/master/Rules/Align_bwameth_rules.py

# WGBS pipeline
#
# Copyright © 2018 Katarzyna Wreczycka katarzyna.wreczycka@mdc-berlin.de
# This pipeline is heavily based on the PiGx BSseq pipeline github.com/BIMSBbioinfo/pigx_bsseq

# # # ==========================================================================================
# # # Mapping bisulfite reads with bwa-meth:
# # 

# Notes:
# bwa meth paper https://arxiv.org/pdf/1401.1129.pdf
# https://www.biostars.org/p/93726/
# bwa meth and snp http://seqanswers.com/forums/showthread.php?t=40562

# 
# # ==========================================================================================
# # Generate methyl-converted version of the reference genome:
#       

rule bwameth_genome_preparation:
    input:
        ancient(GENOMEFILE)
    output:
        GENOMEFILE+".bwameth.c2t.sa",
        GENOMEFILE+".bwameth.c2t.amb",
        GENOMEFILE+".bwameth.c2t.ann",
        GENOMEFILE+".bwameth.c2t.pac",
        GENOMEFILE+".bwameth.c2t.bwt",
        GENOMEFILE+".bwameth.c2t"
    log:
        os.path.join(GENOMEPATH,'bwameth_genome_preparation.output_'+ASSEMBLY+'.log')
    message: "Converting {ASSEMBLY} Genome into Bisulfite analogue with bwa-meth"
    shell:
        nice("bwameth", ["index {input}"],"{log}")

rule bwameth_touch_index:
    input:
        GENOMEFILE+".bwameth.c2t.sa",
        GENOMEFILE+".bwameth.c2t.amb",
        GENOMEFILE+".bwameth.c2t.ann",
        GENOMEFILE+".bwameth.c2t.pac",
        GENOMEFILE+".bwameth.c2t.bwt"
    output:
        GENOMEFILE+".bwameth.c2t_was.touched"
    message: "Update timestamp for {ASSEMBLY} Genome Index"
    shell:
        "sleep 60; touch {input};touch {output}"


# # ==========================================================================================
# # Align whole data sets
#    


def bwameth_input(sample):
    files = list_files_TG(samplesheet(sample,'files'),sample, '')
    return(files)

rule bwameth_align_trimmed:
    input:
        rules.bwameth_touch_index.output,
        index = rules.bwameth_genome_preparation.output,
        files = lambda wc: bwameth_input(wc.sample)
    output:
        bam = DIR_mapped+"{sample}.bwameth.bam"
    params:
      # bwa-meth parameters
        threads = config['execution']['rules']['bwameth_align_trimmed']['threads']
    log:
        DIR_mapped+"{sample}_bwameth_mapping.log"
    message: fmt("Mapping reads to genome using bwa-meth for sample {wildcards.sample}.")
    shell:
      nice("bwameth",["--reference {GENOMEFILE}","-t {params.threads}",
          "{input.files}","2> {log}","|",
          tool('samtools'), "view -bS - > {output.bam}"])
