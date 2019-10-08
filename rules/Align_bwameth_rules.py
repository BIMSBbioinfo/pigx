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
        GENOMEPATH+'bwameth_genome_preparation.output_'+ASSEMBLY+'.log'
    message: "Converting {ASSEMBLY} Genome into Bisulfite analogue with bwa-meth"
    shell:
        nice("bwameth", ["index {input}"],"{log}")

# # ==========================================================================================
# # Align whole data sets
#    


def bwameth_input(sample):
    files = list_files_TG(samples(sample,'files'),sample, '')
    return(files)

rule bwameth_align_trimmed:
    input:
        index = rules.bwameth_genome_preparation.output,
        files = lambda wc: bwameth_input(wc.sample) 
    output:
        bam = DIR_mapped+"{sample}.bwameth.bam"
    params:
      # bwa-meth parameters
        threads = config['execution']['rules']['bwameth_align_pe_raw']['threads']
    log:
        DIR_mapped+"{sample}_bwameth_mapping.log"
    message: "Mapping reads to genome using bwa-meth."
    shell:
      nice("bwameth",["--reference {GENOMEFILE}","-t {params.threads}",
          "{input.files}","2> {log}","|",
          tool('samtools'), "view -bS - > {output.bam}"])
