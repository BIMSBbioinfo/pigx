### This file was originally taken and modified from  

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
# # Filter BAM files based on MAPQ score
#  

# rule filtered_MAPQ_stats:
#   input:
#    DIR_mapped_persample_filtered+"{sample}_MAPQ_"+config['args']['MAPQ']+"_sorted.bam"
#   output:
#     DIR_mapped_persample_filtered+"{sample}_MAPQ_"+config['args']['MAPQ']+".idxstats.txt",
#     DIR_mapped_persample_filtered+"{sample}_MAPQ_"+config['args']['MAPQ']+".stats.txt",
#     DIR_mapped_persample_filtered+"{sample}_MAPQ_"+config['args']['MAPQ']+".flagstat.txt"
#   shell:
#     "{tools}/samtools idxstats {input} > {output};{tools}/samtools stats {input} > {output};{tools}/samtools flagstat {input} > {output}"


# rule filtered_MAPQ_sort_index:
#   input:
#     DIR_mapped_persample_filtered+"{sample}_MAPQ_"+config['args']['MAPQ']+".bam"
#   output:
#     DIR_mapped_persample_filtered+"{sample}_MAPQ_"+config['args']['MAPQ']+"_sorted.bam"
#   params:
#     MAPQ = config['args']['MAPQ'],
#     sort_args = config['args']['sambamba_sort'],
#     tmpdir=DIR_mapped+""
#   log:
#     DIR_mapped_persample_filtered+"{sample}_bwameth_sort.log"
#   shell:
#     "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args} > {log} 2> {log}.err"


# rule filter_MAPQ:
#   input:
#     DIR_mapped_sample+"{sample}_sorted.bam"
#   output:
#     DIR_mapped_persample_filtered+"{sample}_MAPQ_"+config['args']['MAPQ']+".bam"
#   params:
#     MAPQ = config['args']['MAPQ']
#   log:
#     DIR_mapped_persample_filtered+"{sample}_filter_MAPQ.log"
#   shell:
#     '{tools}/sambamba view --format=bam -F "mapping_quality >= {params.MAPQ}" {input} > {output} 2> {log}.err'


# ############################### todo [START]


# rule filtered_MAPQ_stats1:
#   input:
#     DIR_mapped_persample_filtered+"{sample}_not.duplicate.failed_quality_control.secondary_alignment..MAPQ_"+config['args']['MAPQ']+"_sorted.bam"
#   output:
#     DIR_mapped_persample_filtered+"{sample}_not.duplicate.failed_quality_control.secondary_alignment..MAPQ_"+config['args']['MAPQ']+".idxstats.txt",
#     DIR_mapped_persample_filtered+"{sample}_not.duplicate.failed_quality_control.secondary_alignment..MAPQ_"+config['args']['MAPQ']+".stats.txt",
#     DIR_mapped_persample_filtered+"{sample}_not.duplicate.failed_quality_control.secondary_alignment..MAPQ_"+config['args']['MAPQ']+".flagstat.txt"
#   shell:
#     "{tools}/samtools idxstats {input} > {output};{tools}/samtools stats {input} > {output};{tools}/samtools flagstat {input} > {output}"


# rule filtered_MAPQ_sort_index1:
#   input:
#     DIR_mapped_persample_filtered+"{sample}_not.duplicate.failed_quality_control.secondary_alignment..MAPQ_"+config['args']['MAPQ']+".bam"
#   output:
#     DIR_mapped_persample_filtered+"{sample}_not.duplicate.failed_quality_control.secondary_alignment..MAPQ_"+config['args']['MAPQ']+"_sorted.bam"
#   params:
#     MAPQ = config['args']['MAPQ'],
#     sort_args = config['args']['sambamba_sort'],
#     tmpdir=DIR_mapped+""
#   log:
#     DIR_mapped_persample_filtered+"{sample}_bwameth_sort.log"
#   shell:
#     "{tools}/sambamba sort -t 20 {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args} > {log} 2> {log}.err"

# rule filter_MAPQ1:
#   input:
#     DIR_mapped_sample+"{sample}_sorted.bam"
#   output:
#     DIR_mapped_persample_filtered+"{sample}_not.duplicate.failed_quality_control.secondary_alignment..MAPQ_"+config['args']['MAPQ']+".bam"
#   params:
#     MAPQ = config['args']['MAPQ']
#   log:
#     DIR_mapped_persample_filtered+"{sample}_filter_MAPQ.log"
#   shell:
#     '{tools}/sambamba view --format=bam -F "not (duplicate or failed_quality_control or secondary_alignment) and mapping_quality >= 40" {input} > {output} 2> {log}.err'


# # rule filter_bam:
# #   input:
# #     DIR_mapped_sample+"{sample}_sorted.bam"
# #   output:
# #     DIR_mapped_persample_filtered+"{sample}_MAPQ_"+config['args']['MAPQ']+".bam"
# #   params:
# #     #MAPQ = config['args']['MAPQ']
# #     pattern = config['args']['filter_bam']
# #   log:
# #     DIR_mapped_persample_filtered+"{sample}_filter_MAPQ.log"
# #   shell:
# #     'echo {tools}/sambamba view --format=bam -F "{params.pattern}" {input} > {output} 2> {log}.err'

# ############################### todo [END]




# # ==========================================================================================
# # Merge bam files from the same sample but different lanes
#  

#rule merge_lanes_idxstats:
#  input:
#    DIR_mapped_sample+"{sample}_sorted.bam"
#  output:
#    DIR_mapped_sample+"{sample}.idxstats.txt"
#  shell:
#    "{tools}/samtools idxstats {input} > {output}"


#rule merge_lanes_stat:
#  input:
#    DIR_mapped_sample+"{sample}_sorted.bam"
#  output:
#    DIR_mapped_sample+"{sample}.stats.txt"
#  shell:
#    "{tools}/samtools stats {input} > {output}"

#rule merge_lanes_flagstat:
#  input:
#    DIR_mapped_sample+"{sample}_sorted.bam"
#  output:
#    DIR_mapped_sample+"{sample}.flagstat.txt"
#  shell:
#    "{tools}/samtools flagstat {input} > {output}"
    
    
#rule merge_lanes_bwameth:
#     input:
#         #config['lanes_file'],
#         infiles=lambda sample: [DIR_mapped+x+"/"+x+".bwameth_sorted.bam" for x in SAMPLES_LANES[sample[0]] ]
#     output:
#         DIR_mapped_sample+"{sample}_sorted.bam"
#     shell:
#         "{tools}/sambamba merge -t 5 {output} {input.infiles}"


# # ==========================================================================================
# # Align. stats
#  

rule idxstats_bwameth:
  input:
    DIR_mapped+"{sample}.bwameth_sorted.bam"
  output:
    DIR_mapped+"{sample}.bwameth.idxstats.txt"
  shell:
    nice("samtools" ,["idxstats {input} > {output}"])


rule stat_bwameth:
  input:
    DIR_mapped+"{sample}.bwameth_sorted.bam"
  output:
    DIR_mapped+"{sample}.bwameth.stats.txt"
  shell:
    nice("samtools", ['stats {input} > {output}'])


rule flagstat_bwameth:
  input:
    DIR_mapped+"{sample}.bwameth_sorted.bam"
  output:
    DIR_mapped+"{sample}.bwameth.flagstat.txt"
  shell:
    nice("samtools" ,["flagstat {input} > {output}"])


rule sort_index_bam_bwameth:
  input:
    DIR_mapped+"{sample}.bwameth.bam"
  output:
    DIR_mapped+"{sample}.bwameth_sorted.bam"
  params:
    tmpdir=DIR_mapped+""
  log:
    DIR_mapped+"{sample}_bwameth_sort.log"
  shell:
    nice('sambamba_sort', ["sort {input} --tmpdir={params.tmpdir} -o {output}"],"{log}")


# # ==========================================================================================
# # Align whole data sets
#    

if not SUBSET_READS and not NOTRIMMING:
  rule bwameth_align_pe_trimmed:
     input:
         index = rules.bwameth_genome_preparation.output,
         fin1 = DIR_trimmed+"{sample}_1_val_1.fq.gz",
         fin2 = DIR_trimmed+"{sample}_2_val_2.fq.gz"
     output:
         bam = DIR_mapped+"{sample}.bwameth.bam"
     log:
         DIR_mapped+"{sample}_bwameth_pe_mapping.log"
     message: "Mapping paired-end reads to genome using bwa-meth."
     shell:
        nice("bwameth",["--reference {GENOMEFILE}",
            "{input.fin1} {input.fin2}","|",
            tool('samtools'), "view -bS - > {output.bam}"])

if not SUBSET_READS and NOTRIMMING:
  rule bwameth_align_pe_notrimming:
     input:
         index = rules.bwameth_genome_preparation.output,
         fin1 = PATHIN+"{sample}_1.fq.gz",
         fin2 = PATHIN+"{sample}_2.fq.gz"
     output:
         bam = DIR_mapped+"{sample}.bwameth.bam"
     params:
        # bwa-meth parameters
         bwameth_args = config['args']['bwameth']
     log:
         DIR_mapped+"{sample}_bwameth_pe_mapping.log"
     message: "Mapping paired-end reads to genome using bwa-meth."
     shell:
        nice("bwameth",["--reference {GENOMEFILE}",
            "{input.fin1} {input.fin2}","|",
            tool('samtools'), "view -bS - > {output.bam}"])


# # ==========================================================================================
# # Align subset of data sets
#    

if SUBSET_READS:
    # rule bwameth_align_pe:
  #    input:
  #        genomefile+".bwameth.c2t.sa",
  #        genomefile+".bwameth.c2t.amb",
  #        genomefile+".bwameth.c2t.ann",
  #        genomefile+".bwameth.c2t.pac",
  #        genomefile+".bwameth.c2t.bwt",
  #        genomefile+".bwameth.c2t",
  #        fin1 = DIR_trimmed_subset+"{sample}_1_val_1.fq.gz", ###########
  #        fin2 = DIR_trimmed_subset+"{sample}_2_val_2.fq.gz", ##############
  #    output:
  #        bam = DIR_mapped+"{sample}.bwameth.bam"
  #    params:
  #       # bwa-meth parameters
  #        bwameth_args = config['args']['bwameth']
  #    log:
  #        DIR_mapped+"{sample}_bwameth_pe_mapping.log"
  #    message: "Mapping paired-end reads to genome using bwa-meth."
  #    shell:
  #        """
  #        set -o pipefail
  #        {tools}/bwameth \\
  #        {params.bwameth_args} \\
  #        --reference {genomefile} \\
  #        {input.fin1} {input.fin2} | {tools}/samtools view -bS - > {output.bam}
  #       """
  
  rule bwameth_align_pe_raw:
     input:
         index = rules.bwameth_genome_preparation.output,
         fin1 = DIR_trimmed_subset+"{sample}_1.fq.gz", #####################
         fin2 = DIR_trimmed_subset+"{sample}_2.fq.gz", ####################
     output:
         bam = DIR_mapped+"{sample}.bwameth.bam"
     params:
        # bwa-meth parameters
         bwameth_args = config['args']['bwameth']
     log:
         DIR_mapped+"{sample}_bwameth_pe_mapping.log"
     message: "Mapping paired-end reads to genome using bwa-meth."
     shell:
        nice("bwameth",["--reference {GENOMEFILE}",
            "{input.fin1} {input.fin2}","|",
            tool('samtools'), "view -bS - > {output.bam}"],"{log}")
