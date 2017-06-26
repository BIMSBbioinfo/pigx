configfile: "config.yaml"

#rule all

#rule fastqc

#rule trimmomatic:

#rule bbmap_indexgenome

#rule bbmap_map
#    input:
#        expand("sample_data/{sample}.fastq.gz", sample=config["samples"])
#    output:
#        "mapped_reads/"

#rule samtools_sam2bam

#rule samtools_rmdup

#rule samtools_sort

#rule samtools_index

#rule samtools_mpileup

#rule parse_mpileup #find indels

#rule extractPerBaseDeletionScores

#rule getDeletions

#rule report
