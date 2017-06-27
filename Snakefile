configfile: "config.yaml"

nodeN = config["nodeN"]
adapters = config["adapters"]

rule all:
    input:
        expand("fastqc/{sample}_fastqc.html", sample = config["samples"]),
        "ref",
        expand("mpileup/{sample}.mpileup.counts.tsv", sample = config["samples"]),
        expand("aln/{sample}.deduped.bam.bai", sample = config["samples"])

rule fastqc:
    input:
        "sample_data/raw_reads/{sample}.fastq.gz"
    output:
        "fastqc/{sample}_fastqc.html"
    shell:
        "fastqc {input} -o fastqc"

rule trimmomatic:
    input:
        "sample_data/raw_reads/{sample}.fastq.gz"
    output:
        "sample_data/filtered_reads/{sample}.fastq.gz"
    shell:
       "trimmomatic SE -threads {nodeN} {input} {output} \
       ILLUMINACLIP:{adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

rule bbmap_indexgenome:
    input:
        fa=config["genome"]
    output:
        "ref"
    shell:
        "bbmap.sh ref={input.fa}"

rule bbmap_map:
    input:
        "sample_data/filtered_reads/{sample}.fastq.gz",
        "ref"
    output:
        "aln/{sample}.sam"
    shell:
        "bbmap.sh in={input} outm={output} t={nodeN} sam=1.3"

rule samtools_sam2bam:
    input:
        "aln/{sample}.sam"
    output:
        "aln/{sample}.bam"
    shell:
        "samtools view -bh {input} > {output}"

rule samtools_sort:
    input:
        "aln/{sample}.bam"
    output:
        "aln/{sample}.sorted.bam"
    shell:
        "samtools sort {input} > {output}"

rule samtools_rmdup_SE:
    input:
        "aln/{sample}.sorted.bam"
    output:
        "aln/{sample}.deduped.bam"
    shell:
        "samtools rmdup -s {input} {output}"

rule samtools_index:
    input:
        "aln/{sample}.deduped.bam"
    output:
        "aln/{sample}.deduped.bam.bai"
    shell:
        "samtools index {input}"

rule samtools_mpileup:
    input:
        "aln/{sample}.deduped.bam"
    output:
        "mpileup/{sample}.mpileup.tsv"
    shell:
        "samtools mpileup {input} > {output}"

rule parse_mpileup:
    input:
        "mpileup/{sample}.mpileup.tsv"
    output:
        "mpileup/{sample}.mpileup.counts.tsv"
    shell:
        "python src/parse_mpileup.py {input} > {output}"


#rule extractPerBaseDeletionScores

#rule getDeletions

#rule report
