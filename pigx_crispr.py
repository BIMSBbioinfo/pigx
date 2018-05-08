"""
Snakefile for pigx crispr pipeline
"""

import os
import yaml
import csv
import inspect


# tools
RSCRIPT = config['tools']['Rscript']

# input locations
SRC_DIR = config['source-dir']
READS_DIR = config['reads-dir']
ADAPTERS = config['adapters']
SAMPLE_SHEET_FILE = config['sample_sheet']

#output locations
OUTPUT_DIR = config['output-dir']
TRIMMED_READS_DIR = os.path.join(OUTPUT_DIR, 'trimmed_reads')
LOG_DIR           = os.path.join(OUTPUT_DIR, 'logs')
FASTQC_DIR        = os.path.join(OUTPUT_DIR, 'fastqc')
MULTIQC_DIR       = os.path.join(OUTPUT_DIR, 'multiqc')
MAPPED_READS_DIR  = os.path.join(OUTPUT_DIR, 'aln')
BEDGRAPH_DIR      = os.path.join(OUTPUT_DIR, 'bedgraph')
BBMAP_INDEX_DIR   = os.path.join(OUTPUT_DIR, 'bbmap_indexes')
BED_DIR           = os.path.join(OUTPUT_DIR, 'bed') 
REPORT_DIR        = os.path.join(OUTPUT_DIR, 'reports')

#other parameters
AMPLICONS = config.get('amplicons', {})
nodeN = config['nodeN']


## Load sample sheet
with open(SAMPLE_SHEET_FILE, 'r') as fp:
  rows =  [row for row in csv.reader(fp, delimiter=',')]
  header = rows[0]; rows = rows[1:]
  SAMPLE_SHEET = [dict(zip(header, row)) for row in rows]

# Convenience function to access fields of sample sheet columns that
# match the predicate.  The predicate may be a string.
def lookup(column, predicate, fields=[]):
  if inspect.isfunction(predicate):
    records = [line for line in SAMPLE_SHEET if predicate(line[column])]
  else:
    records = [line for line in SAMPLE_SHEET if line[column]==predicate]
  return [record[field] for record in records for field in fields]

SAMPLES = [line['sample_name'] for line in SAMPLE_SHEET]

def reads_input(wc):
  sample = wc.sample
  return [os.path.join(READS_DIR, f) for f in lookup('sample_name', sample, ['reads']) if f]

def get_amplicon_file(wc, file_type):
    path = AMPLICONS[wc.amplicon][file_type]
    return(path)

def get_output_file_list(DIR, ext):
    paths = list()
    for sample in SAMPLES:
        amplicon = lookup('sample_name', sample, ['amplicon'])[0]
        paths.append(os.path.join(DIR, amplicon, ".".join([sample, ext])))
    return(paths)


rule all:
    input:
        get_output_file_list(FASTQC_DIR, "fastqc.done"),
        get_output_file_list(TRIMMED_READS_DIR, "fastq.gz"),
        get_output_file_list(MAPPED_READS_DIR, "sam"),
        get_output_file_list(MAPPED_READS_DIR, "bam"), 
        get_output_file_list(MAPPED_READS_DIR, "bam.bai"), 
        get_output_file_list(MAPPED_READS_DIR, "samtools.stats.txt"),
        get_output_file_list(os.path.join(MAPPED_READS_DIR, "mpileup"), "mpileup.tsv"),                  
        get_output_file_list(os.path.join(MAPPED_READS_DIR, "mpileup"), "mpileup.counts.tsv"),  
        get_output_file_list(BEDGRAPH_DIR, "deletionScores.bedgraph"),
        get_output_file_list(BED_DIR, "deletions.bed"),
        os.path.join(OUTPUT_DIR, "multiqc", "multiqc_report.html"),
        expand(os.path.join(BBMAP_INDEX_DIR, "{amplicon}"), amplicon=AMPLICONS.keys()),
        expand(os.path.join(REPORT_DIR, "{amplicon}.report.html"), amplicon=AMPLICONS.keys())

rule fastqc:
    input: reads_input      
    output: os.path.join(FASTQC_DIR, "{amplicon}", "{sample}.fastqc.done")
    params:
        outdir=os.path.join(FASTQC_DIR, "{amplicon}"), 
        outfile=os.path.join(FASTQC_DIR, "{amplicon}", "{sample}.fastqc.done")
    log: os.path.join(LOG_DIR, "{amplicon}", "{sample}.fastqc.log")
    shell: "fastqc -o {params.outdir} {input} > {log} 2>&1; touch {params.outfile}"

rule trimmomatic:
    input: reads_input
    output: os.path.join(TRIMMED_READS_DIR, "{amplicon}", "{sample}.fastq.gz")
    log: os.path.join(LOG_DIR, "{amplicon}", "trimmomatic.{sample}.log")
    shell: "trimmomatic SE -threads {nodeN} {input} {output} \
       ILLUMINACLIP:{ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > {log} 2>&1"

rule bbmap_indexgenome:
    input: lambda wildcards: get_amplicon_file(wildcards, 'fasta')
    output: os.path.join(BBMAP_INDEX_DIR, "{amplicon}")
    log: os.path.join(LOG_DIR, 'bbmap_index_{amplicon}.log')
    shell: "bbmap.sh ref={input} path={output} > {log} 2>&1"

rule bbmap_map:
    input: 
        reads = os.path.join(TRIMMED_READS_DIR, "{amplicon}", "{sample}.fastq.gz"),
        ref = lambda wildcards: os.path.join(BBMAP_INDEX_DIR, lookup('sample_name', wildcards.sample, ['amplicon'])[0])
    output:
        os.path.join(MAPPED_READS_DIR, "{amplicon}", "{sample}.sam")
    log: os.path.join(LOG_DIR, "{amplicon}", "bbmap_{sample}.log")
    shell:
        "bbmap.sh path={input.ref} in={input.reads} outm={output} t={nodeN} sam=1.3 > {log} 2>&1"

rule samtools_sam2bam:
    input: os.path.join(MAPPED_READS_DIR, "{amplicon}", "{sample}.sam")
    output: os.path.join(MAPPED_READS_DIR, "{amplicon}", "{sample}.bam")
    log: os.path.join(LOG_DIR, "{amplicon}", "sam2bam_{sample}.log")
    shell: "samtools view -bh {input} | samtools sort -o {output} > {log} 2>&1"
        #"samtools view -bh {input} | samtools sort | samtools rmdup -s - {output} > {log} 2>&1"
                           
rule samtools_mpileup:
    input: os.path.join(MAPPED_READS_DIR, "{amplicon}", "{sample}.bam")
    output: os.path.join(MAPPED_READS_DIR, "mpileup", "{amplicon}", "{sample}.mpileup.tsv")
    log: os.path.join(LOG_DIR, "{amplicon}", "mpileup_{sample}.log")
    shell: "samtools mpileup {input} -o {output} > {log} 2>&1"

rule samtools_indexbam:
    input: os.path.join(MAPPED_READS_DIR, "{amplicon}", "{sample}.bam")
    output: os.path.join(MAPPED_READS_DIR, "{amplicon}", "{sample}.bam.bai")
    log: os.path.join(LOG_DIR, "{amplicon}", "samtools_index_{sample}.log")
    shell: "samtools index {input} > {log} 2>&1"

rule samtools_stats:
    input: 
        bamfile = os.path.join(MAPPED_READS_DIR, "{amplicon}", "{sample}.bam"),
        ref = lambda wildcards: get_amplicon_file(wildcards, 'fasta')
    output: os.path.join(MAPPED_READS_DIR, "{amplicon}", "{sample}.samtools.stats.txt")
    log: os.path.join(LOG_DIR, "{amplicon}", "samtools_stats.{sample}.log")
    shell: "samtools stats --reference {input.ref} {input.bamfile} > {output} 2> {log}"


rule parse_mpileup:
    input: os.path.join(MAPPED_READS_DIR, "mpileup", "{amplicon}", "{sample}.mpileup.tsv")
    output: os.path.join(MAPPED_READS_DIR, "mpileup", "{amplicon}", "{sample}.mpileup.counts.tsv")
    log: os.path.join(LOG_DIR, "{amplicon}", 'parse_mpileup.{sample}.log')
    params:
        script=os.path.join(SRC_DIR, "src", "parse_mpileup.py")
    shell: "python {params.script} {input} > {output} 2> {log}"

rule multiqc:
    input:
        fastqc = get_output_file_list(FASTQC_DIR, "fastqc.done"),
        trimmomatic = get_output_file_list(TRIMMED_READS_DIR, "fastq.gz"),
        samtools = get_output_file_list(MAPPED_READS_DIR, "samtools.stats.txt")
    output:
        os.path.join(OUTPUT_DIR, "multiqc", "multiqc_report.html")
    params: 
        analysis_folder = OUTPUT_DIR,
        output_folder = os.path.join(OUTPUT_DIR, "multiqc")
    log: os.path.join(LOG_DIR, 'multiqc.log')
    shell: "multiqc -o {params.output_folder} {params.analysis_folder} > {log} 2>&1"

rule extractDeletionProfiles:
    input: 
        bamIndex = os.path.join(MAPPED_READS_DIR, "{amplicon}", "{sample}.bam.bai"),
        bamFile = os.path.join(MAPPED_READS_DIR, "{amplicon}", "{sample}.bam"),
        mpileupOutput = os.path.join(MAPPED_READS_DIR, "mpileup", "{amplicon}", "{sample}.mpileup.tsv"), 
        parsedMpileupOutput = os.path.join(MAPPED_READS_DIR, "mpileup", "{amplicon}", "{sample}.mpileup.counts.tsv"),
        cutSitesFile = lambda wildcards: get_amplicon_file(wildcards, 'cutsites'),
    output: 
        os.path.join(BEDGRAPH_DIR, "{amplicon}", "{sample}.deletionScores.bedgraph")
    params:
        outdir=os.path.join(BEDGRAPH_DIR, "{amplicon}"),
        sgRNA_list = lambda wildcards: lookup('sample_name', wildcards.sample, ['sgRNA_ids'])[0],
        script=os.path.join(SRC_DIR, "src", "extractDeletionProfiles.R")
    log: os.path.join(LOG_DIR, "{amplicon}", "extractDeletionProfiles_{sample}.log")
    shell: "{RSCRIPT} {params.script} {input.bamFile} {input.mpileupOutput} {input.parsedMpileupOutput} {wildcards.sample} {params.outdir} {input.cutSitesFile} {params.sgRNA_list} > {log} 2>&1"
        
rule extractDeletionCoordinates:
    input:
        bamIndex = os.path.join(MAPPED_READS_DIR, "{amplicon}", "{sample}.bam.bai"),
        bamFile = os.path.join(MAPPED_READS_DIR, "{amplicon}", "{sample}.bam")
    params:
        outdir = os.path.join(BED_DIR, "{amplicon}"),
        script= os.path.join(SRC_DIR, "src", "extractDeletionCoordinates.R")
    output:
        os.path.join(BED_DIR, "{amplicon}", "{sample}.deletions.bed")
    log: os.path.join(LOG_DIR, "{amplicon}", "extractDeletionCoordinates_{sample}.log")
    shell:
        "{RSCRIPT} {params.script} {input.bamFile} {wildcards.sample} {params.outdir} > {log} 2>&1"
        
         
rule report:
  input:
    deletionScores = get_output_file_list(BEDGRAPH_DIR, "deletionScores.bedgraph")
  params:
    fasta = lambda wildcards: get_amplicon_file(wildcards, 'fasta'),
    cutSitesFile = lambda wildcards: get_amplicon_file(wildcards, 'cutsites'),
    reportR = os.path.join(SRC_DIR, "src", "runReport.R"),
    reportRmd = os.path.join(SRC_DIR, "src", "report.Rmd"),
    bedgraphFolder = os.path.join(BEDGRAPH_DIR, "{amplicon}")
  log: os.path.join(LOG_DIR, "{amplicon}.report.log")
  output:
    os.path.join(REPORT_DIR, '{amplicon}.report.html')
  shell:
    "{RSCRIPT} {params.reportR}  --reportFile={params.reportRmd} --ampliconFastaFile={params.fasta} --ampliconName={wildcards.amplicon} --cutSitesFile={params.cutSitesFile} --sampleSheetFile={SAMPLE_SHEET_FILE} --bedgraphFolder={params.bedgraphFolder} --workdir={REPORT_DIR} --prefix={wildcards.amplicon} > {log} 2>&1"        
    
    
    