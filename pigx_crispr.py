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
INDELS_DIR      = os.path.join(OUTPUT_DIR, 'indels')
BBMAP_INDEX_DIR   = os.path.join(OUTPUT_DIR, 'bbmap_indexes')
REPORT_DIR        = os.path.join(OUTPUT_DIR, 'reports')

#other parameters
AMPLICONS = config.get('amplicons', {})

COMPARISONS = config.get('comparisons', {})

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
        ancient(os.path.join(OUTPUT_DIR, 'comparisons.case-control.tsv')),
        ancient(os.path.join(OUTPUT_DIR, 'comparisons.time-series.tsv')),
        ancient(get_output_file_list(FASTQC_DIR, "fastqc.done")),
        ancient(get_output_file_list(TRIMMED_READS_DIR, "fastq.gz")),
        #get_output_file_list(MAPPED_READS_DIR, "sam"),
        ancient(get_output_file_list(MAPPED_READS_DIR, "bam")), 
        ancient(get_output_file_list(MAPPED_READS_DIR, "bam.bai")), 
        ancient(get_output_file_list(MAPPED_READS_DIR, "samtools.stats.txt")),
        ancient(get_output_file_list(INDELS_DIR, "indelScores.bedgraph")),
        ancient(get_output_file_list(INDELS_DIR, "coverageStats.tsv")),
        ancient(get_output_file_list(INDELS_DIR, "indel_stats_at_cutsites.tsv")),
        ancient(get_output_file_list(INDELS_DIR, "deletions.bed")),
        ancient(get_output_file_list(INDELS_DIR, "insertions.bed")),
        ancient(get_output_file_list(INDELS_DIR, "indels.unfiltered.tsv")),
        ancient(os.path.join(OUTPUT_DIR, "multiqc", "multiqc_report.html")),
        ancient(expand(os.path.join(BBMAP_INDEX_DIR, "{amplicon}"), amplicon=AMPLICONS.keys())),
        ancient(expand(os.path.join(REPORT_DIR, "{amplicon}.report.html"), amplicon=AMPLICONS.keys())),
        ancient(expand(os.path.join(REPORT_DIR, 'comparisons', '{amplicon}.report.comparisons.html'), amplicon=COMPARISONS.keys()))

rule get_comparisons:
    output: 
        os.path.join(OUTPUT_DIR, 'comparisons.case-control.tsv'),
        os.path.join(OUTPUT_DIR, 'comparisons.time-series.tsv')
    run:        
        f1 = open(os.path.join(OUTPUT_DIR, 'comparisons.case-control.tsv'), 'w')
        f1.write(''.join(['\t'.join(['amplicon', 'comparison', 'case_samples', 'control_samples']), '\n']))
        
        f2 = open(os.path.join(OUTPUT_DIR, 'comparisons.time-series.tsv'), 'w') 
        f2.write(''.join(['\t'.join(['amplicon', 'comparison', 'sample_series']), '\n']))
        
        for amplicon, comp_dict in COMPARISONS.items():
            for comp, d  in comp_dict.items():
                if('case_samples' in d.keys() and 'control_samples' in d.keys()):
                    f1.write(''.join(['\t'.join([amplicon, comp, d['case_samples'], d['control_samples']]), '\n']))
                if('time_series' in d.keys()):
                    f2.write(''.join(['\t'.join([amplicon, comp, d['time_series']]), '\n']))
        f1.close()
        f2.close()


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

rule getIndelStats:
    input: 
        bamIndex = os.path.join(MAPPED_READS_DIR, "{amplicon}", "{sample}.bam.bai"),
        bamFile = os.path.join(MAPPED_READS_DIR, "{amplicon}", "{sample}.bam"),
        cutSitesFile = lambda wildcards: get_amplicon_file(wildcards, 'cutsites'),
    output: 
        os.path.join(INDELS_DIR, "{amplicon}", "{sample}.indelScores.bedgraph"),
        os.path.join(INDELS_DIR, "{amplicon}", "{sample}.coverageStats.tsv"),
        os.path.join(INDELS_DIR, "{amplicon}", "{sample}.indel_stats_at_cutsites.tsv"),
        os.path.join(INDELS_DIR, "{amplicon}", "{sample}.deletions.bed"),
        os.path.join(INDELS_DIR, "{amplicon}", "{sample}.insertions.bed"),
        os.path.join(INDELS_DIR, "{amplicon}", "{sample}.indels.unfiltered.tsv")
    params:
        outdir=os.path.join(INDELS_DIR, "{amplicon}"),
        sgRNA_list = lambda wildcards: lookup('sample_name', wildcards.sample, ['sgRNA_ids'])[0],
        script=os.path.join(SRC_DIR, "src", "getIndelStats.R")
    log: os.path.join(LOG_DIR, "{amplicon}", "getIndelStats_{sample}.log")
    shell: "{RSCRIPT} {params.script} {input.bamFile} {wildcards.sample} {params.outdir} {input.cutSitesFile} {params.sgRNA_list} > {log} 2>&1"
        
          
rule report:
    input:
        coverageStats = get_output_file_list(INDELS_DIR, "coverageStats.tsv"),
        cutsiteStats = get_output_file_list(INDELS_DIR, "indel_stats_at_cutsites.tsv"),
        indels = get_output_file_list(INDELS_DIR, "indels.unfiltered.tsv")
    params:
      fasta = lambda wildcards: get_amplicon_file(wildcards, 'fasta'),
      cutSitesFile = lambda wildcards: get_amplicon_file(wildcards, 'cutsites'),
      reportR = os.path.join(SRC_DIR, "src", "runReport.R"),
      reportRmd = os.path.join(SRC_DIR, "src", "report.Rmd"),
      indelsFolder = os.path.join(INDELS_DIR, "{amplicon}")
    log: os.path.join(LOG_DIR, "{amplicon}.report.log")
    output:
        os.path.join(REPORT_DIR, '{amplicon}.report.html')
    shell:
        "{RSCRIPT} {params.reportR}  --reportFile={params.reportRmd} --ampliconFastaFile={params.fasta} --ampliconName={wildcards.amplicon} --cutSitesFile={params.cutSitesFile} --sampleSheetFile={SAMPLE_SHEET_FILE} --indelsFolder={params.indelsFolder} --workdir={REPORT_DIR} --prefix={wildcards.amplicon} > {log} 2>&1"        
    

rule report_comparisons:
    input:
        coverageStats = get_output_file_list(INDELS_DIR, "coverageStats.tsv"),
        cutsiteStats = get_output_file_list(INDELS_DIR, "indel_stats_at_cutsites.tsv"),
        indels = get_output_file_list(INDELS_DIR, "indels.unfiltered.tsv"),
        comp = os.path.join(OUTPUT_DIR, 'comparisons.case-control.tsv')
    params:
        reportR = os.path.join(SRC_DIR, "src", "runReport.comparisons.R"),
        reportRmd = os.path.join(SRC_DIR, "src", "report.comparisons.Rmd"),
        indelsFolder = os.path.join(INDELS_DIR, "{amplicon}"),
        outDir = os.path.join(REPORT_DIR, 'comparisons'),
    log: os.path.join(LOG_DIR, "{amplicon}.report.comparisons.log")
    output:
        os.path.join(REPORT_DIR, 'comparisons', '{amplicon}.report.comparisons.html')
    shell:
        "{RSCRIPT} {params.reportR}  --reportFile={params.reportRmd} --ampliconName={wildcards.amplicon} --comparisonsFile={input.comp} --indelsFolder={params.indelsFolder} --workdir={params.outDir} --prefix={wildcards.amplicon} > {log} 2>&1"        
    

       
         
 
    
        
    