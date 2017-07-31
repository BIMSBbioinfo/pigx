"""
Snakefile for pigx rnaseq pipeline
"""

import os
import yaml
import pandas as pd


## Load Settings
with open("settings.yaml") as f:
  SETTINGS = yaml.load(f)

GENOME_DIR = SETTINGS['locations']['genome-folder']
READS_DIR = SETTINGS['locations']['reads-folder']
OUTPUT_DIR = SETTINGS['locations']['output-folder']

TRIMMED_READS_DIR = os.path.join(OUTPUT_DIR, 'trimmed_reads')
LOG_DIR = os.path.join(OUTPUT_DIR, 'logs')
FASTQC_DIR = os.path.join(OUTPUT_DIR, 'fastqc')
MULTIQC_DIR = os.path.join(OUTPUT_DIR, 'multiqc')
MAPPED_READS_DIR = os.path.join(OUTPUT_DIR, 'mapped_reads')


FASTQC_EXEC  = SETTINGS['tools']['fastqc']['executable']
MULTIQC_EXEC = SETTINGS['tools']['multiqc']['executable']
STAR_EXEC    = SETTINGS['tools']['star']['executable']
STAR_THREADS = SETTINGS['tools']['star']['n-threads']
TRIM_GALORE_EXEC = SETTINGS['tools']['trim-galore']['executable']
TRIM_GALORE_ARGS = SETTINGS['tools']['trim-galore']['args']


## Load sample sheet
SAMPLE_SHEET = pd.read_csv("sample_sheet.csv")
SAMPLES = SAMPLE_SHEET['name']

rule all:
  input: os.path.join(MULTIQC_DIR, 'multiqc_report.html')

# TODO
# 1. - add bigwig file generation
# 2. - more counting (featureCounts, HTSeq)
# 3. - make counts matrix
# 4. - integrate Bora's report

## Target rule - multiqc report after everything has been mapped
rule multiqc:
  input:
    star_output=expand(os.path.join(MAPPED_READS_DIR, '{sample}_Aligned.sortedByCoord.out.bam'), sample=SAMPLES),
    fastqc_output=expand(os.path.join(FASTQC_DIR, '{sample}_Aligned.sortedByCoord.out_fastqc.zip'), sample=SAMPLES),
  output: os.path.join(MULTIQC_DIR, 'multiqc_report.html')
  log: os.path.join(LOG_DIR, 'multiqc.log')
  shell: "{MULTIQC_EXEC} -o {MULTIQC_DIR} {OUTPUT_DIR} >> {log} 2>&1"

rule star_map:
  input:
    r1=os.path.join(TRIMMED_READS_DIR, "{sample}_R1.fastq.gz"),
    r2=os.path.join(TRIMMED_READS_DIR, "{sample}_R2.fastq.gz")
  output: os.path.join(MAPPED_READS_DIR, '{sample}_Aligned.sortedByCoord.out.bam')
  params:
    output_prefix=os.path.join(MAPPED_READS_DIR, '{sample}_')
  log: os.path.join(LOG_DIR, 'star_map_{sample}.log')
  shell: "{STAR_EXEC} --runThreadN {STAR_THREADS} --genomeDir {GENOME_DIR} --readFilesIn {input} --readFilesCommand 'gunzip -c' --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.output_prefix} --quantMode TranscriptomeSAM GeneCounts >> {log} 2>&1"

rule fastqc:
  input: os.path.join(MAPPED_READS_DIR, '{sample}_Aligned.sortedByCoord.out.bam')
  output: os.path.join(FASTQC_DIR, '{sample}_Aligned.sortedByCoord.out_fastqc.zip')
  log: os.path.join(LOG_DIR, 'fastqc_{sample}.log')
  shell: "{FASTQC_EXEC} -o {FASTQC_DIR} -f bam {input} >> {log} 2>&1"


def trim_galore_input(args):
  sample = args[0]
  return [os.path.join(READS_DIR, f) for f in SAMPLE_SHEET[SAMPLE_SHEET['name']==sample][['reads', 'reads2']].iloc[0]]


rule trim_galore:
  input: trim_galore_input
  output:
    r1=os.path.join(TRIMMED_READS_DIR, "{sample}_R1.fastq.gz"),
    r2=os.path.join(TRIMMED_READS_DIR, "{sample}_R2.fastq.gz")
  params:
    tmp1=lambda wildcards, output: os.path.join(TRIMMED_READS_DIR, SAMPLE_SHEET[SAMPLE_SHEET['name']==wildcards[0]]['reads'].iloc[0]).replace('.fastq.gz','_val_1.fq.gz'),
    tmp2=lambda wildcards, output: os.path.join(TRIMMED_READS_DIR, SAMPLE_SHEET[SAMPLE_SHEET['name']==wildcards[0]]['reads2'].iloc[0]).replace('.fastq.gz','_val_2.fq.gz')
  log: os.path.join(LOG_DIR, 'trim_galore_{sample}.log')
  shell: "{TRIM_GALORE_EXEC} -o {TRIMMED_READS_DIR} {TRIM_GALORE_ARGS} {input[0]} {input[1]} >> {log} 2>&1 && sleep 10 && mv {params.tmp1} {output.r1} && mv {params.tmp2} {output.r2}"

