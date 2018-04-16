"""
Snakefile for pigx crispr pipeline
"""

import os
import yaml
import csv
import inspect


AMPLICONS = config.get('amplicons', {})
AMPLICONS_FASTA = [AMPLICONS[amp]['fasta'] for amp in AMPLICONS.keys()]
READS_DIR = config['reads-dir']
OUTPUT_DIR = config['output-dir']
SRC_DIR = config['source-dir']
print('SRC_DIR:',SRC_DIR)

TRIMMED_READS_DIR = os.path.join(OUTPUT_DIR, 'trimmed_reads')
LOG_DIR           = os.path.join(OUTPUT_DIR, 'logs')
FASTQC_DIR        = os.path.join(OUTPUT_DIR, 'fastqc')
MULTIQC_DIR       = os.path.join(OUTPUT_DIR, 'multiqc')
MAPPED_READS_DIR  = os.path.join(OUTPUT_DIR, 'aln')
BEDGRAPH_DIR      = os.path.join(OUTPUT_DIR, 'bedgraph')
BBMAP_INDEX_DIR   = os.path.join(OUTPUT_DIR, 'bbmap_indexes') 

SAMPLE_SHEET_FILE = config['sample_sheet']

nodeN = config['nodeN']
ADAPTERS = config['adapters']


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
print('SAMPLES', SAMPLES)

def reads_input(args):
  sample = args[0]
  return [os.path.join(READS_DIR, f) for f in lookup('sample_name', sample, ['reads']) if f]

def get_amplicon_fasta(args):
  amp = args[0]
  print("looking for ",amp)
  return [AMPLICONS[amp]['fasta']]

rule all:
    input:
        #expand(os.path.join(OUTPUT_DIR, "fastqc", "{sample}_fastqc.html"), sample = SAMPLES),
        expand(os.path.join(MAPPED_READS_DIR, "mpileup", "{sample}.mpileup.counts.tsv"), sample = SAMPLES),
        expand(os.path.join(BBMAP_INDEX_DIR, "{amplicon}"), amplicon=AMPLICONS.keys())
        #expand(os.path.join(MAPPED_READS_DIR, "{sample}.rmdup.bam"), sample = SAMPLES)

rule fastqc:
    input: reads_input
    output: os.path.join(OUTPUT_DIR, "fastqc", "{sample}_fastqc.html")
    log: os.path.join(LOG_DIR, 'fastqc_{sample}.log')
    shell: "fastqc -o {FASTQC_DIR} {input} >> {log} 2>&1"

#rule trimmomatic:
#    input:
#        "sample_data/raw_reads/{sample}.fastq.gz"
#    output:
#        "filtered_reads/{sample}.fastq.gz"
#    shell:
#       "trimmomatic SE -threads {nodeN} {input} {output} \
#       ILLUMINACLIP:{adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

rule bbmap_indexgenome:
    input: get_amplicon_fasta
    output: os.path.join(BBMAP_INDEX_DIR, "{amplicon}")
    #params: 
     #   path=os.path.join(BBMAP_INDEX_DIR, "{wildcards.amplicon}")
    log: os.path.join(LOG_DIR, 'bbmap_index_{amplicon}.log')
#    shell: "bbmap.sh ref={input} path={params.path} >> {log} 2>&1"
    shell: "bbmap.sh ref={input} path={output} >> {log} 2>&1"

rule bbmap_map:
    input: 
        reads = reads_input,
        ref = lambda wildcards: os.path.join(BBMAP_INDEX_DIR, lookup('sample_name', wildcards[0], ['amplicon'])[0])
    output:
        os.path.join(MAPPED_READS_DIR, "{sample}.sam")
    log: os.path.join(LOG_DIR, 'bbmap_{sample}.log')
    shell:
        "bbmap.sh path={input.ref} in={input.reads} outm={output} t={nodeN} sam=1.3 >> {log} 2>&1"

rule samtools_sam2bam:
    input: os.path.join(MAPPED_READS_DIR, "{sample}.sam")
    output: os.path.join(MAPPED_READS_DIR, "{sample}.bam")
    log: os.path.join(LOG_DIR, 'sam2bam_{sample}.log')
    shell:
        "samtools view -bh {input} | samtools sort | samtools rmdup -s - {output} >> {log} 2>&1"

rule samtools_mpileup:
    input: os.path.join(MAPPED_READS_DIR, "{sample}.bam")
    output: os.path.join(MAPPED_READS_DIR, "mpileup", "{sample}.mpileup.tsv")
    log: os.path.join(LOG_DIR, 'mpileup_{sample}.log')
    shell: "samtools mpileup {input} -o {output} >> {log} 2>&1"

rule parse_mpileup:
    input: os.path.join(MAPPED_READS_DIR, "mpileup", "{sample}.mpileup.tsv")
    output: os.path.join(MAPPED_READS_DIR, "mpileup", "{sample}.mpileup.counts.tsv")
    log: os.path.join(LOG_DIR, 'parse_mpileup.{sample}.log')
    params:
        script=os.path.join(SRC_DIR, "src", "parse_mpileup.py")
    shell: "python {params.script} {input} > {output} 2> {log}"
#
#
##rule extractPerBaseDeletionScores
#
##rule getDeletions
#
##rule report
