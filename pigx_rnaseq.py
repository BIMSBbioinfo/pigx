"""
Snakefile for pigx rnaseq pipeline
"""

import os
import yaml
import csv
import inspect

GENOME_FASTA = config['locations']['genome-fasta']
READS_DIR = config['locations']['reads-folder']
OUTPUT_DIR = config['locations']['output-folder']
ORGANISM = config['organism']

print("ORGANISM:",ORGANISM)

SCRIPTS_DIR = os.path.join(config['locations']['pkglibexecdir'], 'scripts/')

TRIMMED_READS_DIR = os.path.join(OUTPUT_DIR, 'trimmed_reads')
LOG_DIR           = os.path.join(OUTPUT_DIR, 'logs')
FASTQC_DIR        = os.path.join(OUTPUT_DIR, 'fastqc')
MULTIQC_DIR       = os.path.join(OUTPUT_DIR, 'multiqc')
MAPPED_READS_DIR  = os.path.join(OUTPUT_DIR, 'mapped_reads')
BIGWIG_DIR        = os.path.join(OUTPUT_DIR, 'bigwig_files')
HTSEQ_COUNTS_DIR  = os.path.join(OUTPUT_DIR, 'feature_counts')
PREPROCESSED_OUT  = os.path.join(OUTPUT_DIR, 'preprocessed_data')

FASTQC_EXEC  = config['tools']['fastqc']['executable']
MULTIQC_EXEC = config['tools']['multiqc']['executable']
STAR_EXEC    = config['tools']['star']['executable']
STAR_THREADS = config['tools']['star']['n-threads']
TRIM_GALORE_EXEC = config['tools']['trim-galore']['executable']
BAMCOVERAGE_EXEC = config['tools']['bamCoverage']['executable']
SAMTOOLS_EXEC    = config['tools']['samtools']['executable']
HTSEQ_COUNT_EXEC = config['tools']['htseq-count']['executable']
RSCRIPT_EXEC     = config['tools']['R']['Rscript']


GTF_FILE = config['locations']['gtf-file']
SAMPLE_SHEET_FILE = config['locations']['sample-sheet']

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

SAMPLES = [line['name'] for line in SAMPLE_SHEET]

rule all:
  input: 
      report = os.path.join(OUTPUT_DIR, "report", "comparison1.deseq.report.html"),
      star_index_file = os.path.join(OUTPUT_DIR, 'star_index', "SAindex"),
      multiqc_report = os.path.join(MULTIQC_DIR, 'multiqc_report.html')
      
def trim_galore_input(args):
  sample = args[0]
  return [os.path.join(READS_DIR, f) for f in lookup('name', sample, ['reads', 'reads2']) if f]

rule trim_galore_pe:
  input: trim_galore_input
  output:
    r1=os.path.join(TRIMMED_READS_DIR, "{sample}_R1.fastq.gz"),
    r2=os.path.join(TRIMMED_READS_DIR, "{sample}_R2.fastq.gz")
  params:
    tmp1=lambda wildcards, output: os.path.join(TRIMMED_READS_DIR, lookup('name', wildcards[0], ['reads'])[0]).replace('.fastq.gz','_val_1.fq.gz'),
    tmp2=lambda wildcards, output: os.path.join(TRIMMED_READS_DIR, lookup('name', wildcards[0], ['reads2'])[0]).replace('.fastq.gz','_val_2.fq.gz')
  log: os.path.join(LOG_DIR, 'trim_galore_{sample}.log')
  shell: "{TRIM_GALORE_EXEC} -o {TRIMMED_READS_DIR} --paired {input[0]} {input[1]} >> {log} 2>&1 && sleep 10 && mv {params.tmp1} {output.r1} && mv {params.tmp2} {output.r2}"

rule trim_galore_se:
  input: trim_galore_input
  output: os.path.join(TRIMMED_READS_DIR, "{sample}_R.fastq.gz"),
  params: tmp=lambda wildcards, output: os.path.join(TRIMMED_READS_DIR, lookup('name', wildcards[0], ['reads'])[0]).replace('.fastq.gz','_trimmed.fq.gz'),
  log: os.path.join(LOG_DIR, 'trim_galore_{sample}.log')
  shell: "{TRIM_GALORE_EXEC} -o {TRIMMED_READS_DIR} {input[0]} >> {log} 2>&1 && sleep 10 && mv {params.tmp} {output}"

rule fastqc:
  input: os.path.join(MAPPED_READS_DIR, '{sample}_Aligned.sortedByCoord.out.bam')
  output: os.path.join(FASTQC_DIR, '{sample}_Aligned.sortedByCoord.out_fastqc.zip')
  log: os.path.join(LOG_DIR, 'fastqc_{sample}.log')
  shell: "{FASTQC_EXEC} -o {FASTQC_DIR} -f bam {input} >> {log} 2>&1"

rule star_index:
    input: GENOME_FASTA
    output: 
        star_index_dir = os.path.join(OUTPUT_DIR, 'star_index'),
        star_index_file = os.path.join(OUTPUT_DIR, 'star_index', "SAindex")
    log: os.path.join(LOG_DIR, 'star_index.log')
    shell: "{STAR_EXEC} --runThreadN {STAR_THREADS} --runMode genomeGenerate --genomeDir {output.star_index_dir} --genomeFastaFiles {input} --sjdbGTFfile {GTF_FILE} >> {log} 2>&1"

def star_map_input(args):
  sample = args[0]
  reads_files = [os.path.join(READS_DIR, f) for f in lookup('name', sample, ['reads', 'reads2']) if f]
  if len(reads_files) > 1:
    return [os.path.join(TRIMMED_READS_DIR, "{sample}_R1.fastq.gz".format(sample=sample)), os.path.join(TRIMMED_READS_DIR, "{sample}_R2.fastq.gz".format(sample=sample))]
  elif len(reads_files) == 1:
    return [os.path.join(TRIMMED_READS_DIR, "{sample}_R.fastq.gz".format(sample=sample))]

rule star_map:
  input:
    index_dir = rules.star_index.output.star_index_dir,
    reads = star_map_input
  output:
    os.path.join(MAPPED_READS_DIR, '{sample}_Aligned.sortedByCoord.out.bam'),
    os.path.join(MAPPED_READS_DIR, "{sample}_ReadsPerGene.out.tab")            
  params:
    output_prefix=os.path.join(MAPPED_READS_DIR, '{sample}_'),
  log: os.path.join(LOG_DIR, 'star_map_{sample}.log')
  shell: "{STAR_EXEC} --runThreadN {STAR_THREADS} --genomeDir {input.index_dir} --readFilesIn {input.reads} --readFilesCommand 'gunzip -c' --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.output_prefix} --quantMode TranscriptomeSAM GeneCounts >> {log} 2>&1"

rule index_bam:
  input: rules.star_map.output
  output: os.path.join(MAPPED_READS_DIR, '{sample}_Aligned.sortedByCoord.out.bai')
  log: os.path.join(LOG_DIR, 'samtools_index_{sample}.log')
  shell: "{SAMTOOLS_EXEC} index {input} {output} >> {log} 2>&1"

rule bamCoverage:
  input:
    bam=rules.star_map.output,
    bai=rules.index_bam.output
  output: os.path.join(BIGWIG_DIR, '{sample}.bw')
  log: os.path.join(LOG_DIR, 'bamCoverage_{sample}.log')
  shell: "{BAMCOVERAGE_EXEC} --bam {input.bam} -o {output} >> {log} 2>&1"

rule multiqc:
  input:
    star_output=expand(os.path.join(MAPPED_READS_DIR, '{sample}_Aligned.sortedByCoord.out.bam'), sample=SAMPLES),
    fastqc_output=expand(os.path.join(FASTQC_DIR, '{sample}_Aligned.sortedByCoord.out_fastqc.zip'), sample=SAMPLES),
  output: os.path.join(MULTIQC_DIR, 'multiqc_report.html')
  log: os.path.join(LOG_DIR, 'multiqc.log')
  shell: "{MULTIQC_EXEC} -o {MULTIQC_DIR} {OUTPUT_DIR} >> {log} 2>&1"

rule htseq_count:
  input: rules.star_map.output
  output: os.path.join(HTSEQ_COUNTS_DIR, "{sample}_counts.txt")
  log: os.path.join(LOG_DIR, "htseq-count_{sample}.log")
  shell: "{HTSEQ_COUNT_EXEC} -f bam -t exon -i gene_id {input} {GTF_FILE} 1> {output} 2>> {log}"

rule counts_from_STAR:
  input: expand(os.path.join(MAPPED_READS_DIR, "{sample}_ReadsPerGene.out.tab"), sample=SAMPLES)
  output: os.path.join(PREPROCESSED_OUT, "counts_from_STAR.tsv")
  shell: "{RSCRIPT_EXEC} {SCRIPTS_DIR}/counts_matrix_from_STAR.R {MAPPED_READS_DIR} {output}"

rule translate_sample_sheet_for_report:
  input: SAMPLE_SHEET_FILE
  output: os.path.join(os.getcwd(), "colData.tsv")
  shell: "{RSCRIPT_EXEC} {SCRIPTS_DIR}/translate_sample_sheet_for_report.R {input}"

  
CASE_SAMPLE_GROUPS = ','.join(set(lookup('comparison_factor', lambda x: x=='1', ['sample_type'])))
CASE_CONTROL_GROUPS = ','.join(set(lookup('comparison_factor', lambda x: x!='1', ['sample_type'])))
rule report:
  input:
    counts=str(rules.counts_from_STAR.output),
    coldata=str(rules.translate_sample_sheet_for_report.output)
  params:
    outdir=os.path.join(OUTPUT_DIR, "report"),
    reportR=os.path.join(SCRIPTS_DIR, "runDeseqReport.R"),
    reportRmd=os.path.join(SCRIPTS_DIR, "deseqReport.Rmd")
  log: os.path.join(LOG_DIR, "report.log")
  output: os.path.join(OUTPUT_DIR, "report", "comparison1.deseq.report.html")
  shell: "{RSCRIPT_EXEC} {params.reportR} --reportFile={params.reportRmd} --countDataFile={input.counts} --colDataFile={input.coldata} --caseSampleGroups='{CASE_SAMPLE_GROUPS}' --controlSampleGroups='{CASE_CONTROL_GROUPS}' --workdir={params.outdir} --species='{ORGANISM}' --geneSetsFolder='' >> {log} 2>&1"
      
