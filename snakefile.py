# PiGx SARS CoV2 wastewater sequencing pipeline
#
# Copyright © 2021 Akalin lab.
#
# This file is part of the PiGx SARS-CoV2 wastewater sequencing pipeline.
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

"""
Snakefile for PiGx SARS CoV2 wastewater sequencing pipeline
"""

import os
import csv
from itertools import chain

# TODO: 
# 1. obtain things from config file
# 2. make dummy samples and put them to tests/sample_data/reads
READS_DIR = 'tests/sample_data/reads'
REFERENCE_FASTA = 'tests/sample_data/NC_045512.2.fasta'
KRAKEN_DB = 'tests/sample_data/kraken_db.db' # download first, where?
SAMPLE_SHEET_CSV = 'tests/sample_sheet.csv'

OUTPUT_DIR = 'output'
TRIMMED_READS_DIR = os.path.join(OUTPUT_DIR, 'trimmed_reads')
LOG_DIR = os.path.join(OUTPUT_DIR, 'logs')
FASTQC_DIR = os.path.join(OUTPUT_DIR, 'fastqc')
MULTIQC_DIR = os.path.join(OUTPUT_DIR, 'multiqc')
MAPPED_READS_DIR = os.path.join(OUTPUT_DIR, 'mapped_reads')
VARIANTS_DIR = os.path.join(OUTPUT_DIR, 'variants')
KRAKEN_DIR = os.path.join(OUTPUT_DIR, 'kraken')
REPORT_DIR = os.path.join(OUTPUT_DIR, 'report')

SCRIPTS_DIR = 'scripts'


# Load sample sheet
with open(SAMPLE_SHEET_CSV, 'r') as fp:
  sample_sheet = list(csv.DictReader(fp))

SAMPLES = [row['name'] for row in sample_sheet]

targets = {
    'help': {
        'description': "Print all rules and their descriptions.",
        'files': []
    },
    'final_reports': {
        'description': "Produce a comprehensive report. This is the default target.",
        'files': (
            expand(os.path.join(VARIANTS_DIR, '{sample}_snv.csv'), sample=SAMPLES) + 
            expand(os.path.join(REPORT_DIR, '{sample}_variant_report.html'), sample=SAMPLES) + 
            expand(os.path.join(REPORT_DIR, '{sample}_kraken_report.html'), sample=SAMPLES)
        )
    }
}
selected_targets = ['final_reports']
OUTPUT_FILES = list(chain.from_iterable([targets[name]['files'] for name in selected_targets]))


rule all:
    input: OUTPUT_FILES


# Record any existing output files, so that we can detect if they have
# changed.
expected_files = {}
onstart:
    if OUTPUT_FILES:
        for name in OUTPUT_FILES:
            if os.path.exists(name):
                expected_files[name] = os.path.getmtime(name)


# Print generated target files.
onsuccess:
    if OUTPUT_FILES:
        # check if any existing files have been modified
        generated = []
        for name in OUTPUT_FILES:
            if name not in expected_files or os.path.getmtime(name) != expected_files[name]:
                generated.append(name)
        if generated:
            print("The following files have been generated:")
            for name in generated:
                print("  - {}".format(name))


# TODO: get read length and specify cut off
rule prinseq:
    input:
        r1 = os.path.join(READS_DIR, "{sample}_R1.fastq"),
        r2 = os.path.join(READS_DIR, "{sample}_R2.fastq")
    output:
        r1 = os.path.join(TRIMMED_READS_DIR, "{sample}_R1.fastq"),
        r2 = os.path.join(TRIMMED_READS_DIR, "{sample}_R2.fastq")
    params:
        len_cutoff = int(150 * 0.8) # read length * pct_cutoff
    log: os.path.join(LOG_DIR, 'prinseq_{sample}.log')
    shell: "prinseq-lite.pl -fastq {input.r1} -fastq2 {input.r2} -ns_max_n 4 -min_qual_mean 30 -trim_qual_left 30 -trim_qual_right 30 -trim_qual_window 10 -out_format 3 -out_good {TRIMMED_READS_DIR} -out_bad null -min_len {params.len_cutoff} >> {log} 2>&1"


rule bwa_index:
    input: REFERENCE_FASTA
    output: "{}.bwt".format(REFERENCE_FASTA)
    log: os.path.join(LOG_DIR, 'bwa_index.log')
    shell: "bwa index {input} >> {log} 2>&1"


rule bwa_align:
    input:
        fastq = [os.path.join(TRIMMED_READS_DIR, "{sample}_R1.fastq"), os.path.join(TRIMMED_READS_DIR, "{sample}_R2.fastq")]
        ref = REFERENCE_FASTA
        index = "{}.bwt".format(REFERENCE_FASTA)
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.sam')
    params:
        threads = 4
    log: os.path.join(LOG_DIR, 'bwa_align_{sample}.log')
    shell: "bwa mem -t {params.threads} {input.ref} {input.fastq} > {output} >> {log} 2>&1"


# TODO: also get unaligned .sam --> so far tmp.sam is used
# check what tmp actually is -- are that really the unaligned ones?
rule samtools_filter:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.sam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned.sam')
    log: os.path.join(LOG_DIR, 'samtools_filter_{sample}.log')
    shell: "samtools view -h -f 2 -F 2048 {input} > {output} >> {log} 2>&1"


rule samtools_sam2bam_A:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned.sam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned.bam')
    log: os.path.join(LOG_DIR, 'samtools_sam2bam_{sample}.log')
    shell: "samtools view -bS {input} > {output} >> {log} 2>&1"


# TODO: merge rule A and B and possibly? use unaligned instead of tmp
rule samtools_sam2bam_B:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.sam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.bam')
    log: os.path.join(LOG_DIR, 'samtools_sam2bam_B_{sample}.log')
    shell: "samtools view -bS {input} > {output} >> {log} 2>&1"


rule samtools_sort:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned.bam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam')
    log: os.path.join(LOG_DIR, 'samtools_sort_{sample}.log')
    shell: "samtools sort -o {output} {input} >> {log} 2>&1"


rule samtools_index:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bai')
    log: os.path.join(LOG_DIR, 'samtools_index_{sample}.log')
    shell: "samtools index {input} {output} >> {log} 2>&1"


rule fastqc:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam')
    output: os.path.join(FASTQC_DIR, '{sample}_aligned_sorted.fastqc.zip')
    log: os.path.join(LOG_DIR, 'fastqc_{sample}.log')
    params:
        out_dir = os.path.join(FASTQC_DIR)
    shell: "fastqc -o {params.out_dir} -f bam {input} >> {log} 2>&1"


# TODO: check if --call-indels is needed?
rule lofreq:
    input:
        aligned_bam = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam')
        ref = REFERENCE_FASTA
    output: os.path.join(VARIANTS_DIR, '{sample}_snv.vcf')
    log: os.path.join(LOG_DIR, 'lofreq_{sample}.log')
    shell: "lofreq call -f {input.ref} -o {output} --verbose {input.aligned_bam} >> {log} 2>&1"


rule vcf2csv:
    input: os.path.join(VARIANTS_DIR, '{sample}_snv.vcf')
    output: os.path.join(VARIANTS_DIR, '{sample}_snv.csv')
    log: os.path.join(LOG_DIR, 'vcf2csv_{sample}.log')
    shell: "python {SCRIPTS_DIR}/vcfTocsv.py {input} >> {log} 2>&1"


# TODO: add vep rule
# rule vep:
#     input: 
#     output: 
#     log: 
#     shell: ""


# TODO: fill in blanks, why sample_dir? where are sigmuts coming from?
# rule variant_report:
#     input:
#         vep_txt = os.path.join(VARIANTS_DIR, '{sample}_vep.txt')
#         snv_csv = os.path.join(VARIANTS_DIR, '{sample}_snv.csv')
#     output: os.path.join(REPORT_DIR, '{sample}_variant_report.html')
#     log: os.path.join(LOG_DIR, 'vcf2csv_{sample}.log')
#     shell: "Rscript {SCRIPTS_DIR}/run_variant_report.R --reportFile={SCRIPTS_DIR}/variantreport_p_sample.rmd --vep_txt_file={input.vep_txt} --snv_csv_file={input.snv_csv} --location_sigmuts=??? --sample_dir=??? --sample_name={sample} >> {log} 2>&1"


# TODO: check command/output as in miro board R1 & R2 is mentioned
# rule bam2fastq:
#     input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.bam')
#     output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.fastq')
#     log: os.path.join(LOG_DIR, 'bam2fastq_{sample}.log')
#     shell: "samtools fastq {input} > {output} >> {log} 2>&1"


# TODO: check command as in miro board there where R1 & R2 inputs
# rule kraken:
#     input:
#         unaligned_fastq = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.fastq')
#         database = KRAKEN_DB
#     output: os.path.join(KRAKEN_DIR, '{sample}_classified_unaligned_reads.txt')
#     log: os.path.join(LOG_DIR, 'kraken_{sample}.log')
#     shell: "kraken2 --output {output} --db {input.database} --paired {input.unaligned_fastq} >> {log} 2>&1"


# TODO: adapt command and possibly input due to 'taxonomy_folder'
# What is meant by taxonomy_folder??
# rule kraken_report:
#     input: os.path.join(KRAKEN_DIR, '{sample}_classified_unaligned_reads.txt')
#     output: os.path.join(REPORT_DIR, '{sample}_kraken_report.html')
#     log: os.path.join(LOG_DIR, 'kraken_report_{sample}.log')
#     shell: "ktImportTaxonomy -m 3 -t 5 {input} -tax [taxonomy folder]??? -o {output}"