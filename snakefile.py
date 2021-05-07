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
# 3. align overall appearance
READS_DIR = os.path.join(os.getcwd(), 'tests/sample_data/reads')
REFERENCE_FASTA = os.path.join(os.getcwd(), 'tests/sample_data/NC_045512.2.fasta')
AMPLICONS_BED = os.path.join(os.getcwd(), 'tests/sample_data/nCoV-2019_NCref.bed')
KRAKEN_DB = os.path.join(os.getcwd(), 'tests/databases/kraken_db')
KRONA_DB = os.path.join(os.getcwd(), 'tests/databases/krona_db')
SIGMUT_DB = os.path.join(os.getcwd(), 'tests/databases/sigmut_db')
VEP_DB = os.path.join(os.getcwd(), 'tests/databases/vep_db')
SAMPLE_SHEET_CSV = os.path.join(os.getcwd(), 'tests/sample_sheet.csv')

OUTPUT_DIR = os.path.join(os.getcwd(), 'output')

TRIMMED_READS_DIR = os.path.join(OUTPUT_DIR, 'trimmed_reads')
LOG_DIR = os.path.join(OUTPUT_DIR, 'logs')
FASTQC_DIR = os.path.join(OUTPUT_DIR, 'fastqc')
MULTIQC_DIR = os.path.join(OUTPUT_DIR, 'multiqc')
MAPPED_READS_DIR = os.path.join(OUTPUT_DIR, 'mapped_reads')
VARIANTS_DIR = os.path.join(OUTPUT_DIR, 'variants')
KRAKEN_DIR = os.path.join(OUTPUT_DIR, 'kraken')
COVERAGE_DIR = os.path.join(OUTPUT_DIR, 'coverage')
REPORT_DIR = os.path.join(OUTPUT_DIR, 'report')

SCRIPTS_DIR = os.path.join(os.getcwd(), 'scripts')


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
            expand(os.path.join(REPORT_DIR, '{sample}_qc_report.html'), sample=SAMPLES) + 
            expand(os.path.join(REPORT_DIR, '{sample}_variant_report.html'), sample=SAMPLES) + 
            expand(os.path.join(REPORT_DIR, '{sample}_kraken_report.html'), sample=SAMPLES)
        )
    },
    'lofreq': {
        'description': "Call variants and produce .vcf file and overview .csv file.",
        'files': (
            expand(os.path.join(VARIANTS_DIR, '{sample}_snv.csv'), sample=SAMPLES)
        )
    },
    'variant_reports': {
        'description': "Make variants reports.",
        'files': (
            expand(os.path.join(REPORT_DIR, '{sample}_variant_report.html'), sample=SAMPLES)
        )
    },
    'qc_reports': {
        'description': "Make QC reports.",
        'files': (
            expand(os.path.join(REPORT_DIR, '{sample}_qc_report.html'), sample=SAMPLES)
        )
    },
    'kraken_reports': {
        'description': "Make Kraken reports.",
        'files': (
            expand(os.path.join(REPORT_DIR, '{sample}_kraken_report.html'), sample=SAMPLES) +
            expand(os.path.join(REPORT_DIR, '{sample}_krona_report.html'), sample=SAMPLES)
        )
    }
}
selected_targets = ['variant_reports']
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
        r1 = os.path.join(TRIMMED_READS_DIR, "{sample}_1.fastq"),
        r2 = os.path.join(TRIMMED_READS_DIR, "{sample}_2.fastq")
    params:
        len_cutoff = int(150 * 0.8), # read length * pct_cutoff
        output = os.path.join(TRIMMED_READS_DIR, "{sample}")
    log: os.path.join(LOG_DIR, 'prinseq_{sample}.log')
    shell: "prinseq-lite.pl -fastq {input.r1} -fastq2 {input.r2} -ns_max_n 4 -min_qual_mean 30 -trim_qual_left 30 -trim_qual_right 30 -trim_qual_window 10 -out_good {params.output} -out_bad null -min_len {params.len_cutoff} >> {log} 2>&1"


rule bwa_index:
    input: REFERENCE_FASTA
    output: "{}.bwt".format(REFERENCE_FASTA)
    log: os.path.join(LOG_DIR, 'bwa_index.log')
    shell: "bwa index {input} >> {log} 2>&1"


rule bwa_align:
    input:
        fastq = [os.path.join(TRIMMED_READS_DIR, "{sample}_1.fastq"), os.path.join(TRIMMED_READS_DIR, "{sample}_2.fastq")],
        ref = REFERENCE_FASTA,
        index = "{}.bwt".format(REFERENCE_FASTA)
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.sam')
    params:
        threads = 4
    log: os.path.join(LOG_DIR, 'bwa_align_{sample}.log')
    shell: "bwa mem -t {params.threads} {input.ref} {input.fastq} > {output} 2>> {log} 3>&2"


rule samtools_filter_aligned:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.sam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned.bam')
    log: os.path.join(LOG_DIR, 'samtools_filter_aligned_{sample}.log')
    shell: "samtools view -bh -f 2 -F 2048 {input} > {output} 2>> {log} 3>&2"


# TODO: check command to get unaligned ones
rule samtools_filter_unaligned:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.sam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_unaligned.bam')
    log: os.path.join(LOG_DIR, 'samtools_filter_unaligned_{sample}.log')
    shell: "samtools view -bh -F 2 {input} > {output} 2>> {log} 3>&2"


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
    output: os.path.join(FASTQC_DIR, '{sample}_aligned_sorted_fastqc.zip')
    log: os.path.join(LOG_DIR, 'fastqc_{sample}.log')
    shell: "fastqc -o {FASTQC_DIR} -f bam {input} >> {log} 2>&1"


# TODO: check if --call-indels is needed?
rule lofreq:
    input:
        aligned_bam = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam'),
        aligned_bai = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bai'),
        ref = REFERENCE_FASTA
    output: os.path.join(VARIANTS_DIR, '{sample}_snv.vcf')
    log: os.path.join(LOG_DIR, 'lofreq_{sample}.log')
    shell: "lofreq call -f {input.ref} -o {output} --verbose {input.aligned_bam} >> {log} 2>&1"


rule vcf2csv:
    input: os.path.join(VARIANTS_DIR, '{sample}_snv.vcf')
    output: os.path.join(VARIANTS_DIR, '{sample}_snv.csv')
    params:
        script = os.path.join(SCRIPTS_DIR, 'vcfTocsv.py')
    log: os.path.join(LOG_DIR, 'vcf2csv_{sample}.log')
    shell: "python {params.script} {input} >> {log} 2>&1"


rule vep:
    input: os.path.join(VARIANTS_DIR, '{sample}_snv.vcf')
    output: os.path.join(VARIANTS_DIR, '{sample}_vep_sarscov2.txt')
    params:
        species = "sars_cov_2"
    log: os.path.join(LOG_DIR, 'vep_{sample}.log')
    shell: "vep --verbose --offline --dir_cache {VEP_DB} --DB_VERSION 101 --appris --biotype --buffer_size 5000 --check_existing --distance 5000 --mane --protein --species {params.species} --symbol --transcript_version --tsl --input_file {input} --output_file {output} >> {log} 2>&1"


rule parse_vep:
    input: os.path.join(VARIANTS_DIR, '{sample}_vep_sarscov2.txt')
    output: os.path.join(VARIANTS_DIR, '{sample}_vep_sarscov2_parsed.txt')
    params:
        script = os.path.join(SCRIPTS_DIR, 'parse_vep.py')
    log: os.path.join(LOG_DIR, 'parse_vep_{sample}.log')
    shell: "python {params.script} {VARIANTS_DIR} {input} {output} >> {log} 2>&1"


rule variant_report:
    input:
        vep_txt = os.path.join(VARIANTS_DIR, '{sample}_vep_sarscov2_parsed.txt'),
        snv_csv = os.path.join(VARIANTS_DIR, '{sample}_snv.csv')
    output: os.path.join(REPORT_DIR, '{sample}_variant_report.html')
    params:
        run_report = os.path.join(SCRIPTS_DIR, 'run_variant_report.R'),
        report_rmd = os.path.join(SCRIPTS_DIR, 'variant_report.Rmd')
    log: os.path.join(LOG_DIR, 'variant_report_{sample}.log')
    shell: "Rscript {params.run_report} --reportFile={params.report_rmd} --vep_txt_file={input.vep_txt} --snv_csv_file={input.snv_csv} --location_sigmuts={SIGMUT_DB} --sample_dir={VARIANTS_DIR} --sample_name={wildcards.sample} --outFile={output} >> {log} 2>&1"


rule bam2fastq:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_unaligned.bam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_unaligned.fastq')
    log: os.path.join(LOG_DIR, 'bam2fastq_{sample}.log')
    shell: "samtools fastq {input} > {output} 2>> {log} 3>&2"


rule kraken:
    input:
        unaligned_fastq = os.path.join(MAPPED_READS_DIR, '{sample}_unaligned.fastq'),
        database = KRAKEN_DB
    output: os.path.join(KRAKEN_DIR, '{sample}_classified_unaligned_reads.txt')
    log: os.path.join(LOG_DIR, 'kraken_{sample}.log')
    shell: "kraken2 --report {output} --db {input.database} {input.unaligned_fastq} >> {log} 2>&1"


rule kraken_report:
    input: os.path.join(KRAKEN_DIR, '{sample}_classified_unaligned_reads.txt')
    output: os.path.join(REPORT_DIR, '{sample}_kraken_report.html')
    params:
        run_report = os.path.join(SCRIPTS_DIR, 'run_kraken_report.R'),
        report_rmd = os.path.join(SCRIPTS_DIR, 'kraken_report.Rmd')
    log: os.path.join(LOG_DIR, 'kraken_report_{sample}.log')
    shell: "Rscript {params.run_report} --reportFile={params.report_rmd} --kraken_output={input} --sample_name={wildcards.sample} --output_file={output} >> {log} 2>&1"


rule krona_report:
    input: 
        kraken_output = os.path.join(KRAKEN_DIR, '{sample}_classified_unaligned_reads.txt'),
        database = KRONA_DB
    output: os.path.join(REPORT_DIR, '{sample}_krona_report.html')
    log: os.path.join(LOG_DIR, 'krona_report_{sample}.log')
    shell: "ktImportTaxonomy -m 3 -t 5 {input.kraken_output} -tax {input.database} -o {output} >> {log} 2>&1"


# TODO: does it have to be the 'unsorted' bam file? If so, we have to do indexing on those, too
rule samtools_bedcov:
    input:
        amplicons_bed = AMPLICONS_BED,
        aligned_bam = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam'),
        aligned_bai = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bai')
    output: os.path.join(COVERAGE_DIR, '{sample}_amplicons.csv')
    log: os.path.join(LOG_DIR, 'samtools_bedcov_{sample}.log')
    shell:"samtools bedcov {input.amplicons_bed} {input.aligned_bam} > {output} 2>> {log} 3>&2"


# TODO: does it have to be the 'unsorted' bam file? If so, we have to do indexing on those, too
rule samtools_coverage:
    input:
        aligned_bam = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam'),
        aligned_bai = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bai')
    output: os.path.join(COVERAGE_DIR, '{sample}_coverage.csv')
    log: os.path.join(LOG_DIR, 'samtools_coverage_{sample}.log')
    shell: "samtools coverage {input.aligned_bam} > {output} 2>> {log} 3>&2"


rule get_qc_table:
    input:
        coverage_csv = os.path.join(COVERAGE_DIR, '{sample}_coverage.csv'),
        amplicon_csv = os.path.join(COVERAGE_DIR, '{sample}_amplicons.csv')
    output: os.path.join(COVERAGE_DIR, '{sample}_merged_covs.csv')
    params:
        script = os.path.join(SCRIPTS_DIR, 'get_qc_table.py')
    log: os.path.join(LOG_DIR, 'get_qc_table_{sample}.log')
    shell: "python {params.script} {input.coverage_csv} {input.amplicon_csv} {output} >> {log} 2>&1"
        

# TODO: fix Error: unexpected end of input
# rule qc_report:
#     input: os.path.join(COVERAGE_DIR, '{sample}_merged_covs.csv')
#     output: os.path.join(REPORT_DIR, '{sample}_qc_report.html')
#     params:
#         run_report = os.path.join(SCRIPTS_DIR, 'run_qc_report.R'),
#         report_rmd = os.path.join(SCRIPTS_DIR, 'qc_report.Rmd')
#     log: os.path.join(LOG_DIR, 'qc_report_{sample}.log')
#     shell: "Rscript {params.run_report} --reportFile={params.report_rmd} --coverage_file={input} --sample_name={wildcards.sample} --outFile={output} >> {log} 2>&1"