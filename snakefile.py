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
import yaml
from itertools import chain
import re
import inspect

SAMPLE_SHEET_CSV = config['locations']['sample-sheet']
READS_DIR        = config['locations']['reads-dir']
REFERENCE_FASTA  = config['locations']['reference-fasta']
AMPLICONS_BED    = config['locations']['amplicons-bed']
KRAKEN_DB        = config['locations']['kraken-db-dir']
KRONA_DB         = config['locations']['krona-db-dir']
SIGMUT_DB        = config['locations']['sigmut-db-dir']
VEP_DB           = config['locations']['vep-db-dir']
OUTPUT_DIR       = config['locations']['output-dir']

# TODO: get default read length from multiqc
READ_LENGTH      = config['trimming']['read-length']
CUT_OFF          = config['trimming']['cut-off']

INDEX_DIR         = os.path.join(OUTPUT_DIR, 'index')
TRIMMED_READS_DIR = os.path.join(OUTPUT_DIR, 'trimmed_reads')
LOG_DIR           = os.path.join(OUTPUT_DIR, 'logs')
MAPPED_READS_DIR  = os.path.join(OUTPUT_DIR, 'mapped_reads')
VARIANTS_DIR      = os.path.join(OUTPUT_DIR, 'variants')
KRAKEN_DIR        = os.path.join(OUTPUT_DIR, 'kraken')
COVERAGE_DIR      = os.path.join(OUTPUT_DIR, 'coverage')
REPORT_DIR        = os.path.join(OUTPUT_DIR, 'report')
FASTQC_DIR        = os.path.join(REPORT_DIR, 'fastqc')
MULTIQC_DIR       = os.path.join(REPORT_DIR, 'multiqc')
SCRIPTS_DIR       = os.path.join(config['locations']['pkglibexecdir'], 'scripts/')
TMP_DIR           = os.path.join(config['locations']['output-dir'], 'pigx_work')

if os.getenv("PIGX_UNINSTALLED"):
    LOGO = os.path.join(config['locations']['pkgdatadir'], "images/Logo_PiGx.png")
else:
    LOGO = os.path.join(config['locations']['pkgdatadir'], "Logo_PiGx.png")


def toolArgs(name):
    if 'args' in config['tools'][name]:
        return config['tools'][name]['args']
    else:
        return ""

def tool(name):
    cmd = config['tools'][name]['executable']
    return cmd + " " + toolArgs(name)

BWA_EXEC             = tool("bwa")
FASTQC_EXEC          = tool("fastqc")
MULTIQC_EXEC         = tool("multiqc")
IMPORT_TAXONOMY_EXEC = tool("import_taxonomy")
KRAKEN2_EXEC         = tool("kraken2")
LOFREQ_EXEC          = tool("lofreq")
PRINSEQ_EXEC         = tool("prinseq")
PYTHON_EXEC          = tool("python")
RSCRIPT_EXEC         = tool("Rscript")
SAMTOOLS_EXEC        = tool("samtools")
VEP_EXEC             = tool("vep")

## Load sample sheet
with open(SAMPLE_SHEET_CSV, 'r') as fp:
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

targets = {
    'help': {
        'description': "Print all rules and their descriptions.",
        'files': []
    },
    'final_reports': {
        'description': "Produce a comprehensive report. This is the default target.",
        'files': (
            expand(os.path.join(REPORT_DIR, '{sample}.qc_report_per_sample.html'), sample=SAMPLES) +
            expand(os.path.join(REPORT_DIR, '{sample}.variantreport_p_sample.html'), sample=SAMPLES) +
            expand(os.path.join(REPORT_DIR, '{sample}.taxonomic_classification.html'), sample=SAMPLES) +
            expand(os.path.join(REPORT_DIR, '{sample}.Krona_report.html'), sample=SAMPLES) +
            [os.path.join(REPORT_DIR, 'index.html')]
        )
    },
    'lofreq': {
        'description': "Call variants and produce .vcf file and overview .csv file.",
        'files': (
            expand(os.path.join(VARIANTS_DIR, '{sample}_snv.csv'), sample=SAMPLES)
        )
    },
    'multiqc': {
        'description': "Create MultiQC reports for including raw and trimmed reads.",
        'files': (
            expand(os.path.join(MULTIQC_DIR, '{sample}', 'multiqc_report.html'), sample=SAMPLES)
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

# function to pass read files to trim/filter/qc improvement
def trim_reads_input(args):
  sample = args[0]
  return [os.path.join(READS_DIR, f) for f in lookup('name', sample, ['reads', 'reads2']) if f]

 # TODO: add extraction rule gunzip Output should be read.fastq
rule prinseq:
    input: trim_reads_input
    output:
        r1 = os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed_R1.fastq"),
        r2 = os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed_R2.fastq")
    params:
        len_cutoff = int(READ_LENGTH * CUT_OFF),
        output = os.path.join(TRIMMED_READS_DIR, "{sample}"),
        tmp_r1 = os.path.join(TRIMMED_READS_DIR, "{sample}_1.fastq"),
        tmp_r2 = os.path.join(TRIMMED_READS_DIR, "{sample}_2.fastq")
    log: os.path.join(LOG_DIR, 'prinseq_{sample}.log')
    shell: "{PRINSEQ_EXEC} -fastq {input[0]} -fastq2 {input[1]} -ns_max_n 4 -min_qual_mean 30 -trim_qual_left 30 -trim_qual_right 30 -trim_qual_window 10 -out_good {params.output} -out_bad null -min_len {params.len_cutoff} >> {log} 2>&1 && mv {params.tmp_r1} {output.r1} && mv {params.tmp_r2} {output.r2}"

# TODO add gzip command after prinseq, output of prinseq rule should be fastq.gz

rule bwa_index:
    input: REFERENCE_FASTA
    output:
      ref=os.path.join(INDEX_DIR, os.path.basename(REFERENCE_FASTA)),
      index=os.path.join(INDEX_DIR, "{}.bwt".format(os.path.basename(REFERENCE_FASTA)))
    log: os.path.join(LOG_DIR, 'bwa_index.log')
    shell: """mkdir -p {INDEX_DIR};
ln -sf {input} {INDEX_DIR};
cd {INDEX_DIR};
{BWA_EXEC} index {output.ref} >> {log} 2>&1"""


rule bwa_align:
    input:
        fastq = [os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed_R1.fastq"), os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed_R2.fastq")], # TODO: match to prinseq output
        ref = os.path.join(INDEX_DIR, "{}".format(os.path.basename(REFERENCE_FASTA))),
        index = os.path.join(INDEX_DIR, "{}.bwt".format(os.path.basename(REFERENCE_FASTA)))
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.sam')
    params:
        threads = 4
    log: os.path.join(LOG_DIR, 'bwa_align_{sample}.log')
    shell: "{BWA_EXEC} mem -t {params.threads} {input.ref} {input.fastq} > {output} 2>> {log} 3>&2"


rule samtools_filter_aligned:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.sam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned.bam')
    log: os.path.join(LOG_DIR, 'samtools_filter_aligned_{sample}.log')
    shell: "{SAMTOOLS_EXEC} view -bh -f 2 -F 2048 {input} > {output} 2>> {log} 3>&2"

rule samtools_filter_unaligned:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.sam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_unaligned.bam')
    log: os.path.join(LOG_DIR, 'samtools_filter_unaligned_{sample}.log')
    shell: "{SAMTOOLS_EXEC} view -bh -F 2 {input} > {output} 2>> {log} 3>&2"


rule samtools_sort:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned.bam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam')
    log: os.path.join(LOG_DIR, 'samtools_sort_{sample}.log')
    shell: "{SAMTOOLS_EXEC} sort -o {output} {input} >> {log} 2>&1"


rule samtools_index:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bai')
    log: os.path.join(LOG_DIR, 'samtools_index_{sample}.log')
    shell: "{SAMTOOLS_EXEC} index {input} {output} >> {log} 2>&1"


rule fastqc_raw:
    input: trim_reads_input
    output:
        r1_zip = os.path.join(FASTQC_DIR, '{sample}', '{sample}_R1_fastqc.zip'),
        r2_zip = os.path.join(FASTQC_DIR, '{sample}', '{sample}_R2_fastqc.zip')
    log: [os.path.join(LOG_DIR, 'fastqc_{sample}_raw_R1.log'), os.path.join(LOG_DIR, 'fastqc_{sample}_raw_R2.log')]
    params:
        output_dir = os.path.join(FASTQC_DIR, '{sample}')
    run:
        tmp_R1_output = os.path.join(FASTQC_DIR, '{sample}', os.path.basename({input[0]}).removesuffix('.fastq') + ".zip")
        tmp_R2_output = os.path.join(FASTQC_DIR, '{sample}', os.path.basename({input[1]}).removesuffix('.fastq') + ".zip")
        shell("{FASTQC_EXEC} -o {params.output_dir} {input} >> {log} 2>&1 && mv {tmp_R1_output} {output.r1_zip} && mv {tmp_R2_output} {output.r2_zip}")


rule fastqc_trimmed:
    input: os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed_R{read_num}.fastq")
    output:
        os.path.join(FASTQC_DIR, '{sample}', '{sample}_trimmed_R{read_num}_fastqc.zip')
    log: os.path.join(LOG_DIR, 'fastqc_{sample}_trimmed_R{read_num}.log')
    params:
        output_dir = os.path.join(FASTQC_DIR, '{sample}')
    shell: "{FASTQC_EXEC} -o {params.output_dir} {input} >> {log} 2>&1"


rule multiqc:
  input:
    fastqc_raw_output = expand(os.path.join(FASTQC_DIR, '{sample}', '{sample}_R{read_num}_fastqc.zip'), sample=SAMPLES, read_num=[1, 2]),
    fastqc_trimmed_output = expand(os.path.join(FASTQC_DIR, '{sample}', '{sample}_trimmed_R{read_num}_fastqc.zip'), sample=SAMPLES, read_num=[1, 2])
  output: os.path.join(MULTIQC_DIR, '{sample}', 'multiqc_report.html')
  params:
    fastqc_dir = os.path.join(FASTQC_DIR, '{sample}'),
    output_dir = os.path.join(MULTIQC_DIR, '{sample}')
  log: os.path.join(LOG_DIR, 'multiqc_{sample}.log')
  shell: "{MULTIQC_EXEC} -o {params.output_dir} {params.fastqc_dir} >> {log} 2>&1"


# TODO: check if --call-indels is needed?
rule lofreq:
    input:
        aligned_bam = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam'),
        aligned_bai = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bai'),
        ref = os.path.join(INDEX_DIR, "{}".format(os.path.basename(REFERENCE_FASTA)))
    output: os.path.join(VARIANTS_DIR, '{sample}_snv.vcf')
    log: os.path.join(LOG_DIR, 'lofreq_{sample}.log')
    shell: "{LOFREQ_EXEC} call -f {input.ref} -o {output} --verbose {input.aligned_bam} >> {log} 2>&1"


rule vcf2csv:
    input: os.path.join(VARIANTS_DIR, '{sample}_snv.vcf')
    output: os.path.join(VARIANTS_DIR, '{sample}_snv.csv')
    params:
        script = os.path.join(SCRIPTS_DIR, 'vcfTocsv.py')
    log: os.path.join(LOG_DIR, 'vcf2csv_{sample}.log')
    shell: "{PYTHON_EXEC} {params.script} {input} >> {log} 2>&1"

rule vep:
    input: os.path.join(VARIANTS_DIR, '{sample}_snv.vcf')
    output: os.path.join(VARIANTS_DIR, '{sample}_vep_sarscov2.txt')
    params:
        species = "sars_cov_2"
    log: os.path.join(LOG_DIR, 'vep_{sample}.log')
    shell:
      """
      {VEP_EXEC} --verbose --offline --dir_cache {VEP_DB} --DB_VERSION 101 --appris --biotype --buffer_size 5000 --check_existing\
      --distance 5000 --mane --protein --species {params.species} --symbol --transcript_version --tsl\
      --input_file {input} --output_file {output} >> {log} 2>&1
      """

rule parse_vep:
    input: os.path.join(VARIANTS_DIR, '{sample}_vep_sarscov2.txt')
    output: os.path.join(VARIANTS_DIR, '{sample}_vep_sarscov2_parsed.txt')
    params:
        script = os.path.join(SCRIPTS_DIR, 'parse_vep.py')
    log: os.path.join(LOG_DIR, 'parse_vep_{sample}.log')
    shell: "{PYTHON_EXEC} {params.script} {input} {output} >> {log} 2>&1"


rule bam2fastq:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_unaligned.bam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_unaligned.fastq')
    log: os.path.join(LOG_DIR, 'bam2fastq_{sample}.log')
    shell: "{SAMTOOLS_EXEC} fastq {input} > {output} 2>> {log} 3>&2"


rule kraken:
    input:
        unaligned_fastq = os.path.join(MAPPED_READS_DIR, '{sample}_unaligned.fastq'),
        database = KRAKEN_DB
    output: os.path.join(KRAKEN_DIR, '{sample}_classified_unaligned_reads.txt')
    log: os.path.join(LOG_DIR, 'kraken_{sample}.log')
    shell: "{KRAKEN2_EXEC} --report {output} --db {input.database} {input.unaligned_fastq} >> {log} 2>&1"


rule krona_report:
    input:
        kraken_output = os.path.join(KRAKEN_DIR, '{sample}_classified_unaligned_reads.txt'),
        database = KRONA_DB
    output: os.path.join(REPORT_DIR, '{sample}.Krona_report.html')
    log: os.path.join(LOG_DIR, 'krona_report_{sample}.log')
    shell: "{IMPORT_TAXONOMY_EXEC} -m 3 -t 5 {input.kraken_output} -tax {input.database} -o {output} >> {log} 2>&1"


rule samtools_bedcov:
    input:
        amplicons_bed = AMPLICONS_BED,
        aligned_bam = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam'),
        aligned_bai = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bai')
    output: os.path.join(COVERAGE_DIR, '{sample}_amplicons.csv')
    log: os.path.join(LOG_DIR, 'samtools_bedcov_{sample}.log')
    shell: "{SAMTOOLS_EXEC} bedcov {input.amplicons_bed} {input.aligned_bam} > {output} 2>> {log} 3>&2"


rule samtools_coverage:
    input:
        aligned_bam = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam'),
        aligned_bai = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bai')
    output: os.path.join(COVERAGE_DIR, '{sample}_coverage.csv')
    log: os.path.join(LOG_DIR, 'samtools_coverage_{sample}.log')
    shell: "{SAMTOOLS_EXEC} coverage {input.aligned_bam} > {output} 2>> {log} 3>&2"


rule get_qc_table:
    input:
        coverage_csv = os.path.join(COVERAGE_DIR, '{sample}_coverage.csv'),
        amplicon_csv = os.path.join(COVERAGE_DIR, '{sample}_amplicons.csv')
    output: os.path.join(COVERAGE_DIR, '{sample}_merged_covs.csv')
    params:
        script = os.path.join(SCRIPTS_DIR, 'get_qc_table.py')
    log: os.path.join(LOG_DIR, 'get_qc_table_{sample}.log')
    shell: "{PYTHON_EXEC} {params.script} {input.coverage_csv} {input.amplicon_csv} {output} >> {log} 2>&1"


rule generate_navbar:
    input:
      script = os.path.join(SCRIPTS_DIR, "generateNavigation.R")
    output:
      os.path.join(REPORT_DIR, "_navbar.html")
    params:
      report_scripts_dir = os.path.join(SCRIPTS_DIR, "report_scripts")
    log: os.path.join(LOG_DIR, "generate_navigation.log")
    shell: "{RSCRIPT_EXEC} {input.script} \
{params.report_scripts_dir} {SAMPLE_SHEET_CSV} {output} > {log} 2>&1"


rule render_kraken2_report:
    input:
      script=os.path.join(SCRIPTS_DIR, "renderReport.R"),
      report=os.path.join(SCRIPTS_DIR, "report_scripts", "taxonomic_classification.Rmd"),
      header=os.path.join(REPORT_DIR, "_navbar.html"),
      kraken=os.path.join(KRAKEN_DIR, "{sample}_classified_unaligned_reads.txt"),
      krona=os.path.join(REPORT_DIR, "{sample}.Krona_report.html")
    output: os.path.join(REPORT_DIR, "{sample}.taxonomic_classification.html")
    log: os.path.join(LOG_DIR, "reports", "{sample}_taxonomic_classification.log")
    shell: """{RSCRIPT_EXEC} {input.script} \
{input.report} {output} {input.header} \
'{{\
  "sample_name": "{wildcards.sample}",  \
  "site_dir":    "{REPORT_DIR}",        \
  "krona_file":  "{input.krona}",       \
  "kraken_file": "{input.kraken}",      \
  "logo": "{LOGO}" \
}}' > {log} 2>&1"""


rule render_variant_report:
    input:
      script=os.path.join(SCRIPTS_DIR, "renderReport.R"),
      report=os.path.join(SCRIPTS_DIR,"report_scripts", "variantreport_p_sample.Rmd"),
      header=os.path.join(REPORT_DIR, "_navbar.html"),
      vep=os.path.join(VARIANTS_DIR, "{sample}_vep_sarscov2_parsed.txt"),
      snv=os.path.join(VARIANTS_DIR, "{sample}_snv.csv")
    output: os.path.join(REPORT_DIR, "{sample}.variantreport_p_sample.html")
    log: os.path.join(LOG_DIR, "reports", "{sample}_variant_report.log")
    shell: """{RSCRIPT_EXEC} {input.script} \
{input.report} {output} {input.header} \
'{{\
  "sample_name":  "{wildcards.sample}",  \
  "sigmut_db":    "{SIGMUT_DB}",         \
  "variants_dir": "{VARIANTS_DIR}",      \
  "vep_file":     "{input.vep}",         \
  "snv_file":     "{input.snv}",         \
  "sample_sheet": "{SAMPLE_SHEET_CSV}",  \
  "logo": "{LOGO}" \
}}' > {log} 2>&1"""


rule render_qc_report:
    input:
      script=os.path.join(SCRIPTS_DIR, "renderReport.R"),
      report=os.path.join(SCRIPTS_DIR, "report_scripts", "qc_report_per_sample.Rmd"),
      header=os.path.join(REPORT_DIR, "_navbar.html"),
      coverage=os.path.join(COVERAGE_DIR, "{sample}_merged_covs.csv"),
      multiqc=os.path.join(MULTIQC_DIR, '{sample}', 'multiqc_report.html')
    output:
      os.path.join(REPORT_DIR, "{sample}.qc_report_per_sample.html")
    params:
      multiqc_rel_path=lambda wildcards, input: input.multiqc[len(REPORT_DIR)+1:]
    log: os.path.join(LOG_DIR, "reports", "{sample}_qc_report.log")
    shell: """{RSCRIPT_EXEC} {input.script} \
{input.report} {output} {input.header} \
'{{\
  "sample_name": "{wildcards.sample}",  \
  "coverage_file": "{input.coverage}",   \
  "multiqc_report": "{params.multiqc_rel_path}", \
  "logo": "{LOGO}" \
}}' > {log} 2>&1"""


rule render_index:
    input:
      script=os.path.join(SCRIPTS_DIR, "renderReport.R"),
      report=os.path.join(SCRIPTS_DIR, "report_scripts", "index.Rmd"),
      header=os.path.join(REPORT_DIR, "_navbar.html"),
      # TODO: see comment below
      side_effects=expand(os.path.join(REPORT_DIR, "{sample}.variantreport_p_sample.html"), sample = SAMPLES),
      # This can only be done after all other reports have been built,
      # because these reports are referenced in the overview section.
      taxonomy=expand(os.path.join(REPORT_DIR, "{sample}.taxonomic_classification.html"), sample = SAMPLES),
      krona=expand(os.path.join(REPORT_DIR, "{sample}.Krona_report.html"), sample = SAMPLES),
      qc=expand(os.path.join(REPORT_DIR, "{sample}.qc_report_per_sample.html"), sample = SAMPLES),
      variant=expand(os.path.join(REPORT_DIR, "{sample}.variantreport_p_sample.html"), sample = SAMPLES),
    # TODO: these CSV files should be declared as inputs!  Due to
    # https://github.com/BIMSBbioinfo/pigx_sarscov2_ww/issues/19 we
    # cannot do this yet, so we just add the variant reports for all
    # samples as inputs.
    params:
      variants=os.path.join(VARIANTS_DIR, 'data_variant_plot.csv'),
      mutations=os.path.join(VARIANTS_DIR, 'data_mutation_plot.csv')
    output: os.path.join(REPORT_DIR, "index.html")
    log: os.path.join(LOG_DIR, "reports", "index.log")
    shell: """{RSCRIPT_EXEC} {input.script} \
{input.report} {output} {input.header}   \
'{{                                      \
  "variants_csv": "{params.variants}",   \
  "mutations_csv": "{params.mutations}", \
  "sample_sheet": "{SAMPLE_SHEET_CSV}",  \
  "logo": "{LOGO}" \
}}' > {log} 2>&1"""
