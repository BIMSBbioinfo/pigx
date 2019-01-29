"""
Snakefile for pigx crispr pipeline
"""

import os
import yaml
import csv
import inspect


# tools
RSCRIPT = config['tools']['Rscript']
JAVA = config['tools']['java']['path']
JAVA_MEM = config['tools']['java']['mem']
JAVA_THREADS = config['tools']['java']['threads']
GATK = config['tools']['gatk']

# input locations
SRC_DIR = config['source-dir']
READS_DIR = config['reads-dir']
ADAPTERS = config['adapters']
SAMPLE_SHEET_FILE = config['sample_sheet']
CUT_SITES_FILE = config['cutsites']
COMPARISONS_FILE = config.get('comparisonsFile', {})
REFERENCE_FASTA = config['reference_fasta']

#output locations
OUTPUT_DIR = config['output-dir']
TRIMMED_READS_DIR = os.path.join(OUTPUT_DIR, 'trimmed_reads')
LOG_DIR           = os.path.join(OUTPUT_DIR, 'logs')
FASTA_DIR         = os.path.join(OUTPUT_DIR, 'fasta')
FASTQC_DIR        = os.path.join(OUTPUT_DIR, 'fastqc')
MULTIQC_DIR       = os.path.join(OUTPUT_DIR, 'multiqc')
MAPPED_READS_DIR  = os.path.join(OUTPUT_DIR, 'aln')
INDELS_DIR      = os.path.join(OUTPUT_DIR, 'indels')
BBMAP_INDEX_DIR   = os.path.join(OUTPUT_DIR, 'bbmap_indexes')
REPORT_DIR        = os.path.join(OUTPUT_DIR, 'reports')
GATK_DIR          = os.path.join(OUTPUT_DIR, 'gatk')

#other parameters
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
  #print("sample:",sample)
  return [os.path.join(READS_DIR, f) for f in lookup('sample_name', sample, ['reads']) if f]

def get_output_file_list(DIR, ext):
    paths = list()
    for sample in SAMPLES:
        amplicon = lookup('sample_name', sample, ['amplicon'])[0]
        paths.append(os.path.join(DIR, amplicon, ".".join([sample, ext])))
    return(paths)

def get_bbmap_command(wc):
    sample = wc.sample
    tech = lookup('sample_name', sample, ['tech'])[0]
    if tech == 'illumina':
        return('bbmap.sh')
    elif tech == 'pacbio':
        return('mapPacBio.sh maxlen=6000')

## Check which amplicons are needed for comparisons are to be made
with open(COMPARISONS_FILE, 'r') as fp:
  rows =  [row for row in csv.reader(fp, delimiter='\t')]
  header = rows[0]; rows = rows[1:]
  COMPARISONS = [dict(zip(header, row)) for row in rows]
  COMPARISON_AMPLICONS = [COMPARISONS[i]['amplicon'] for i in range(len(COMPARISONS))]


rule all:
    input:
        #os.path.join(FASTA_DIR, os.path.basename(REFERENCE_FASTA)),
        #expand(os.path.join(FASTQC_DIR, "{sample}.fastqc.done"), sample = SAMPLES),
        expand(os.path.join(MAPPED_READS_DIR, "{sample}.bam.bai"), sample = SAMPLES),
        #expand(os.path.join(GATK_DIR, "{sample}.realigned.bam"), sample = SAMPLES),
        os.path.join(OUTPUT_DIR, "multiqc", "multiqc_report.html"),
        expand(os.path.join(INDELS_DIR, "{sample}", "{sample}.indelScores.bigwig"), sample = SAMPLES),
        expand(os.path.join(INDELS_DIR, "{sample}", "{sample}.deletionScores.bigwig"), sample = SAMPLES),
        expand(os.path.join(INDELS_DIR, "{sample}", "{sample}.insertionScores.bigwig"), sample = SAMPLES),
        expand(os.path.join(INDELS_DIR, "{sample}", "{sample}.alnCoverage.bigwig"), sample = SAMPLES),
        expand(os.path.join(INDELS_DIR, "{sample}", "{sample}.deletions.bed"), sample = SAMPLES),
        expand(os.path.join(INDELS_DIR, "{sample}", "{sample}.insertions.bed"), sample = SAMPLES),
        expand(os.path.join(INDELS_DIR, "{sample}", "{sample}.indels.tsv"), sample = SAMPLES),
        expand(os.path.join(INDELS_DIR, "{sample}", "{sample}.reads_with_indels.tsv"), sample = SAMPLES),
        expand(os.path.join(INDELS_DIR, "{sample}", "{sample}.insertedSequences.tsv"), sample = SAMPLES),
        #get_output_file_list(INDELS_DIR, "freeBayes_variants.vcf"),
        #get_output_file_list(INDELS_DIR, "freeBayes_deletions.bed"),
        #get_output_file_list(INDELS_DIR, "freeBayes_insertions.bed"),
        #expand(os.path.join(REPORT_DIR, "{amplicon}.report.html"), amplicon=AMPLICONS.keys()),
        #expand(os.path.join(REPORT_DIR, 'comparisons', '{amplicon}.report.comparisons.html'), amplicon=COMPARISON_AMPLICONS),
        #expand(os.path.join(REPORT_DIR, 'comparisons', '{amplicon}.comparison.stats.tsv'), amplicon=COMPARISON_AMPLICONS)

#notice that get_amplicon_file function for 'fasta' should only be used once.
#Other rules that need the amplicon fasta sequence as input should use :
#lambda wildcards: os.path.join(FASTA_DIR, ''.join([wildcards.amplicon, ".fasta"]))
rule reformatFasta:
    input: REFERENCE_FASTA
    output: os.path.join(FASTA_DIR, os.path.basename(REFERENCE_FASTA))
    log: os.path.join(LOG_DIR, ".".join(["reformatFasta", os.path.basename(REFERENCE_FASTA), "log"]))
    shell:
        "reformat.sh in={input} out={output} tuc > {log} 2>&1"

rule getFastaIndex:
    input: os.path.join(FASTA_DIR, os.path.basename(REFERENCE_FASTA))
    output: os.path.join(FASTA_DIR, ".".join([os.path.basename(REFERENCE_FASTA), 'fai']))
    log: os.path.join(LOG_DIR, ".".join(["getFastaIndex", os.path.basename(REFERENCE_FASTA), "log"]))
    shell:
        "samtools faidx {input} > {log} 2>&1"

rule getFastaDict:
    input: os.path.join(FASTA_DIR, os.path.basename(REFERENCE_FASTA))
    output: os.path.join(FASTA_DIR, ".".join([os.path.basename(REFERENCE_FASTA), 'dict']))
    params:
        script=os.path.join(SRC_DIR, "src", "getFastaDict.R")
    log: os.path.join(LOG_DIR, ".".join(["getFastaDict", os.path.basename(REFERENCE_FASTA), "log"]))
    shell:
        "{RSCRIPT} {params.script} {input} {output} > {log} 2>&1"

rule fastqc:
    input: reads_input
    output: os.path.join(FASTQC_DIR, "{sample}.fastqc.done")
    log: os.path.join(LOG_DIR, "FASTQC", "{sample}.fastqc.log")
    shell: "fastqc -o {FASTQC_DIR} {input} > {log} 2>&1; touch {output}"

rule trimmomatic:
    input: reads_input
    output: os.path.join(TRIMMED_READS_DIR, "{sample}.fastq.gz")
    log: os.path.join(LOG_DIR, "TRIM", "trimmomatic.{sample}.log")
    params:
        tech = lambda wildcards: lookup('sample_name', wildcards.sample, ['tech'])[0]
    run:
        if params.tech == 'illumina':
            shell("trimmomatic SE -threads 2 {input} {output} \
            ILLUMINACLIP:{ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > {log} 2>&1")
        elif params.tech == 'pacbio':
            shell("ln -s {input} {output}")

rule bbmap_indexgenome:
    input: os.path.join(FASTA_DIR, os.path.basename(REFERENCE_FASTA))
    output: os.path.join(BBMAP_INDEX_DIR, re.sub("\.fa(sta)?$", "", os.path.basename(REFERENCE_FASTA)))
    log: os.path.join(LOG_DIR, ".".join(["bbmap_index", os.path.basename(REFERENCE_FASTA), "log"]))
    shell: "bbmap.sh ref={input} path={output} > {log} 2>&1"

rule bbmap_map:
    input:
        reads = os.path.join(TRIMMED_READS_DIR, '{sample}.fastq.gz'),
        ref = os.path.join(BBMAP_INDEX_DIR, re.sub("\.fa(sta)?$", "", os.path.basename(REFERENCE_FASTA)))
    output:
        os.path.join(MAPPED_READS_DIR, "{sample}.sam")
    params:
        aligner = lambda wildcards: get_bbmap_command(wildcards)
    log: os.path.join(LOG_DIR, "BBMAP", "bbmap_align.{sample}.log")
    shell:
        "{params.aligner} path={input.ref} in={input.reads} outm={output} t=2 > {log} 2>&1"

# GATK requires read groups, so here we add some dummy read group information to the bam files
rule samtools_addReadGroups:
    input: os.path.join(MAPPED_READS_DIR, "{sample}.sam")
    output: os.path.join(MAPPED_READS_DIR, "{sample}.sam_withreadgroups")
    log: os.path.join(LOG_DIR, "SAMTOOLS", "samtools_addReadGroups.{sample}.log")
    shell:
        """
        samtools addreplacerg -m overwrite_all -o {output} -r ID:pigx_crispr -r PL:illumina -r SM:{wildcards.sample} {input} > {log} 2>&1
        rm {input}
        """

rule samtools_sam2bam:
    input: os.path.join(MAPPED_READS_DIR, "{sample}.sam_withreadgroups")
    output: os.path.join(MAPPED_READS_DIR, "{sample}.bam")
    log: os.path.join(LOG_DIR, "SAMTOOLS", "samtools_sam2bam.{sample}.log")
    shell:
        """
        samtools view -bh {input} | samtools sort -o {output} > {log} 2>&1
        rm {input}
        """

rule samtools_indexbam:
    input: os.path.join(MAPPED_READS_DIR, "{sample}.bam")
    output: os.path.join(MAPPED_READS_DIR, "{sample}.bam.bai")
    log: os.path.join(LOG_DIR, "samtools_indexbam.{sample}.log")
    shell: "samtools index {input} > {log} 2>&1"

# we need an interval file that tells gatk to correct indels within those intervals
# we use the whole amplicon start end positions for this purpose, however, in genome-wide
# context, it should be limited to certain regions.
rule get_gatk_realigner_intervals:
    input:
        ref = os.path.join(FASTA_DIR, os.path.basename(REFERENCE_FASTA))
    output:
        os.path.join(GATK_DIR, "_".join([os.path.basename(REFERENCE_FASTA), "realigner.intervals"]))
    params:
        script=os.path.join(SRC_DIR, "src", "get_gatk_realigner_intervals.R")
    shell:
        #find the length of the amplicon sequence and print the start-end positions in a file
        "{RSCRIPT} {params.script} {input.ref} {output}"

rule gatk_indelRealigner:
    input:
        ref = os.path.join(FASTA_DIR, os.path.basename(REFERENCE_FASTA)),
        ref_index = os.path.join(FASTA_DIR, ".".join([os.path.basename(REFERENCE_FASTA), 'fai'])),
        ref_dict = os.path.join(FASTA_DIR, ".".join([os.path.basename(REFERENCE_FASTA), 'dict'])),
        bamIndex = os.path.join(MAPPED_READS_DIR, "{sample}.bam.bai"),
        bamFile = os.path.join(MAPPED_READS_DIR, "{sample}.bam"),
        intervals = os.path.join(GATK_DIR, "_".join([os.path.basename(REFERENCE_FASTA), "realigner.intervals"]))
    output:
        bam = os.path.join(GATK_DIR, "{sample}.realigned.bam"),
        bai = os.path.join(GATK_DIR, "{sample}.realigned.bai")
    log: os.path.join(LOG_DIR, "GATK", "gatk_realigner.{sample}.log")
    shell: "{JAVA} {JAVA_MEM} {JAVA_THREADS} -jar {GATK} -T IndelRealigner -nt 1 -R {input.ref} -I {input.bamFile} -targetIntervals {input.intervals} -o {output.bam} > {log} 2>&1"


rule samtools_stats:
    input:
        bamfile = os.path.join(MAPPED_READS_DIR, "{sample}.bam"),
        ref = os.path.join(FASTA_DIR, os.path.basename(REFERENCE_FASTA))
    output: os.path.join(MAPPED_READS_DIR, "SAMTOOLS", "{sample}.samtools.stats.txt")
    log: os.path.join(LOG_DIR, "SAMTOOLS", "samtools_stats.{sample}.log")
    shell: "samtools stats --reference {input.ref} {input.bamfile} > {output} 2> {log}"

rule multiqc:
    input:
        fastqc = expand(os.path.join(FASTQC_DIR, "{sample}.fastqc.done"), sample = SAMPLES),
        trimmomatic = expand(os.path.join(TRIMMED_READS_DIR, "{sample}.fastq.gz"), sample = SAMPLES),
        samtools = expand(os.path.join(MAPPED_READS_DIR, "SAMTOOLS", "{sample}.samtools.stats.txt"), sample = SAMPLES)
    output:
        os.path.join(OUTPUT_DIR, "multiqc", "multiqc_report.html")
    params:
        analysis_folder = OUTPUT_DIR,
        output_folder = os.path.join(OUTPUT_DIR, "multiqc")
    log: os.path.join(LOG_DIR, 'multiqc.log')
    shell: "multiqc --force -o {params.output_folder} {params.analysis_folder} > {log} 2>&1"

# runs freebayes tool on top of a BAM file and creates a VCF file containing variants
# we assume the sequences may come from multiple individuals and suppress snps.
# we apply a frequency cut-off of 0.1%
# rule getFreeBayesVariants:
#     input:
#         ref = lambda wildcards:  os.path.join(FASTA_DIR, ''.join([wildcards.amplicon, ".fasta"])),
#         bamIndex = os.path.join(GATK_DIR, "{amplicon}", "{sample}.realigned.bai"),
#         bamFile = os.path.join(GATK_DIR, "{amplicon}", "{sample}.realigned.bam")
#     output:
#         os.path.join(INDELS_DIR, "{amplicon}", "{sample}.freeBayes_variants.vcf")
#     log: os.path.join(LOG_DIR, "{amplicon}", "getFreeBayesVariants.{sample}.log")
#     shell: "freebayes -f {input.ref} -F 0.001 -C 1 --no-snps --use-duplicate-reads --pooled-continuous {input.bamFile} > {output} 2> {log}"
#
# #convert VCF output from freebayes to BED files
# rule vcf2bed:
#     input: os.path.join(INDELS_DIR, "{amplicon}", "{sample}.freeBayes_variants.vcf")
#     output:
#         deletions = os.path.join(INDELS_DIR, "{amplicon}", "{sample}.freeBayes_deletions.bed"),
#         insertions = os.path.join(INDELS_DIR, "{amplicon}", "{sample}.freeBayes_insertions.bed")
#     log:
#         deletions = os.path.join(LOG_DIR, "{amplicon}", "vcf2bed.deletions.{sample}.log"),
#         insertions = os.path.join(LOG_DIR, "{amplicon}", "vcf2bed.insertions.{sample}.log")
#     shell:
#         """
#         cut -f 1,2,3,4,5,6,7,8 {input} | vcf2bed --deletions > {output.deletions} 2> {log.deletions}
#         cut -f 1,2,3,4,5,6,7,8 {input} | vcf2bed --insertions > {output.insertions} 2> {log.insertions}
#         """

rule getIndelStats:
    input:
        #bamIndex = os.path.join(GATK_DIR, "{amplicon}", "{sample}.realigned.bai"),
        #bamFile = os.path.join(GATK_DIR, "{amplicon}", "{sample}.realigned.bam"),
        bamIndex = os.path.join(MAPPED_READS_DIR,  "{sample}.bam.bai"),
        bamFile = os.path.join(MAPPED_READS_DIR,  "{sample}.bam")
    output:
        os.path.join(INDELS_DIR, "{sample}", "{sample}.indelScores.bigwig"),
        os.path.join(INDELS_DIR, "{sample}", "{sample}.deletionScores.bigwig"),
        os.path.join(INDELS_DIR, "{sample}", "{sample}.insertionScores.bigwig"),
        os.path.join(INDELS_DIR, "{sample}", "{sample}.alnCoverage.bigwig"),
        # os.path.join(INDELS_DIR, "{amplicon}", "{sample}.coverageStats.tsv"),
        # os.path.join(INDELS_DIR, "{amplicon}", "{sample}.indel_stats_at_cutsites.tsv"),
        os.path.join(INDELS_DIR, "{sample}", "{sample}.deletions.bed"),
        os.path.join(INDELS_DIR, "{sample}", "{sample}.insertions.bed"),
        os.path.join(INDELS_DIR, "{sample}", "{sample}.indels.tsv"),
        os.path.join(INDELS_DIR, "{sample}", "{sample}.reads_with_indels.tsv"),
        os.path.join(INDELS_DIR, "{sample}", "{sample}.insertedSequences.tsv")
    params:
        sgRNA_list = lambda wildcards: lookup('sample_name', wildcards.sample, ['sgRNA_ids'])[0],
        script=os.path.join(SRC_DIR, "src", "getIndelStats.R")
    log: os.path.join(LOG_DIR, "indel_stats", "getIndelStats_{sample}.log")
    shell: "{RSCRIPT} {params.script} {input.bamFile} {wildcards.sample} {INDELS_DIR} {CUT_SITES_FILE} {params.sgRNA_list} > {log} 2>&1"

# rule report:
#     input:
#         coverageStats = get_output_file_list(INDELS_DIR, "coverageStats.tsv"),
#         cutsiteStats = get_output_file_list(INDELS_DIR, "indel_stats_at_cutsites.tsv"),
#         indels = get_output_file_list(INDELS_DIR, "indels.unfiltered.tsv")
#     params:
#       fasta = lambda wildcards:  os.path.join(FASTA_DIR, ''.join([wildcards.amplicon, ".fasta"])),
#       cutSitesFile = lambda wildcards: get_amplicon_file(wildcards, 'cutsites'),
#       reportR = os.path.join(SRC_DIR, "src", "runReport.R"),
#       reportRmd = os.path.join(SRC_DIR, "src", "report.Rmd"),
#       indelsFolder = os.path.join(INDELS_DIR, "{amplicon}")
#     log: os.path.join(LOG_DIR, "{amplicon}.report.log")
#     output:
#         os.path.join(REPORT_DIR, '{amplicon}.report.html')
#     shell:
#         "{RSCRIPT} {params.reportR}  --reportFile={params.reportRmd} --ampliconFastaFile={params.fasta} --ampliconName={wildcards.amplicon} --cutSitesFile={params.cutSitesFile} --sampleSheetFile={SAMPLE_SHEET_FILE} --indelsFolder={params.indelsFolder} --workdir={REPORT_DIR} --prefix={wildcards.amplicon} > {log} 2>&1"
#
#
# rule report_comparisons:
#     input:
#         coverageStats = get_output_file_list(INDELS_DIR, "coverageStats.tsv"),
#         cutsiteStats = get_output_file_list(INDELS_DIR, "indel_stats_at_cutsites.tsv"),
#         indels = get_output_file_list(INDELS_DIR, "indels.unfiltered.tsv"),
#         comp = COMPARISONS_FILE
#     params:
#         reportR = os.path.join(SRC_DIR, "src", "runReport.comparisons.R"),
#         reportRmd = os.path.join(SRC_DIR, "src", "report.comparisons.Rmd"),
#         indelsFolder = os.path.join(INDELS_DIR, "{amplicon}"),
#         outDir = os.path.join(REPORT_DIR, 'comparisons'),
#     log: os.path.join(LOG_DIR, "{amplicon}.report.comparisons.log")
#     output:
#         os.path.join(REPORT_DIR, 'comparisons', '{amplicon}.report.comparisons.html'),
#         os.path.join(REPORT_DIR, 'comparisons', '{amplicon}.comparison.stats.tsv')
#     shell:
#         "{RSCRIPT} {params.reportR}  --reportFile={params.reportRmd} --ampliconName={wildcards.amplicon} --comparisonsFile={COMPARISONS_FILE} --indelsFolder={params.indelsFolder} --workdir={params.outDir} --prefix={wildcards.amplicon} > {log} 2>&1"
