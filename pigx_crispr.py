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
TARGET_NAMES = [line['target_name'] for line in SAMPLE_SHEET]

def reads_input(wc):
  sample = wc.sample
  #print("sample:",sample)
  return [os.path.join(READS_DIR, f) for f in lookup('sample_name', sample, ['reads']) if f]

def get_bbmap_command(wc):
    sample = wc.sample
    tech = lookup('sample_name', sample, ['tech'])[0]
    if tech == 'illumina':
        return('bbmap.sh')
    elif tech == 'pacbio':
        return('mapPacBio.sh maxlen=6000')

rule all:
    input:
        #expand(os.path.join(FASTQC_DIR, "{sample}.fastqc.done"), sample = SAMPLES),
        #expand(os.path.join(MAPPED_READS_DIR, "{sample}.bam.bai"), sample = SAMPLES),
        #expand(os.path.join(GATK_DIR, "{sample}.indels.realigned.bam"), sample = SAMPLES),
        expand(os.path.join(OUTPUT_DIR,  "aln_merged", "{sample}.bam"), sample = SAMPLES),
        #expand(os.path.join(OUTPUT_DIR, "SAMTOOLS", "{sample}.samtools.stats.txt"), sample = SAMPLES),
        os.path.join(OUTPUT_DIR, "multiqc", "multiqc_report.html"),
        #expand(os.path.join(INDELS_DIR, "{sample}", "{sample}.sgRNA_efficiency.tsv"), sample = SAMPLES),
        os.path.join(REPORT_DIR, "index.html")
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
    output: os.path.join(FASTA_DIR, ".".join([os.path.basename(REFERENCE_FASTA).replace(".fa", ""), 'dict']))
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
        samtools view -bh {input} | samtools rmdup -s - - | samtools sort -o {output} > {log} 2>&1
        rm {input}
        """

rule samtools_indexbam:
    input: os.path.join(MAPPED_READS_DIR, "{sample}.bam")
    output: os.path.join(MAPPED_READS_DIR, "{sample}.bam.bai")
    log: os.path.join(LOG_DIR, "samtools_indexbam.{sample}.log")
    shell: "samtools index {input} > {log} 2>&1"


# split the bam file into two: First contains only reads alignments with indels
# second contains alignments without indels
rule split_bam_by_indels:
    input:
        bamFile = os.path.join(MAPPED_READS_DIR, "{sample}.bam"),
        bamIndex = os.path.join(MAPPED_READS_DIR, "{sample}.bam.bai")
    output:
        os.path.join(OUTPUT_DIR, "aln_split_by_indels", "{sample}.with_indels.bam"),
        os.path.join(OUTPUT_DIR, "aln_split_by_indels", "{sample}.with_indels.bam.bai"),
        os.path.join(OUTPUT_DIR, "aln_split_by_indels", "{sample}.without_indels.bam"),
        os.path.join(OUTPUT_DIR, "aln_split_by_indels", "{sample}.without_indels.bam.bai"),
    params:
        script = os.path.join(SRC_DIR, "src", "split_bam_by_indels.R"),
        tech = lambda wildcards: lookup('sample_name', wildcards.sample, ['tech'])[0],
        output_dir = os.path.join(OUTPUT_DIR, "aln_split_by_indels")
    log: os.path.join(LOG_DIR, "split_bam_by_indels.{sample}.log")
    shell: "{RSCRIPT} {params.script} {input.bamFile} {wildcards.sample} {params.output_dir} > {log} 2>&1"


# we need an interval file that tells gatk to correct indels within those intervals
# we provide the target_region that is available
rule get_gatk_realigner_intervals:
    input:
        ref = os.path.join(FASTA_DIR, os.path.basename(REFERENCE_FASTA)),
        ref_index = os.path.join(FASTA_DIR, ".".join([os.path.basename(REFERENCE_FASTA), 'fai'])),
        ref_dict = os.path.join(FASTA_DIR, ".".join([os.path.basename(REFERENCE_FASTA).replace(".fa", ""), 'dict'])),
        bamIndex = os.path.join(OUTPUT_DIR, "aln_split_by_indels", "{sample}.with_indels.bam.bai"),
        bamFile = os.path.join(OUTPUT_DIR, "aln_split_by_indels", "{sample}.with_indels.bam")
    output:
        os.path.join(GATK_DIR, "{sample}.gatk_realigner.intervals")
    params:
        intervals = lambda wildcards: lookup('sample_name', wildcards.sample, ['target_region'])[0]
    log: os.path.join(LOG_DIR, "GATK", "gatk_realigner_target_creator.{sample}.log")
    #shell:
        #syntax for intervals: <chromosome>:start-end
        #(here we remove the strand info that may exist in target_region)
        #"echo {params.intervals} | sed 's/:[-+]$//g' > {output}"
    #to use RealignerTargetCreator from GATK use the line below
    shell: "{JAVA} {JAVA_MEM} {JAVA_THREADS} -jar {GATK} -T RealignerTargetCreator --maxIntervalSize 10000 -R {input.ref} -I {input.bamFile} -o {output} > {log} 2>&1"


rule gatk_indelRealigner:
    input:
        ref = os.path.join(FASTA_DIR, os.path.basename(REFERENCE_FASTA)),
        ref_index = os.path.join(FASTA_DIR, ".".join([os.path.basename(REFERENCE_FASTA), 'fai'])),
        ref_dict = os.path.join(FASTA_DIR, ".".join([os.path.basename(REFERENCE_FASTA).replace(".fa", ""), 'dict'])),
        bamIndex = os.path.join(OUTPUT_DIR, "aln_split_by_indels", "{sample}.with_indels.bam.bai"),
        bamFile = os.path.join(OUTPUT_DIR, "aln_split_by_indels", "{sample}.with_indels.bam"),
        intervals = os.path.join(GATK_DIR, "{sample}.gatk_realigner.intervals")
    output:
        bam = os.path.join(GATK_DIR, "{sample}.indels.realigned.bam"),
        bai = os.path.join(GATK_DIR, "{sample}.indels.realigned.bai")
    log: os.path.join(LOG_DIR, "GATK", "gatk_realigner.{sample}.log")
    shell: "{JAVA} {JAVA_MEM} {JAVA_THREADS} -jar {GATK} -T IndelRealigner -maxReads 1000000 -R {input.ref} -I {input.bamFile} -targetIntervals {input.intervals} -o {output.bam} > {log} 2>&1"

# merge bam files that were split for realigning the indels
# merge the realigned bam file with the bam file containing no indels.
rule merge_bam:
    input:
        bam1 = os.path.join(GATK_DIR, "{sample}.indels.realigned.bam"),
        bam2 = os.path.join(OUTPUT_DIR, "aln_split_by_indels", "{sample}.without_indels.bam")
    output:
        bam=os.path.join(OUTPUT_DIR, "aln_merged", "{sample}.bam"),
        bai=os.path.join(OUTPUT_DIR, "aln_merged", "{sample}.bam.bai"),
    log: os.path.join(LOG_DIR, "SAMTOOLS", "merge.realigned.{sample}.log")
    shell:
        """
        samtools merge {output.bam} {input.bam1} {input.bam2}
        samtools index {output.bam}
        """

rule samtools_stats:
    input:
        bamfile = os.path.join(OUTPUT_DIR, "aln_merged", "{sample}.bam"),
        ref = os.path.join(FASTA_DIR, os.path.basename(REFERENCE_FASTA))
    output: os.path.join(OUTPUT_DIR, "SAMTOOLS", "{sample}.samtools.stats.txt")
    log: os.path.join(LOG_DIR, "SAMTOOLS", "samtools_stats.{sample}.log")
    shell: "samtools stats --reference {input.ref} {input.bamfile} > {output} 2> {log} 2>&1"

rule multiqc:
    input:
        fastqc = expand(os.path.join(FASTQC_DIR, "{sample}.fastqc.done"), sample = SAMPLES),
        trimmomatic = expand(os.path.join(TRIMMED_READS_DIR, "{sample}.fastq.gz"), sample = SAMPLES),
        samtools = expand(os.path.join(OUTPUT_DIR, "SAMTOOLS", "{sample}.samtools.stats.txt"), sample = SAMPLES)
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
        bamIndex = os.path.join(OUTPUT_DIR, "aln_merged",  "{sample}.bam.bai"),
        bamFile = os.path.join(OUTPUT_DIR,  "aln_merged", "{sample}.bam"),
        #bamIndex = os.path.join(MAPPED_READS_DIR,  "{sample}.bam.bai"),
        #bamFile = os.path.join(MAPPED_READS_DIR,  "{sample}.bam")
    output:
        # os.path.join(INDELS_DIR, "{sample}", "{sample}.indelScores.bigwig"),
        # os.path.join(INDELS_DIR, "{sample}", "{sample}.deletionScores.bigwig"),
        # os.path.join(INDELS_DIR, "{sample}", "{sample}.insertionScores.bigwig"),
        # os.path.join(INDELS_DIR, "{sample}", "{sample}.alnCoverage.bigwig"),
        # os.path.join(INDELS_DIR, "{sample}", "{sample}.deletions.bed"),
        # os.path.join(INDELS_DIR, "{sample}", "{sample}.insertions.bed"),
        # os.path.join(INDELS_DIR, "{sample}", "{sample}.indels.tsv"),
        # os.path.join(INDELS_DIR, "{sample}", "{sample}.reads_with_indels.tsv"),
        # os.path.join(INDELS_DIR, "{sample}", "{sample}.insertedSequences.tsv"),
        os.path.join(INDELS_DIR, "{sample}", "{sample}.sgRNA_efficiency.tsv")
    params:
        script = os.path.join(SRC_DIR, "src", "getIndelStats.R"),
        tech = lambda wildcards: lookup('sample_name', wildcards.sample, ['tech'])[0]
    log: os.path.join(LOG_DIR, "indel_stats", "getIndelStats.{sample}.log")
    shell: "{RSCRIPT} {params.script} {input.bamFile} {wildcards.sample} {INDELS_DIR} {CUT_SITES_FILE} {params.tech} > {log} 2>&1"

#prepare _site.yml and other Rmd files to be rendered into a html report (see renderSite rule)
rule generateSiteFiles:
    input: #TODO: print a file that flags all getIndelStats jobs are finished.
        expand(os.path.join(INDELS_DIR, "{sample}", "{sample}.sgRNA_efficiency.tsv"), sample = SAMPLES)
    output:
        os.path.join(REPORT_DIR, "_site.yml"),
        os.path.join(REPORT_DIR, "index.Rmd"),
        os.path.join(REPORT_DIR, "config.yml")
    params:
        report_scripts_dir = os.path.join(SRC_DIR, "src", "report_scripts"),
        script = os.path.join(SRC_DIR, "src", "generateSiteFiles.R")
    log: os.path.join(LOG_DIR, "generateSiteFiles.log")
    shell:
        "{RSCRIPT} {params.script} {params.report_scripts_dir} {SAMPLE_SHEET_FILE} {CUT_SITES_FILE} {COMPARISONS_FILE} {OUTPUT_DIR} {REPORT_DIR} {RSCRIPT} > {log} 2>&1"

rule renderSite:
    input:
        os.path.join(REPORT_DIR, "_site.yml"),
        os.path.join(REPORT_DIR, "index.Rmd"),
        os.path.join(REPORT_DIR, "config.yml")
    output:
        os.path.join(REPORT_DIR, "index.html")
    params:
        render_script = os.path.join(SRC_DIR, "src", "render_site.sh"),
        report_scripts_dir = os.path.join(SRC_DIR, "src", "report_scripts")
    log: os.path.join(LOG_DIR, "renderSite.log")
    shell:
        "{RSCRIPT} -e \"library(rmarkdown); rmarkdown::render_site(\'{REPORT_DIR}\')\" > {log} 2>&1"
#        "bash {params.render_script} {REPORT_DIR} {RSCRIPT} > {log} 2>&1"

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
