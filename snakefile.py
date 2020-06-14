"""
Snakefile for crispr-DART pipeline
"""
import sys
#print(sys.version)
import os
import yaml
import pandas as pd
from itertools import chain


# tools
RSCRIPT = config['tools']['Rscript']
JAVA = config['tools']['java']['path']
JAVA_MEM = config['tools']['java']['mem']
JAVA_THREADS = config['tools']['java']['threads']
GATK = config['tools']['gatk']

# input locations
SRC_DIR = os.path.abspath(config['source-dir'])
READS_DIR = os.path.abspath(config['reads-dir'])
SAMPLE_SHEET_FILE = os.path.abspath(config['sample_sheet'])
CUT_SITES_FILE = os.path.abspath(config['cutsites'])
COMPARISONS_FILE = os.path.abspath(config.get('comparisonsFile', {}))
REFERENCE_FASTA = os.path.abspath(config['reference_fasta'])

#output locations
OUTPUT_DIR = os.path.abspath(config['output-dir'])
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
SAMPLE_SHEET = pd.read_csv(SAMPLE_SHEET_FILE)
TARGET_NAMES = list(set(SAMPLE_SHEET['target_name'].tolist()))
#print(TARGET_NAMES)
# get unique rows for only considering sample id/name/read files
# (one sample may contain multiple target region/name fields, but we
# don't need to process the same read files for each target region)
SAMPLE_SHEET_tmp = SAMPLE_SHEET[['sample_name', 'reads', 'reads2', 'tech']]
SAMPLE_SHEET = SAMPLE_SHEET_tmp.drop_duplicates()

SAMPLES = SAMPLE_SHEET['sample_name'].tolist()
#print(SAMPLES)

# look up values from sample sheet (assuming it is a data frame) for multiple fields
def lookup(column, predicate, fields=[]):
    values = [SAMPLE_SHEET[SAMPLE_SHEET[column] == predicate][f].tolist() for f in fields]
    values = list(chain.from_iterable(values))
    #remove nan values
    return([f for f in values if str(f) != 'nan'])

def reads_input(wc):
  sample = wc.sample
  files = [os.path.join(READS_DIR, f) for f in lookup('sample_name', sample, ['reads', 'reads2']) if f]
  #print(files)
  return files

def get_bbmap_command(wc):
    sample = wc.sample
    tech = lookup('sample_name', sample, ['tech'])[0]
    if tech == 'illumina':
        return('bbmap.sh')
    elif tech == 'pacbio':
        return('mapPacBio.sh maxlen=6000')

# determine if the sample library is single end or paired end
def libType(wc):
  sample = wc.sample
  files = lookup('sample_name', sample, ['reads', 'reads2'])
  #print(files)
  count = sum(1 for f in files if f)
  if count == 2:
      return 'paired'
  elif count == 1:
      return 'single'

def map_input(wc):
  sample = wc.sample
  if libType(wc) == 'paired':
    return [os.path.join(TRIMMED_READS_DIR, "{sample}_val_1.fq.gz".format(sample=sample)),
            os.path.join(TRIMMED_READS_DIR, "{sample}_val_2.fq.gz".format(sample=sample))]
  elif libType(wc) == 'single':
    return [os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed.fq.gz".format(sample=sample))]


rule all:
    input:
        #expand(os.path.join(FASTQC_DIR, "{sample}.fastqc.done"), sample = SAMPLES),
        #expand(os.path.join(MAPPED_READS_DIR, "{sample}.bam.bai"), sample = SAMPLES),
        #expand(os.path.join(GATK_DIR, "{sample}.indels.realigned.bam"), sample = SAMPLES),
        expand(os.path.join(OUTPUT_DIR,  "aln_merged", "{sample}.bam"), sample = SAMPLES),
        #expand(os.path.join(OUTPUT_DIR, "SAMTOOLS", "{sample}.samtools.stats.txt"), sample = SAMPLES),
        #os.path.join(OUTPUT_DIR, "multiqc", "multiqc_report.html"),
        #expand(os.path.join(INDELS_DIR, "{sample}", "{sample}.sgRNA_efficiency.tsv"), sample = SAMPLES),
        os.path.join(REPORT_DIR, "index.html"),
        #expand(os.path.join(REPORT_DIR, "{target}.CoverageProfiles.html"), target = TARGET_NAMES),
        #expand(os.path.join(REPORT_DIR, "{target}.SampleComparisons.html"), target = TARGET_NAMES),
        #expand(os.path.join(REPORT_DIR, "{target}.sgRNA_efficiency_stats.html"), target = TARGET_NAMES),
        #expand(os.path.join(REPORT_DIR, "{target}.Indel_Diversity.html"), target = TARGET_NAMES)
        #get_output_file_list(INDELS_DIR, "freeBayes_variants.vcf"),
        #get_output_file_list(INDELS_DIR, "freeBayes_deletions.bed"),
        #get_output_file_list(INDELS_DIR, "freeBayes_insertions.bed"),

#notice that get_amplicon_file function for 'fasta' should only be used once.
#Other rules that need the amplicon fasta sequence as input should use :
#lambda wildcards: os.path.join(FASTA_DIR, ''.join([wildcards.amplicon, ".fasta"]))
rule reformatFasta:
    input: REFERENCE_FASTA
    output: os.path.join(FASTA_DIR, os.path.basename(REFERENCE_FASTA))
    log: os.path.join(LOG_DIR, ".".join(["reformatFasta", os.path.basename(REFERENCE_FASTA), "log"]))
    params:
        memory = config['tools']['reformat']['memory']
    shell:
        "reformat.sh {params.memory} in={input} out={output} tuc > {log} 2>&1"

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

rule trim_galore_pe:
  input: reads_input
  output:
    r1=os.path.join(TRIMMED_READS_DIR, "{sample}_val_1.fq.gz"),
    r2=os.path.join(TRIMMED_READS_DIR, "{sample}_val_2.fq.gz")
  log: os.path.join(LOG_DIR, "TRIM", "trimgalore.{sample}.log")
  shell: "trim_galore -o {TRIMMED_READS_DIR} --cores 2 --basename {wildcards.sample} --paired {input[0]} {input[1]} >> {log} 2>&1"

rule trim_galore_se:
  input: reads_input
  output: os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed.fq.gz"),
  log: os.path.join(LOG_DIR, "TRIM", "trimgalore.{sample}.log")
  params:
    tech = lambda wildcards: lookup('sample_name', wildcards[0], ['tech'])[0]
  run:
    if params.tech == 'illumina':
        shell("trim_galore -o {TRIMMED_READS_DIR} --cores 2 --basename {wildcards.sample} {input[0]} >> {log} 2>&1")
    elif params.tech == 'pacbio':
        shell("cp {input} {output}")

rule bbmap_indexgenome:
    input: os.path.join(FASTA_DIR, os.path.basename(REFERENCE_FASTA))
    output: os.path.join(BBMAP_INDEX_DIR, re.sub("\.fa(sta)?$", "", os.path.basename(REFERENCE_FASTA)))
    log: os.path.join(LOG_DIR, ".".join(["bbmap_index", os.path.basename(REFERENCE_FASTA), "log"]))
    shell: "bbmap.sh t=10 ref={input} path={output} > {log} 2>&1"

rule bbmap_map:
    input:
        reads = map_input,
        ref = os.path.join(BBMAP_INDEX_DIR, re.sub("\.fa(sta)?$", "", os.path.basename(REFERENCE_FASTA)))
    output:
        os.path.join(MAPPED_READS_DIR, "{sample}", "{sample}.sam")
    params:
        aligner = lambda wildcards: get_bbmap_command(wildcards),
        memory = config['tools']['bbmap']['memory'],
        options = config['tools']['bbmap']['options'],
        libtype = lambda wildcards: libType(wildcards)
    log: os.path.join(LOG_DIR, "BBMAP", "bbmap_align.{sample}.log")
    run:
        if params.libtype == 'single':
            shell("{params.aligner} {params.memory} {params.options} path={input.ref} in={input.reads} outm={output} > {log} 2>&1")
        elif params.libtype == 'paired':
            shell("{params.aligner} {params.memory} {params.options} keepnames=t path={input.ref} in1={input.reads[0]} in2={input.reads[1]} outm={output}> {log} 2>&1")


# GATK requires read groups, so here we add some dummy read group information to the bam files
rule samtools_addReadGroups:
    input: os.path.join(MAPPED_READS_DIR,  "{sample}", "{sample}.sam")
    output: os.path.join(MAPPED_READS_DIR, "{sample}", "{sample}.sam_withreadgroups")
    log: os.path.join(LOG_DIR, "SAMTOOLS", "samtools_addReadGroups.{sample}.log")
    shell:
        """
        samtools addreplacerg -m overwrite_all -o {output} -r ID:crips_dart -r PL:illumina -r SM:{wildcards.sample} {input} > {log} 2>&1
        rm {input}
        """

rule samtools_sam2bam:
    input: os.path.join(MAPPED_READS_DIR,  "{sample}", "{sample}.sam_withreadgroups")
    output: os.path.join(MAPPED_READS_DIR,  "{sample}", "{sample}.bam")
    log: os.path.join(LOG_DIR, "SAMTOOLS", "samtools_sam2bam.{sample}.log")
    shell:
        """
        samtools view -bh {input} | samtools sort -o {output} > {log} 2>&1
        rm {input}
        """

rule samtools_indexbam:
    input: os.path.join(MAPPED_READS_DIR,  "{sample}", "{sample}.bam")
    output: os.path.join(MAPPED_READS_DIR,  "{sample}", "{sample}.bam.bai")
    log: os.path.join(LOG_DIR, "samtools_indexbam.{sample}.log")
    shell: "samtools index {input} > {log} 2>&1"


# split the bam file into two: First contains only reads alignments with indels
# second contains alignments without indels
rule split_bam_by_indels:
    input:
        bamFile = os.path.join(MAPPED_READS_DIR,  "{sample}", "{sample}.bam"),
        bamIndex = os.path.join(MAPPED_READS_DIR,  "{sample}", "{sample}.bam.bai")
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
    log: os.path.join(LOG_DIR, "GATK", "gatk_realigner_target_creator.{sample}.log")
    shell: "{JAVA} {JAVA_MEM} {JAVA_THREADS} -jar {GATK} -T RealignerTargetCreator --maxIntervalSize 10000 -R {input.ref} -I {input.bamFile} -o {output} > {log} 2>&1"

# only consider realigner intervals that overlap target region
rule subset_gatk_realigner_intervals:
    input:
        os.path.join(GATK_DIR, "{sample}.gatk_realigner.intervals")
    output:
        os.path.join(GATK_DIR, "{sample}.gatk_realigner.target.intervals")
    params:
        script = os.path.join(SRC_DIR, "src", "subset_gatk_realigner_intervals.R"),
    log: os.path.join(LOG_DIR, "GATK", "subset_gatk_realigner_intervals.{sample}.log")
    shell: "{RSCRIPT} {params.script} {SAMPLE_SHEET_FILE} {wildcards.sample} {input} > {log} 2>&1"



rule gatk_indelRealigner:
    input:
        ref = os.path.join(FASTA_DIR, os.path.basename(REFERENCE_FASTA)),
        ref_index = os.path.join(FASTA_DIR, ".".join([os.path.basename(REFERENCE_FASTA), 'fai'])),
        ref_dict = os.path.join(FASTA_DIR, ".".join([os.path.basename(REFERENCE_FASTA).replace(".fa", ""), 'dict'])),
        bamIndex = os.path.join(OUTPUT_DIR, "aln_split_by_indels", "{sample}.with_indels.bam.bai"),
        bamFile = os.path.join(OUTPUT_DIR, "aln_split_by_indels", "{sample}.with_indels.bam"),
        intervals = os.path.join(GATK_DIR, "{sample}.gatk_realigner.target.intervals")
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
        samtools merge {output.bam} {input.bam1} {input.bam2} 2> {log} 2>&1
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
        trim = expand(os.path.join(TRIMMED_READS_DIR, "{sample}.fastq.gz"), sample = SAMPLES),
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
        os.path.join(REPORT_DIR, "config.yml"),
        #expand(os.path.join(REPORT_DIR, "{target}.CoverageProfiles.Rmd"), target = TARGET_NAMES),
        #expand(os.path.join(REPORT_DIR, "{target}.SampleComparisons.Rmd"), target = TARGET_NAMES),
        #expand(os.path.join(REPORT_DIR, "{target}.sgRNA_efficiency_stats.Rmd"), target = TARGET_NAMES),
        #expand(os.path.join(REPORT_DIR, "{target}.Indel_Diversity.Rmd"), target = TARGET_NAMES)
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
        os.path.join(REPORT_DIR, "config.yml"),
        #expand(os.path.join(REPORT_DIR, "{target}.CoverageProfiles.Rmd"), target = TARGET_NAMES),
        #expand(os.path.join(REPORT_DIR, "{target}.SampleComparisons.Rmd"), target = TARGET_NAMES),
        #expand(os.path.join(REPORT_DIR, "{target}.sgRNA_efficiency_stats.Rmd"), target = TARGET_NAMES),
        #expand(os.path.join(REPORT_DIR, "{target}.Indel_Diversity.Rmd"), target = TARGET_NAMES)
    output:
        os.path.join(REPORT_DIR, "index.html"),
        #expand(os.path.join(REPORT_DIR, "{target}.CoverageProfiles.html"), target = TARGET_NAMES),
        #expand(os.path.join(REPORT_DIR, "{target}.SampleComparisons.html"), target = TARGET_NAMES),
        #expand(os.path.join(REPORT_DIR, "{target}.sgRNA_efficiency_stats.html"), target = TARGET_NAMES),
        #expand(os.path.join(REPORT_DIR, "{target}.Indel_Diversity.html"), target = TARGET_NAMES)
    params:
        render_script = os.path.join(SRC_DIR, "src", "render_site.sh"),
        report_scripts_dir = os.path.join(SRC_DIR, "src", "report_scripts")
    log: os.path.join(LOG_DIR, "renderSite.log")
    shell:
        "{RSCRIPT} -e \"library(rmarkdown); rmarkdown::render_site(\'{REPORT_DIR}\')\" > {log} 2>&1"
#        "bash {params.render_script} {REPORT_DIR} {RSCRIPT} > {log} 2>&1"
