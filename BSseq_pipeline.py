#============================================================================================================
# SNAKEMAKE FILE WRITTEN BY THE AKALIN GROUP AT MDC, BERLIN, 2017
# Alexander Gosdschan, Katarzyna Wreczycka, Bren Osberg, Ricardo Wurmus
# To process bisulfite sequencing data from raw fastq files to performing integrated bioinformatics analysis.

# SUBMIT THIS JOB INTERACTIVELY WITH:
# import IPython;
# IPython.embed()
 
#============================================================================================================

#------ set config file, include function definitions, and set os:
import os
include   : "./rules/post_mapping.rules"
include   : "./scripts/func_defs.py"

#---------------------------     LIST THE OUTPUT DIRECTORIED AND SUBDIRECTORIED TO BE PRODUCED     ------------------------------

DIR_xmethed     = '07_xmethed/'
DIR_sorted      = '06_sorted/'
DIR_deduped     = '05_deduped/'
DIR_mapped      = '04_mapped/'
DIR_posttrim_QC = '03_posttrim_QC/'
DIR_trimmed     = '02_trimmed/'
DIR_rawqc       = '01_rawqc/'
DIR_annot       = 'annotation/'
DIR_diffmeth    = 'differential_methylation/'
DIR_final       = "final_Report/"


#---------------------------------     DEFINE PATHS AND FILE NAMES:  ----------------------------------

PATHIN          = "path_links/input/"       #--- location of the data files to be imported (script creates symbolic link) 
GENOMEPATH      = "path_links/refGenome/"   #--- where the reference genome being mapped to is stored
VERSION         = config["GENOME_VERSION"]  #--- version of the genome being mapped to.

bismark_cores=config["bismark_cores"]       #--- from config file. Gets passed to bismark multicore argument.  

#-------------------------------      DEFINE PROGRAMS TO BE EXECUTED: ---------------------------------

programs = config['programs']
FASTQC                         =  programs["FASTQC"]            #--- self-explanatory program names.
TRIMGALORE                     =  programs["TRIMGALORE"]
CUTADAPT                       =  programs["CUTADAPT"]
BISMARK_GENOME_PREPARATION     =  programs["BISMARK_GENOME_PREPARATION"]
BISMARK                        =  programs["BISMARK"]
BOWTIE2                        =  programs["BOWTIE2"]
DEDUPLICATE_BISMARK            =  programs["DEDUPLICATE_BISMARK"]
SAMTOOLS                       =  programs["SAMTOOLS"]


#---------------------------     LIST THE OUTPUT FILES TO BE PRODUCED     ------------------------------

# Below is a mapping of rule names to the expected output files they
# produce.  The desired output files are specified in
# "OUTPUT_FILES".  A different set of output files can be
# selected to run fewer rules.

all_output_files = {
    # These are expensive one-time rules to prepare the genome.
    'genome-prep-CT': GENOMEPATH+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
    'genome-prep-GA': GENOMEPATH+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",

    '01-raw-qc': [
        expand (list_files_rawQC(DIR_rawqc,
                                 config["SAMPLES"][sample]["files"],
                                 config["SAMPLES"][sample]["SampleID"]))
        for sample in config["SAMPLES"]],

    # This rule is always executed, as trimming is a prerequisite for
    # subsequent rules
    '02-trimgalore': [
        expand (list_files_TG(DIR_trimmed,
                              config["SAMPLES"][sample]["files"],
                              config["SAMPLES"][sample]["SampleID"]))
        for sample in config["SAMPLES"]],

    # fastQC output files are not needed downstream and need to be
    # called explicitly.
    '03-posttrim-qc': [
        expand (list_files_posttrim_QC(DIR_posttrim_QC,
                                       config["SAMPLES"][sample]["files"],
                                       config["SAMPLES"][sample]["SampleID"]))
        for sample in config["SAMPLES"]],

    '04-mapping': [
        expand (list_files_bismark(DIR_mapped,
                                   config["SAMPLES"][sample]["files"],
                                   config["SAMPLES"][sample]["SampleID"]))
        for sample in config["SAMPLES"]],

    '05-deduplication': [
        expand (list_files_dedupe(DIR_deduped,
                                  config["SAMPLES"][sample]["files"],
                                  config["SAMPLES"][sample]["SampleID"]))
        for sample in config["SAMPLES"]],

    '06-sorting': [
        expand (list_files_sortbam(DIR_sorted,
                                   config["SAMPLES"][sample]["files"],
                                   config["SAMPLES"][sample]["SampleID"]))
        for sample in config["SAMPLES"]],

     # TODO: had to add this part to call bam_methCall for diff meth rule
    'bam-processing': [
        expand (bam_processing(METHCALLDIR,
                               config["SAMPLES"][sample]["files"],
                               sample))
        for sample in config["SAMPLES"]],

    # differential methylation calling
    'diffmeth': [ DIR_diffmeth+"_".join(x)+".sorted_diffmeth.nb.html" for x in config["DIFF_METH"]],
		            
    # annotation diff meth cytosines
    'annotation': [ DIR_annot+"_".join(x)+".sorted_"+config["GENOME_VERSION"]+"_annotation.diff.meth.nb.html" for x in config["DIFF_METH"]],

    # final report
    # TODO: This needs to be editted once we determine what final reports we want to export!
    'final-report': [
        expand (Final(DIR_final,
                      config["SAMPLES"][sample]["files"],
                      VERSION,
                      config["SAMPLES"][sample]["SampleID"]))
        for sample in config["SAMPLES"]]
}

# Selected output files from the above set.
selected_rules = ['final-report']
OUTPUT_FILES = [all_output_files[rule] for rule in selected_rules]

#--- NICE gauges the computational burden, ranging from -19 to +19.
#--- The more "nice" you are, the more you allow other processes to jump ahead of you
#--- (like in traffic). Generally set to maximally nice=19 to avoid interference with other users.
def nice(cmd):
    return "nice -" + str(config["NICE"]) + " " + cmd

#--- In case you want to debug the code with interactive commands:
# import IPython;
# IPython.embed()
# print("Executing job to produce the following files: ")
# print("OUTPUT_FILES=")
# for x in OUTPUT_FILES: print( x)
#-------

# ==============================================================================================================
#
#                                         BEGIN RULES    
#
# rules are separated by "==" bars into pairs for paired-end and single-end (subdivided by smaller "--" dividers)
# ===============================================================================================================

rule all:
    input:
        OUTPUT_FILES


# ==========================================================================================
# sort the bam file:

rule sortbam_se:
    input:
        DIR_deduped+"{sample}_se_bt2.deduped.bam"
    output:
        DIR_sorted+"{sample}_se_bt2.deduped.sorted.bam"
    message: fmt("Sorting bam file {input}")
    shell:
        nice("{SAMTOOLS} sort {input} -o {output}")
#-----------------------
rule sortbam_pe:
    input:
        DIR_deduped+"{sample}_1_val_1_bt2.deduped.bam"
    output:
        DIR_sorted+"{sample}_1_val_1_bt2.deduped.sorted.bam"
    message: fmt("Sorting bam file {input}")
    shell:
        nice("{SAMTOOLS} sort {input} -o {output}")

# ==========================================================================================
# deduplicate the bam file:

rule deduplication_se:
    input:
        DIR_mapped+"{sample}_trimmed_bismark_bt2.bam"
    output:
        DIR_deduped+"{sample}_se_bt2.deduped.bam"
    params:
        bam="--bam ",
        sampath="--samtools_path "+SAMTOOLS
    log:
        DIR_deduped+"{sample}_deduplication.log"
    message: fmt("Deduplicating single-end aligned reads from {input}")
    shell:
        nice("{SAMTOOLS} rmdup {input}  {output} > {log} 2>&1 ")
#-----------------------
rule deduplication_pe:
    input:
        DIR_mapped+"{sample}_1_val_1_bismark_bt2_pe.bam"
    output:
        DIR_deduped+"{sample}_1_val_1_bt2.deduped.bam"
    log:
        DIR_deduped+"{sample}_deduplication.log"
    message: fmt("Deduplicating paired-end aligned reads from {input}")
    shell:
        nice("{SAMTOOLS} fixmate {input}  {output} > {log} 2>&1 ")

# ==========================================================================================
# align and map:
 
rule bismark_se:
    input:
        refconvert_CT = GENOMEPATH+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
	refconvert_GA = GENOMEPATH+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
        fqfile = DIR_trimmed+"{sample}_trimmed.fq.gz",
        qc     = DIR_posttrim_QC+"{sample}_trimmed_fastqc.html"
    output:
        DIR_mapped+"{sample}_trimmed_bismark_bt2.bam",
        DIR_mapped+"{sample}_trimmed_bismark_bt2_SE_report.txt"
    params:
        bismark_args = config.get("bismark_args",""),
        genomeFolder = "--genome_folder " + GENOMEPATH,
        outdir = "--output_dir  "+DIR_mapped,
        nucCov = "--nucleotide_coverage",
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(BOWTIE2) ,
        useBowtie2  = "--bowtie2 ",
        samtools    = "--samtools_path "+ os.path.dirname(SAMTOOLS),
        tempdir     = "--temp_dir "+DIR_mapped
    log:
        DIR_mapped+"{sample}_bismark_se_mapping.log"
    message: fmt("Mapping single-end reads to genome {VERSION}")
    shell:
        nice("{BISMARK} {params} --multicore "+bismark_cores+" {input.fqfile} > {log}  2>&1 ")

#-----------------------

rule bismark_pe:
    input:
        refconvert_CT = GENOMEPATH+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
	refconvert_GA = GENOMEPATH+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
        fin1 = DIR_trimmed+"{sample}_1_val_1.fq.gz",
        fin2 = DIR_trimmed+"{sample}_2_val_2.fq.gz",
        qc   = [ DIR_posttrim_QC+"{sample}_1_val_1_fastqc.html",
                 DIR_posttrim_QC+"{sample}_2_val_2_fastqc.html"]
    output:
        DIR_mapped+"{sample}_1_val_1_bismark_bt2_pe.bam",
        DIR_mapped+"{sample}_1_val_1_bismark_bt2_PE_report.txt"
    params:
        bismark_args = config.get("bismark_args",""),
        genomeFolder = "--genome_folder " + GENOMEPATH,
        outdir = "--output_dir  "+DIR_mapped,
        nucCov = "--nucleotide_coverage",
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(BOWTIE2) ,
        useBowtie2  = "--bowtie2 ",
        samtools    = "--samtools_path "+ os.path.dirname(SAMTOOLS),
        tempdir     = "--temp_dir "+DIR_mapped
    log:
        DIR_mapped+"{sample}_bismark_pe_mapping.log"
    message: fmt("Mapping paired-end reads to genome {VERSION}.")
    shell:
        nice("{BISMARK} {params} --multicore "+bismark_cores+" -1 {input.fin1} -2 {input.fin2} > {log} 2>&1 ")


# ==========================================================================================
# generate reference genome:

rule bismark_genome_preparation:
    input:
        GENOMEPATH
    output:
        GENOMEPATH+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
        GENOMEPATH+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
    params:
        bismark_genome_preparation_args = config.get("bismark_genome_preparation",""),
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(BOWTIE2) ,
        useBowtie2 = "--bowtie2 ",
        verbose = "--verbose "
    log:
        'bismark_genome_preparation_'+VERSION+'.log'
    message: fmt("Converting {VERSION} Genome into Bisulfite analogue")
    shell:
        nice("{BISMARK_GENOME_PREPARATION} {params} {input} > {log} 2>&1 ")

# ==========================================================================================
# post-trimming quality control

rule fastqc_after_trimming_se:
    input:
        DIR_trimmed+"{sample}_trimmed.fq.gz",
    output:
    	DIR_posttrim_QC+"{sample}_trimmed_fastqc.html",
    	DIR_posttrim_QC+"{sample}_trimmed_fastqc.zip"
    params:
        fastqc_args = config.get("fastqc_args", ""),
        outdir = "--outdir "+DIR_posttrim_QC
    log:
   	    DIR_posttrim_QC+"{sample}_trimmed_fastqc.log"
    message: fmt("Quality checking trimmmed single-end data from {input}")
    shell:
        nice("{FASTQC} {params.outdir} {input} > {log} 2>&1 ")
#--------
rule fastqc_after_trimming_pe:
    input:
        DIR_trimmed+"{sample}_1_val_1.fq.gz",
        DIR_trimmed+"{sample}_2_val_2.fq.gz"
    output:
    	DIR_posttrim_QC+"{sample}_1_val_1_fastqc.html",
    	DIR_posttrim_QC+"{sample}_1_val_1_fastqc.zip",
    	DIR_posttrim_QC+"{sample}_2_val_2_fastqc.zip",
        DIR_posttrim_QC+"{sample}_2_val_2_fastqc.html"
    params:
        fastqc_args = config.get("fastqc_args", ""),
        outdir = "--outdir "+DIR_posttrim_QC
    log:
   	    DIR_posttrim_QC+"{sample}_trimmed_fastqc.log"
    message: fmt("Quality checking trimmmed paired-end data from {input}")
    shell:
        nice("{FASTQC} {params.outdir} {input} > {log} 2>&1 ")

# ==========================================================================================
# trim the reads

rule trimgalore_se:
    input:
       qc   = DIR_rawqc+"{sample}_fastqc.html",
       file = PATHIN+"{sample}.fq.gz"
    output:
       DIR_trimmed+"{sample}_trimmed.fq.gz" #---- this ALWAYS outputs .fq.qz format.
    params:
       extra          = config.get("trim_galore_args", ""),
       outdir = "--output_dir "+DIR_trimmed,
       phred = "--phred33",
       gz = "--gzip",
       cutadapt = "--path_to_cutadapt "+CUTADAPT,
    log:
       DIR_trimmed+"{sample}.trimgalore.log"
    message: fmt("Trimming raw single-end read data from {input}")
    shell:
       nice("{TRIMGALORE} {params} {input.file} > {log} 2>&1 ")

#-----------------------
rule trimgalore_pe:
    input:
        qc    = [ DIR_rawqc+"{sample}_1_fastqc.html",
                  DIR_rawqc+"{sample}_2_fastqc.html"],
        files = [ PATHIN+"{sample}_1.fq.gz",
                  PATHIN+"{sample}_2.fq.gz"]
    output:
        DIR_trimmed+"{sample}_1_val_1.fq.gz", #---- this ALWAYS outputs .fq.qz format.
        DIR_trimmed+"{sample}_2_val_2.fq.gz",
    params:
        extra          = config.get("trim_galore_args", ""),
        outdir         = "--output_dir "+DIR_trimmed,
        phred          = "--phred33",
        gz             = "--gzip",
        cutadapt       = "--path_to_cutadapt "+CUTADAPT,
        paired         = "--paired"
    log:
        DIR_trimmed+"{sample}.trimgalore.log"
    message:
        fmt("Trimming raw paired-end read data from {input}")
    shell:
        nice("{TRIMGALORE} {params} {input.files} > {log} 2>&1 ")

# ==========================================================================================
# raw quality control

rule fastqc_raw: #----only need one: covers BOTH pe and se cases.
    input:
        PATHIN+"{sample}.fq.gz"
    output:
        DIR_rawqc+"{sample}_fastqc.html",
        DIR_rawqc+"{sample}_fastqc.zip"
    params:
        fastqc_args = config.get("fastqc_args", ""),
        outdir = "--outdir "+ DIR_rawqc     # usually pass params as strings instead of wildcards.
    log:
        DIR_rawqc+"{sample}_fastqc.log"
    message: fmt("Quality checking raw read data from {input}")
    shell:
        nice("{FASTQC} {params}  {input} > {log} 2>&1 ")
