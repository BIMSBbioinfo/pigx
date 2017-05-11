#!/usr/bin/env python3.5
# ---last updated on  Thu May 11 15:21:47 CEST 2017  by  blosberg  at location  , BrensMB.local

#  changes from  Thu May 11 15:21:47 CEST 2017 : removed some rules that were specific to deconvolution (this branch is now the general one for everyone's use --ends with sorted .bam); also cleaned up some redundant commentary

#============================================================================================================
# SNAKEMAKE FILE WRITTEN BY THE AKALIN GROUP AT MDC, BERLIN, 2017
# Alexander Gosdschan, Katarzyna Wreczycka, Bren Osberg
# To process bisulfite sequencing data from raw fastq files to performing integrated bioinformatics analysis.

# SUBMIT THIS JOB INTERACTIVELY WITH:
# $nohup  snakemake -s [this filename] --jobs [# of jobs to submit] > [logfilename] &
# You can also add the following options for cluster submission: --cluster "qsub -V -l h_vmem={params.mem} -pe smp {params.threads} -l h_rt=36:00:00" &

# import IPython;
# IPython.embed()
 
#============================================================================================================

#------ set config file, include function definitions, and set os:
import os
configfile: "./config.json"
include   : "./func_defs.py"

NICE=config["NICE"]
#--- NICE is an option to gauge the burden on computational resources, ranges from -19 to +19. 
#--- The more "nice" you are, the more you allow other processes to jump ahead of you 
#--- (like in traffic). Generally set to maximally nice=19 to avoid interference with system processes.

#---------------------------------     DEFINE PATHS AND FILE NAMES:  ----------------------------------

RCODE           = config["RCODE"]                   #--- the string that denotes which read the file corresponds to (only relevent for paired-end experiments.)
PATHIN          = config["paths"]["PATHIN"]         #--- location of the data files to be imported
PATHOUT         = config["paths"]["PATHOUT"]        #--- where to send the output
GTOOLBOX        = config["paths"]["GTOOLBOX"]       #--- where the programs are stored to carry out the necessary operations
GENOMEPATH      = config["paths"]["GENOMEPATH"]     #--- where the reference genome being mapped to is stored
LOGS            = config["paths"]["LOGS"]           #--- subfolder name for the logs of some programs

INEXT           = config["INEXT"]                   #--- input file extension; usually .fq.gz, but can also be .bz2 among other possibilities.
VERSION         = config["genomedat"]["VERSION"]        #--- version of the genome being mapped to.

CHROM_INFO      = config["genomedat"]["CHROM_INFO"]     #--- details of the reference genome (length, etc.) haploid chroms have been removed.
NUMTHREADS      = config["NUMTHREADS"]


# --- Directional adapters/reads are assumed by default:
if ( config["directional"] ):
    NON_DIR_FLAG=""
else:
    NON_DIR_FLAG=" --non_directional "

#-------------------------------      DEFINE PROGRAMS TO BE EXECUTED: ---------------------------------

FASTQC                         =  GTOOLBOX+config["progs"]["FASTQC"]            #--- self-explanatory program names.
TRIMGALORE                     =  GTOOLBOX+config["progs"]["TRIMGALORE"]
CUTADAPT                       =  GTOOLBOX+config["progs"]["CUTADAPT"]
BISMARK_GENOME_PREPARATION     =  GTOOLBOX+config["progs"]["BISMARK_GENOME_PREPARATION"]
BISMARK                        =  GTOOLBOX+config["progs"]["BISMARK"]
BOWTIE2                        =  GTOOLBOX+config["progs"]["BOWTIE2"]
DEDUPLICATE_BISMARK            =  GTOOLBOX+config["progs"]["DEDUPLICATE"]
BISMARK_METHYLATION_EXTRACTOR  =  GTOOLBOX+config["progs"]["BISMARK_METHYLATION_EXTRACTOR"]
BISMARK2REPORT                 =  GTOOLBOX+config["progs"]["BISMARK2REPORT"]

SAMTOOLS                       =  GTOOLBOX+config["progs"]["SAMTOOLS"] 

#---------------------------     LIST THE OUTPUT FILES TO BE PRODUCED     ------------------------------

# --- Below is the list of expected output files. They are enumerated by their sequence in the rules of the processing pipeline
# --- the process can be terminated earlier by expressing (i.e. uncommenting) only the [expand] commands corresponding to the 
# --- last rule that you wish to have executed.

OUTPUT_FILES = [
                #               ======  rule 01 raw QC    =========
                [ expand (list_files(PATHOUT+"01_rawqc/", config["SAMPLES"][sampleID]["files"], "_fastqc.html")  ) for sampleID in config["SAMPLES"]  ],
                # [ expand (list_files(PATHOUT+"01_rawqc/", config["SAMPLES"][sampleID]["files"], "_fastqc.zip" )  ) for sampleID in config["SAMPLES"]  ],

                #----RULE 2 IS ALWAYS EXECUTED, TRIMMING IS A PREREQUISITE FOR SUBSEQUENT RULES ----
                
                #               ======  rule 03 posttrim_QC_ ======
                [ expand ( list_files_posttrim_QC(PATHOUT+"03_posttrim_QC/", config["SAMPLES"][sampleID]["files"],".html")  ) for sampleID in config["SAMPLES"]  ],
                #--- fastQC output files are not needed downstream and need to be called explicitly.

                
                #-----  here is a list of intermediary files to be uncommented back into execution if you want to stop part-way along the process -----------------
                #               ====rule 02 trimgalore ======
                #               [ expand ( list_files_TG( PATHOUT+"02_trimmed/", config["SAMPLES"][sampleID]["files"] )  ) for sampleID in config["SAMPLES"]  ],
                #               [ expand ( list_files_TG( PATHOUT+"02_trimmed/", config["SAMPLES"][sampleID]["files"] )  ) for sampleID in config["SAMPLES"]  ],                
                #               ====rule 04 Mapping ======
                [ expand ( list_files_bismark(PATHOUT+"04_mapped/", config["SAMPLES"][sampleID]["files"] )  ) for sampleID in config["SAMPLES"]  ],
              
                #               ====formerly rule 05 Deduplication ======
                [ expand ( list_files_dedupe(PATHOUT+"05_deduped/", config["SAMPLES"][sampleID]["files"] )  ) for sampleID in config["SAMPLES"]  ],                                

                #               ====rule 06 sorting ======
                [ expand ( list_files_sortbam(PATHOUT+"06_sorted/", config["SAMPLES"][sampleID]["files"] )  ) for sampleID in config["SAMPLES"]  ],
                
                # ==================  FINAL REPORT =========================
                # @@@! This needs to be editted once we determine what final reports we want to export!
		# [ expand (PATHOUT+config["SAMPLES"][sampleID]["files"][0]+SEPEstr(config["SAMPLES"][sampleID]["files"] )+"_report.html"  ) for sampleID in config["SAMPLES"]  ],
                
]

#--- In case you want to debug the code with interactive commands:
# import IPython;
# IPython.embed()
# print("Executing job to produce the following files: ")
# for x in OUTPUT_FILES: print( x)
#--- 

# =======================================================================================================
#
#                                         BEGIN RULES    
#
# =======================================================================================================

rule all:
    input:
        OUTPUT_FILES

# --------------------------------------------------------------------------------
# rule clean:
#    shell: "if [ -d {PATHOUT} ]; then rm -r {PATHOUT}; fi"
# ==========================================================================================
rule sortbam_se:
    input:
        PATHOUT+"05_deduped/{sample}_se.deduplicated.bam"
    output:
        PATHOUT+"06_sorted/{sample}_se.deduplicated.sorted.bam"
    shell:
        "nice -"+str(NICE)+" samtools sort {input} -o {output}"

#-------- TWO SORTING STEPS ARE REQUIRED IN PAIRED END (see bismark deduplication doc) -----
rule sortbam_pe:
    input:
        PATHOUT+"05_deduped/{sample}"+RCODE+"1_val_1.deduplicated.bam"
    output:
        PATHOUT+"06_sorted/{sample}"+RCODE+"1_val_1.deduplicated.sorted.bam"
    shell:
        "nice -"+str(NICE)+" samtools sort {input} -o {output}"

# ==========================================================================================
rule bismark_se_deduplication:
    input:
        PATHOUT+"04_mapped/{sample}_trimmed_bismark_bt2.bam"
    output:
        PATHOUT+"05_deduped/{sample}_se.deduplicated.bam"
    params:
        bam="--bam ",
        sampath="--samtools_path "+SAMTOOLS
    log:
        PATHOUT+"05_deduped/{sample}_deduplication.log"
    message: """-----------   Deduplicating single-end read alignments ---------------------- """
    shell:
        "nice -"+str(NICE)+" samtools rmdup {input}  {output} 2> {log}"
#--------
rule bismark_pe_deduplication:
    input:
        PATHOUT+"04_mapped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.bam"
    output:
        PATHOUT+"05_deduped/{sample}"+RCODE+"1_val_1.deduplicated.bam"
    log:
        PATHOUT+"05_deduped/{sample}_deduplication.log"
    message: """-----------   Deduplicating paired-end read alignments ---------------------- """
    shell:
        "nice -"+str(NICE)+" samtools rmdup {input}  {output} 2> {log}"
# ==========================================================================================
# Align and map:

rule bismark_se:
    input:
       PATHOUT+"02_trimmed/{sample}_trimmed.fq.gz"
    output:
        PATHOUT+"04_mapped/{sample}_trimmed_bismark_bt2.bam",
        PATHOUT+"04_mapped/{sample}_trimmed_bismark_bt2.nucleotide_stats.txt",
        PATHOUT+"04_mapped/{sample}_trimmed_bismark_bt2_SE_report.txt"
    threads: 2
    params:
        N = "-N 1",
        L = "-L 20",
        genomeFolder = "--genome_folder " + GENOMEPATH,
        outdir = "--output_dir  "+PATHOUT+"04_mapped/",
        nucCov = "--nucleotide_coverage",
	    nonDir = NON_DIR_FLAG,             #--- THIS IS EMPTY IF DATA IS DIRECTIONAL (DEFAULT) OTHERWISE "non_directional". 
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(BOWTIE2) ,
        useBowtie2  = "--bowtie2 ",
        samtools    = "--samtools_path "+ os.path.dirname(SAMTOOLS),
        tempdir     = "--temp_dir "+PATHOUT
    log:
        PATHOUT+"04_mapped/{sample}_bismark_se_mapping.log"
    message: """-------------   Mapping single-end reads to genome {VERSION}. ------------- """
    shell:
        "nice -"+str(NICE)+" {BISMARK} {params} --multicore {threads} {input} 2> {log}"

#--------
rule bismark_pe:
    input:
        fin1 = PATHOUT+"02_trimmed/{sample}"+RCODE+"1_val_1.fq.gz",
        fin2 = PATHOUT+"02_trimmed/{sample}"+RCODE+"2_val_2.fq.gz"
    output:
        PATHOUT+"04_mapped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.bam",
        PATHOUT+"04_mapped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.nucleotide_stats.txt",
        PATHOUT+"04_mapped/{sample}"+RCODE+"1_val_1_bismark_bt2_PE_report.txt"
    threads: 2
    params:
        N = "-N 1",
        L = "-L 20",
        genomeFolder = "--genome_folder " + GENOMEPATH,
        outdir = "--output_dir  "+PATHOUT+"04_mapped/",
        nucCov = "--nucleotide_coverage",
        nonDir = NON_DIR_FLAG,             #--- THIS IS EMPTY IF DATA IS DIRECTIONAL (DEFAULT) OTHERWISE "non_directional".
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(BOWTIE2) ,
        useBowtie2  = "--bowtie2 ",
        samtools    = "--samtools_path "+ os.path.dirname(SAMTOOLS),
        tempdir     = "--temp_dir "+PATHOUT
    log:
        PATHOUT+"04_mapped/{sample}_bismark_pe_mapping.log"
    message: """-------------   Mapping paired-end reads to genome {VERSION}. ------------- """
    shell:
        "nice -"+str(NICE)+" {BISMARK} {params} --multicore {threads} -1 {input.fin1} -2 {input.fin2} 2> {log}"

# ==========================================================================================
# generate reference genome: ----  THIS ONLY GETS INVOKED WHEN MANUALLY CALLED SPECIFICIALLY ------

rule bismark_genome_preparation:
    input:
        GENOMEPATH
    output:
        GENOMEPATH+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
        GENOMEPATH+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
    params:
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(BOWTIE2) ,
        useBowtie2 = "--bowtie2 ",
        verbose = "--verbose "
    log:
        PATHOUT+'bismark_genome_preparation_'+VERSION+'.log'
    message: """ --------  converting {VERSION} Genome into Bisulfite analogue ------- """
    shell:
        "nice -"+str(NICE)+" {BISMARK_GENOME_PREPARATION} {params} {input} 2> {log}"

# ==========================================================================================
# post-trimming quality control

rule fastqc_after_trimming_se:
    input:
        PATHOUT+"02_trimmed/{sample}_trimmed.fq.gz",
    output:
    	PATHOUT+"03_posttrim_QC/{sample}_trimmed_fastqc.html",
    	PATHOUT+"03_posttrim_QC/{sample}_trimmed_fastqc.zip"
    params:
        outdir = "--outdir "+PATHOUT+"03_posttrim_QC/"
    log:
   	    PATHOUT+"03_posttrim_QC/{sample}_trimmed_fastqc.log"
    message: """ ------------  Quality checking trimmmed single-end data with Fastqc ------------- """
    shell:
        "nice -"+str(NICE)+" {FASTQC} {params.outdir} {input} 2> {log}"
#--------
rule fastqc_after_trimming_pe:
    input:
        PATHOUT+"02_trimmed/{sample}"+RCODE+"1_val_1.fq.gz",
        PATHOUT+"02_trimmed/{sample}"+RCODE+"2_val_2.fq.gz"
    output:
    	PATHOUT+"03_posttrim_QC/{sample}"+RCODE+"1_val_1_fastqc.html",
    	PATHOUT+"03_posttrim_QC/{sample}"+RCODE+"1_val_1_fastqc.zip",
    	PATHOUT+"03_posttrim_QC/{sample}"+RCODE+"2_val_2_fastqc.zip",
        PATHOUT+"03_posttrim_QC/{sample}"+RCODE+"2_val_2_fastqc.html"
    params:
        outdir = "--outdir "+PATHOUT+"03_posttrim_QC/"
    log:
   	    PATHOUT+"03_posttrim_QC/{sample}_trimmed_fastqc.log"
    message: """ ------------  Quality checking trimmmed paired-end data with Fastqc ------------- """
    shell:
        "nice -"+str(NICE)+" {FASTQC} {params.outdir} {input} 2> {log}"
#
# ==========================================================================================
# trim the reads

rule trimgalore_se:
   input:
       PATHIN+"{sample}"+INEXT
   output:
       PATHOUT+"02_trimmed/{sample}_trimmed.fq.gz" #---- this ALWAYS outputs .fq.qz format.
#       PATHOUT+"02_trimmed/{sample}"+INEXT+"_trimming_report.txt" #---@ commented out to avoid ambiguity: should find a way to put this back in while specififying it applies to any trimming report that does *NOT* contain "RCODE"
   params:
       outdir = "--output_dir "+PATHOUT+"02_trimmed/",
       phred = "--phred33",
       gz = "--gzip",
       cutadapt = "--path_to_cutadapt "+CUTADAPT,
       FivePrimeClip = "--clip_R1 10",
       ThreePrimeClip = "--three_prime_clip_R1 10"
   log:
       PATHOUT+"02_trimmed/{sample}.trimgalore.log"
   message:
       " ---------  Trimming raw single-end read data using {TRIMGALORE} -------  "
   shell:
       "nice -"+str(NICE)+" {TRIMGALORE} {params} {input} 2> {log}"

#-----------------------
rule trimgalore_pe:
    input:
        PATHIN+"{sample}"+RCODE+"1"+INEXT,
        PATHIN+"{sample}"+RCODE+"2"+INEXT
    output:
        PATHOUT+"02_trimmed/{sample}"+RCODE+"1_val_1.fq.gz", #---- this ALWAYS outputs .fq.qz format.
        PATHOUT+"02_trimmed/{sample}"+RCODE+"2_val_2.fq.gz",
        #        PATHOUT+"02_trimmed/{sample}"+RCODE+"1"+INEXT+"_trimming_report.txt",
        #        PATHOUT+"02_trimmed/{sample}"+RCODE+"2"+INEXT+"_trimming_report.txt"
    params:
        outdir         = "--output_dir "+PATHOUT+"02_trimmed/",
        phred          = "--phred33",
        gz             = "--gzip",
        cutadapt       = "--path_to_cutadapt "+CUTADAPT,
        FivePrimeClip  = "--clip_R1 10",
        ThreePrimeClip = "--three_prime_clip_R1 10",
        paired         = "--paired"
    log:
        PATHOUT+"02_trimmed/{sample}.trimgalore.log"
    message:
        " ---------  Trimming raw paired-end read data using {TRIMGALORE} -------  "
    shell:
        "nice -"+str(NICE)+" {TRIMGALORE} {params} {input} 2> {log}"

# ==========================================================================================
# raw quality control 

rule fastqc_raw: #----only need one: covers BOTH PE and SE cases.
    input:
        PATHIN+"{sample}"+INEXT
    output:
        PATHOUT+"01_rawqc/{sample}_fastqc.html",
        PATHOUT+"01_rawqc/{sample}_fastqc.zip"
    params:
        outdir = "--outdir "+PATHOUT+"01_rawqc/"     # usually pass params as strings instead of wildcards.

    log:
        PATHOUT+"01_rawqc/{sample}_fastqc.log"
    message: """ ----------  Quality checking raw read data with {FASTQC}.  --------------   """
    shell:
        "nice -"+str(NICE)+" {FASTQC} {params.outdir}  {input} 2> {log}"
