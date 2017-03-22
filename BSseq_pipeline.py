#!/usr/bin/env python3.5
# ---last updated on  Wed Mar 15 14:35:28 CET 2017  by  bosberg  at location  , sl-bimsb-p-ap1

#  changes from  Mon Mar 13 12:25:28 CET 2017 : manually set the bismark_se_ output to .bam with the "--bam" option and removed commented-out sections of the total [OUTPUT_FILE] list to make it cleaner

#  changes from  Mon Mar 6 19:20:15 CET 2017 : made the folder assignments per rule more systematic

#============================================================================================================
# SNAKEMAKE FILE WRITTEN BY THE AKALIN GROUP AT MDC, BERLIN, 2017
# Alexander Gosdschan, Katarzyna Wreczycka, Bren Osberg
# To process bisulfite sequencing data from raw fastq files to performing integrated bioinformatics analysis.

# SUBMIT THIS JOB INTERACTIVELY WITH:
# $nohup  snakemake -s [this filename] --jobs [# of jobs to submit] > [logfilename] &
# You can also add the following options for cluster submission: --cluster "qsub -V -l h_vmem={params.mem} -pe smp {params.threads} -l h_rt=36:00:00" &

# import IPython;
# IPython.embed()
 
#
#============================================================================================================

import os
#------ set config file and include function definitions:
configfile: "./config.json"
include   : "./func_defs.py"


#---------------------------------     DEFINE PATHS AND FILE NAMES:  ----------------------------------

RCODE           = config["RCODE"]                   #--- the string that denotes which read the file corresponds to (only relevent for paired-end experiments.)
PATHIN          = config["paths"]["PATHIN"]         #--- location of the data files to be imported
PATHOUT         = config["paths"]["PATHOUT"]        #--- where to send the output
GTOOLBOX        = config["paths"]["GTOOLBOX"]       #--- where the programs are stored to carry out the necessary operations
GENOMEPATH      = config["paths"]["GENOMEPATH"]     #--- where the reference genome being mapped to is stored
LOGS            = config["paths"]["LOGS"]           #--- subfolder name for the logs of some programs

INEXT           = config["INEXT"]                   #--- input file extension; usually .fq.gz, but can also be .bz2 among other possibilities.
VERSION         = config["files"]["VERSION"]        #--- version of the genome being mapped to.

CHROM_INFO      = config["files"]["CHROM_INFO"]     #--- details of the reference genome (length, etc.) haploid chroms have been removed.


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

OUTPUT_FILES = [

                #               ======  rule 01 raw QC    =========
                [ expand (list_files(PATHOUT+"01_rawqc/", config["SAMPLES"][sampleID]["files"], "_fastqc.html")  ) for sampleID in config["SAMPLES"]  ],
                [ expand (list_files(PATHOUT+"01_rawqc/", config["SAMPLES"][sampleID]["files"], "_fastqc.zip" )  ) for sampleID in config["SAMPLES"]  ],
                #               ======  rule 03 posttrim_QC_ ======
                [ expand ( list_files_posttrim_QC(PATHOUT+"03_posttrim_QC/", config["SAMPLES"][sampleID]["files"],".html")  ) for sampleID in config["SAMPLES"]  ],
                [ expand ( list_files_posttrim_QC(PATHOUT+"03_posttrim_QC/", config["SAMPLES"][sampleID]["files"],".html")  ) for sampleID in config["SAMPLES"]  ],
                #--- fastQC output files are not needed downstream and need to be called explicitly.

                
                #-----  here is a list of intermediary files to be uncommented back into execution if you want to stop part-way along the process -----------------
                #               ====rule 02 trimgalore ======
                #               [ expand ( list_files_TG( PATHOUT+"02_trimmed/", config["SAMPLES"][sampleID]["files"] )  ) for sampleID in config["SAMPLES"]  ],
                #               [ expand ( list_files_TG( PATHOUT+"02_trimmed/", config["SAMPLES"][sampleID]["files"] )  ) for sampleID in config["SAMPLES"]  ],                
                #               ====rule 04 Mapping ======
                #               [ expand ( list_files_bismark(PATHOUT+"04_mapped/", config["SAMPLES"][sampleID]["files"] )  ) for sampleID in config["SAMPLES"]  ],
                #               ====rule 05 Deduplication ======
                #               [ expand ( list_files_Dedupe(PATHOUT+"05_Deduped/", config["SAMPLES"][sampleID]["files"] )  ) for sampleID in config["SAMPLES"]  ],                                
                #               ====rule 06 extract_methylation ======           
                #               [ expand   ( list_files_xmeth( PATHOUT+"06_xmeth/", config["SAMPLES"][sampleID]["files"] )  ) for sampleID in config["SAMPLES"]  ]  
                                            
                # ==================  FINAL REPORT =========================
                [ expand (PATHOUT+config["SAMPLES"][sampleID]["files"][0]+SEPEstr(config["SAMPLES"][sampleID]["files"] )+"_report.html"  ) for sampleID in config["SAMPLES"]  ],
                
]

#--- In case you want to debug the code with interactive commands:
# import IPython;
# IPython.embed()
# print(" \n output files = \n")
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

rule bismark_se_report:
    input:
        aln   = PATHOUT+"04_mapped/{sample}_trimmed_bismark_bt2_SE_report.txt",
        sp    = PATHOUT+"06_xmeth/{sample}_trimmed_bismark_bt2.deduplicated_splitting_report.txt",
        dd    = PATHOUT+"05_Deduped/{sample}_trimmed_bismark_bt2.deduplication_report.txt",
        mbias = PATHOUT+"06_xmeth/{sample}_trimmed_bismark_bt2.deduplicated.M-bias.txt",
        nuc   = PATHOUT+"04_mapped/{sample}_trimmed_bismark_bt2.nucleotide_stats.txt"
    output:
        PATHOUT+"{sample}_trimmed_bismark_bt2_SE_report.html",
    params:
        dir = "--dir "+PATHOUT
    log: 
        PATHOUT+LOGS+"{sample}_SE_final_report.log"
    message: """----------- Generating Bismark report ----------- (with the following command:) \n  """
    shell:
        " {BISMARK2REPORT} {params} \
            --alignment_report {input.aln} \
            --splitting_report {input.sp} \
            --dedup_report {input.dd} \
            --mbias_report {input.mbias} \
            --nucleotide_report {input.nuc} \
            2> {log} "
#----------
rule bismark_pe_report:
    input:
        aln   = PATHOUT+"04_mapped/{sample}"+RCODE+"1_val_1_bismark_bt2_PE_report.txt",
        sp    = PATHOUT+"06_xmeth/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated_splitting_report.txt",
        dd    = PATHOUT+"05_Deduped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplication_report.txt",
        mbias = PATHOUT+"06_xmeth/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated.M-bias.txt",
        nuc   = PATHOUT+"04_mapped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.nucleotide_stats.txt"
    output:
        PATHOUT+"{sample}"+RCODE+"1_val_1_bismark_bt2_PE_report.html"
    params:
        dir = "--dir "+PATHOUT
    log: 
        PATHOUT+LOGS+"{sample}_PE_final_report.log"
    message: """----------- Generating Bismark report ----------- \n  """
    shell:
        " {BISMARK2REPORT} {params} \
        --alignment_report {input.aln} \
        --splitting_report {input.sp} \
        --dedup_report {input.dd} \
        --mbias_report {input.mbias} \
        --nucleotide_report {input.nuc} \
        2> {log} "

# ==========================================================================================
rule bismark_se_methylation_extractor:
    input:
        PATHOUT+"05_Deduped/{sample}_trimmed_bismark_bt2.deduplicated.bam"
    output:
        expand(PATHOUT+"06_xmeth/{{sample}}_trimmed_bismark_bt2.deduplicated.{file}.gz",  file=["bedGraph","bismark.cov","CpG_report.txt"]),
        PATHOUT+"06_xmeth/{sample}_trimmed_bismark_bt2.deduplicated.M-bias.txt",
        PATHOUT+"06_xmeth/{sample}_trimmed_bismark_bt2.deduplicated.M-bias_R1.png",
        PATHOUT+"06_xmeth/{sample}_trimmed_bismark_bt2.deduplicated_splitting_report.txt"
	#      	expand(PATHOUT+"06_xmeth/{type}_{strand}_{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz",type=["CHG","CHH","CpG"],strand=["OT","OB","CTOT","CTOB"]),
	#      	expand(PATHOUT+"06_xmeth/{type}_{strand}_{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz",type=["CHG","CHH","CpG"],strand=["OT","OB"]),
	#----- MANUALLY EXPANDED HERE : -----------
	# PATHOUT+"06_xmeth/CHG_OT_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",       
	# PATHOUT+"06_xmeth/CHG_OB_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",       
	# PATHOUT+"06_xmeth/CHH_OT_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",       
	# PATHOUT+"06_xmeth/CHH_OB_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",       
	# PATHOUT+"06_xmeth/CpG_OT_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",       
	# PATHOUT+"06_xmeth/CpG_OB_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",       
	#-----  AND WITH THE CTOT/CTOB FILES: ----
	#	PATHOUT+"06_xmeth/CpG_CTOT_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",       
	#	PATHOUT+"06_xmeth/CpG_CTOB_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",       
	#	PATHOUT+"06_xmeth/CHH_CTOT_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",       
	#	PATHOUT+"06_xmeth/CHH_CTOB_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",       
	#	PATHOUT+"06_xmeth/CHG_CTOT_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",       
	#	PATHOUT+"06_xmeth/CHG_CTOB_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",       
	#-----------------------------------------

    threads: 4
    params:
        se = "--single-end",
        gz = "--gzip",
        cReport = "--cytosine_report",
        bg = "--bedgraph",

        genomeFolder = "--genome_folder " + GENOMEPATH,
        outdir = "--output "+PATHOUT+"06_xmeth/"
    log: PATHOUT+"06_xmeth/{sample}_bismark_methylation_extraction.log"
    message: """--------------  Extracting  Methylation Information --------------- \n"""
    shell:
        "{BISMARK_METHYLATION_EXTRACTOR} {params} --multicore {threads} {input} 2> {log}"
#-----------------
rule bismark_pe_methylation_extractor:
    input:
        PATHOUT+"05_Deduped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated.bam"
    output:
        expand(PATHOUT+"06_xmeth/{{sample}}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated.{file}.gz",  file=["bedGraph","bismark.cov","CpG_report.txt"]),
        PATHOUT+"06_xmeth/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated.M-bias.txt",
        PATHOUT+"06_xmeth/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated.M-bias_R1.png",
        PATHOUT+"06_xmeth/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated_splitting_report.txt"    
    threads: 4
    params:
        pe = "--paired-end",
        gz = "--gzip",
        cReport = "--cytosine_report",
        bg = "--bedgraph",
        
        genomeFolder = "--genome_folder " + GENOMEPATH,
        outdir = "--output "+PATHOUT+"06_xmeth/"
    log: PATHOUT+"06_xmeth/{sample}_bismark_methylation_extraction.log"
    message: """--------------  Extracting  Methylation Information --------------- \n"""
    shell:
        "{BISMARK_METHYLATION_EXTRACTOR} {params} --multicore {threads} {input} 2> {log}"

# ==========================================================================================

rule bismark_se_deduplication:
    input:
        PATHOUT+"04_mapped/{sample}_trimmed_bismark_bt2.bam"
    output:
        PATHOUT+"05_Deduped/{sample}_trimmed_bismark_bt2.deduplicated.bam",
        PATHOUT+"05_Deduped/{sample}_trimmed_bismark_bt2.deduplication_report.txt"
    params:
        bam="--bam ",
        sampath="--samtools_path "+SAMTOOLS
    log:
        PATHOUT+"04_mapped/{sample}_deduplication.log"
    message: """-----------   Deduplicating read alignments ---------------------- """
    shell:
        """{DEDUPLICATE_BISMARK} {params} {input} 2> {log} ; 
        mv {PATHOUT}04_mapped/*dedupl* {PATHOUT}05_Deduped/"""

rule bismark_pe_deduplication:
    input:
        PATHOUT+"04_mapped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.bam"
    output:
        PATHOUT+"05_Deduped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated.bam",
        PATHOUT+"05_Deduped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplication_report.txt"
    params:
        bam="--bam ",
        sampath="--samtools_path "+SAMTOOLS,
        paired="--paired "
    log:
        PATHOUT+"04_mapped/{sample}_deduplication.log"
    message: """-----------   Deduplicating read alignments ---------------------- """
    shell:
        """{DEDUPLICATE_BISMARK} {params} {input} 2> {log} ; 
            mv {PATHOUT}04_mapped/*dedupl* {PATHOUT}05_Deduped/"""

# ==========================================================================================

rule bismark_se_nondirectional:
    input:
       PATHOUT+"02_trimmed/{sample}_trimmed.fq.gz"
    output:
        PATHOUT+"04_mapped/{sample}_trimmed_bismark_bt2.bam",
        PATHOUT+"04_mapped/{sample}_trimmed_bismark_bt2.nucleotide_stats.txt",
        PATHOUT+"04_mapped/{sample}_trimmed_bismark_bt2_SE_report.txt"
    threads: 4
    params:
        N = "-N 1",
        L = "-L 20",
        genomeFolder = "--genome_folder " + GENOMEPATH,
        outdir = "--output_dir  "+PATHOUT+"04_mapped/",
        nucCov = "--nucleotide_coverage",
        nonDir = "--non_directional ",
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(BOWTIE2) ,
        useBowtie2  = "--bowtie2 ",
        samtools    = "--samtools_path "+ os.path.dirname(SAMTOOLS),
        tempdir     = "--temp_dir "+PATHOUT
    log:
        PATHOUT+"04_mapped/{sample}_bismark_se_mapping.log"
    message: """-------------   Mapping reads to genome {VERSION}. ------------- """
    shell:
        "{BISMARK} {params} --multicore {threads} {input} 2> {log}"

#--------

rule bismark_pe_nondirectional:
    input:
        fin1 = PATHOUT+"02_trimmed/{sample}"+RCODE+"1_val_1.fq.gz",
        fin2 = PATHOUT+"02_trimmed/{sample}"+RCODE+"2_val_2.fq.gz"
    output:
        PATHOUT+"04_mapped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.bam",
        PATHOUT+"04_mapped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.nucleotide_stats.txt",
        PATHOUT+"04_mapped/{sample}"+RCODE+"1_val_1_bismark_bt2_PE_report.txt"
    threads: 4
    params:
        N = "-N 1",
        L = "-L 20",
        genomeFolder = "--genome_folder " + GENOMEPATH,
        outdir = "--output_dir  "+PATHOUT+"04_mapped/",
        nucCov = "--nucleotide_coverage",
        nonDir = "--non_directional ",
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(BOWTIE2) ,
        useBowtie2  = "--bowtie2 ",
        samtools    = "--samtools_path "+ os.path.dirname(SAMTOOLS),
        tempdir     = "--temp_dir "+PATHOUT
    log:
        PATHOUT+"04_mapped/{sample}_bismark_pe_mapping.log"
    message: """-------------   Mapping reads to genome {VERSION}. ------------- """
    shell:
        "{BISMARK} {params} --multicore {threads} -1 {input.fin1} -2 {input.fin2} 2> {log}"

# ==========================================================================================
#### ----  THIS ONLY GETS INVOKED WHEN MANUALLY CALLED SPECIFICIALLY ------

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
        "{BISMARK_GENOME_PREPARATION} {params} {input} 2> {log}"

# ==========================================================================================

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
    message: """ ------------  Quality checkin trimmmed data with Fastqc ------------- """
    shell:
      	"{FASTQC} {params.outdir} {input} 2> {log}"

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
    message: """ ------------  Quality checkin trimmmed data with Fastqc ------------- """
    shell:
      	"{FASTQC} {params.outdir} {input} 2> {log}"
#
# ==========================================================================================

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
       "{TRIMGALORE} {params} {input} 2> {log}"

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
        " ---------  Trimming raw paired read data using {TRIMGALORE} -------  "
    shell:
        "{TRIMGALORE} {params} {input} 2> {log}"

# ==========================================================================================

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
        "{FASTQC} {params.outdir}  {input} 2> {log}"
