#!/usr/bin/env python3.5
# ---last updated on  Wed Mar 29 16:05:58 CEST 2017  by  blosberg  at location  , Bren-Osbergs-MacBook.local

#  changes from  Wed Mar 29 16:05:58 CEST 2017 : mapping alignment rule can now be directional or non-directional

#  changes from  Thu Mar 23 15:44:48 CET 2017 : Eliminated race condition (from file latency) in dedupe rule by no-longer moving output to separate subfolder. decremented folder names accordingly.

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

OUTPUT_FILES = [

                #               ======  rule 01 raw QC    =========
                #3[ expand (list_files(PATHOUT+"01_rawqc/", config["SAMPLES"][sampleID]["files"], "_fastqc.html")  ) for sampleID in config["SAMPLES"]  ],
                # [ expand (list_files(PATHOUT+"01_rawqc/", config["SAMPLES"][sampleID]["files"], "_fastqc.zip" )  ) for sampleID in config["SAMPLES"]  ],
                #               ======  rule 03 posttrim_QC_ ======
                # [ expand ( list_files_posttrim_QC(PATHOUT+"03_posttrim_QC/", config["SAMPLES"][sampleID]["files"],".html")  ) for sampleID in config["SAMPLES"]  ],
                # [ expand ( list_files_posttrim_QC(PATHOUT+"03_posttrim_QC/", config["SAMPLES"][sampleID]["files"],".html")  ) for sampleID in config["SAMPLES"]  ],
                #--- fastQC output files are not needed downstream and need to be called explicitly.

                
                #-----  here is a list of intermediary files to be uncommented back into execution if you want to stop part-way along the process -----------------
                #               ====rule 02 trimgalore ======
                #               [ expand ( list_files_TG( PATHOUT+"02_trimmed/", config["SAMPLES"][sampleID]["files"] )  ) for sampleID in config["SAMPLES"]  ],
                #               [ expand ( list_files_TG( PATHOUT+"02_trimmed/", config["SAMPLES"][sampleID]["files"] )  ) for sampleID in config["SAMPLES"]  ],                
                #               ====rule 04 Mapping ======
                #               [ expand ( list_files_bismark(PATHOUT+"04_mapped_n_deduped/", config["SAMPLES"][sampleID]["files"] )  ) for sampleID in config["SAMPLES"]  ],

                #               ====rule 04a sorting ======
                [ expand ( list_files_sortbam(PATHOUT+"04a_sorted/", config["SAMPLES"][sampleID]["files"] )  ) for sampleID in config["SAMPLES"]  ],

                #               ====rule 05 Deduplication ======
                #               [ expand ( list_files_dedupe(PATHOUT+"04_mapped_n_deduped/", config["SAMPLES"][sampleID]["files"] )  ) for sampleID in config["SAMPLES"]  ],                                
                #               ====rule 06 extract_methylation ======           
                #               [ expand   ( list_files_xmeth( PATHOUT+"05_xmeth/", config["SAMPLES"][sampleID]["files"] )  ) for sampleID in config["SAMPLES"]  ]  
                                            
                # ==================  FINAL REPORT =========================
                # [ expand (PATHOUT+config["SAMPLES"][sampleID]["files"][0]+SEPEstr(config["SAMPLES"][sampleID]["files"] )+"_report.html"  ) for sampleID in config["SAMPLES"]  ],
                
]

#--- In case you want to debug the code with interactive commands:
# import IPython;
# IPython.embed()
# print(" \n debugging strings here \n")
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
        aln   = PATHOUT+"04_mapped_n_deduped/{sample}_trimmed_bismark_bt2_SE_report.txt",
        sp    = PATHOUT+"05_xmeth/{sample}_trimmed_bismark_bt2.deduplicated_splitting_report.txt",
        dd    = PATHOUT+"04_mapped_n_deduped/{sample}_trimmed_bismark_bt2.deduplication_report.txt",
        mbias = PATHOUT+"05_xmeth/{sample}_trimmed_bismark_bt2.deduplicated.M-bias.txt",
        nuc   = PATHOUT+"04_mapped_n_deduped/{sample}_trimmed_bismark_bt2.nucleotide_stats.txt"
    output:
        PATHOUT+"{sample}_trimmed_bismark_bt2_SE_report.html",
    params:
        dir = "--dir "+PATHOUT
    log: 
        PATHOUT+LOGS+"{sample}_SE_final_report.log"
    message: """----------- Generating single-end Bismark report ----------- (with the following command:) \n  """
    shell:
        "nice -"+str(NICE)+" {BISMARK2REPORT} {params} \
            --alignment_report {input.aln} \
            --splitting_report {input.sp} \
            --dedup_report {input.dd} \
            --mbias_report {input.mbias} \
            --nucleotide_report {input.nuc} \
            2> {log} "
#----------
rule bismark_pe_report:
    input:
        aln   = PATHOUT+"04_mapped_n_deduped/{sample}"+RCODE+"1_val_1_bismark_bt2_PE_report.txt",
        sp    = PATHOUT+"05_xmeth/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated_splitting_report.txt",
        dd    = PATHOUT+"04_mapped_n_deduped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplication_report.txt",
        mbias = PATHOUT+"05_xmeth/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated.M-bias.txt",
        nuc   = PATHOUT+"04_mapped_n_deduped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.nucleotide_stats.txt"
    output:
        PATHOUT+"{sample}"+RCODE+"1_val_1_bismark_bt2_PE_report.html"
    params:
        dir = "--dir "+PATHOUT
    log: 
        PATHOUT+LOGS+"{sample}_PE_final_report.log"
    message: """----------- Generating paired-end Bismark report ----------- \n  """
    shell:
        "nice -"+str(NICE)+" {BISMARK2REPORT} {params} \
        --alignment_report {input.aln} \
        --splitting_report {input.sp} \
        --dedup_report {input.dd} \
        --mbias_report {input.mbias} \
        --nucleotide_report {input.nuc} \
        2> {log} "

# ==========================================================================================
rule bismark_se_methylation_extractor:
    input:
        PATHOUT+"04_mapped_n_deduped/{sample}_trimmed_bismark_bt2.deduplicated.bam"
    output:
        expand(PATHOUT+"05_xmeth/{{sample}}_trimmed_bismark_bt2.deduplicated.{file}.gz",  file=["bedGraph","bismark.cov","CpG_report.txt"]),
        PATHOUT+"05_xmeth/{sample}_trimmed_bismark_bt2.deduplicated.M-bias.txt",
        PATHOUT+"05_xmeth/{sample}_trimmed_bismark_bt2.deduplicated_splitting_report.txt"
	#      	expand(PATHOUT+"05_xmeth/{type}_{strand}_{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz",type=["CHG","CHH","CpG"],strand=["OT","OB","CTOT","CTOB"]),
	#      	expand(PATHOUT+"05_xmeth/{type}_{strand}_{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz",type=["CHG","CHH","CpG"],strand=["OT","OB"]),
	#----- MANUALLY EXPANDED HERE : -----------
	# PATHOUT+"05_xmeth/CHG_OT_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",
	# PATHOUT+"05_xmeth/CHG_OB_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",
	# PATHOUT+"05_xmeth/CHH_OT_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",
	# PATHOUT+"05_xmeth/CHH_OB_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",
	# PATHOUT+"05_xmeth/CpG_OT_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",
	# PATHOUT+"05_xmeth/CpG_OB_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",
	#-----  AND WITH THE CTOT/CTOB FILES: ----
	#	PATHOUT+"05_xmeth/CpG_CTOT_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",
	#	PATHOUT+"05_xmeth/CpG_CTOB_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",
	#	PATHOUT+"05_xmeth/CHH_CTOT_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",
	#	PATHOUT+"05_xmeth/CHH_CTOB_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",
	#	PATHOUT+"05_xmeth/CHG_CTOT_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",
	#	PATHOUT+"05_xmeth/CHG_CTOB_{sample}_trimmed_bismark_bt2.deduplicated.txt.gz",
	#-----------------------------------------

    threads: 2 
    params:
        se = "--single-end ",
        gz = "--gzip ",
        cReport = "--cytosine_report ",
        bg = "--bedgraph ",
        genomeFolder = "--genome_folder " + GENOMEPATH,
        outdir = "--output "+PATHOUT+"05_xmeth/"
    log: PATHOUT+"05_xmeth/{sample}_bismark_methylation_extraction.log"
    message: """--------------  Extracting  single-end Methylation Information --------------- \n"""
    shell:
        "nice -"+str(NICE)+" {BISMARK_METHYLATION_EXTRACTOR} {params} --multicore {threads} {input} 2> {log}"
#-----------------
rule bismark_pe_methylation_extractor:
    input:
        PATHOUT+"04_mapped_n_deduped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated.bam"
    output:
        expand(PATHOUT+"05_xmeth/{{sample}}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated.{file}.gz",  file=["bedGraph","bismark.cov","CpG_report.txt"]),
        PATHOUT+"05_xmeth/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated.M-bias.txt",
        PATHOUT+"05_xmeth/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated_splitting_report.txt"
    threads: 2
    params:
        pe = "--paired-end",
        gz = "--gzip",
        cReport = "--cytosine_report",
        bg = "--bedgraph",
        
        genomeFolder = "--genome_folder " + GENOMEPATH,
        outdir = "--output "+PATHOUT+"05_xmeth/"
    log: PATHOUT+"05_xmeth/{sample}_bismark_methylation_extraction.log"
    message: """--------------  Extracting  paired-end Methylation Information --------------- \n"""
    shell:
        "nice -"+str(NICE)+" {BISMARK_METHYLATION_EXTRACTOR} {params} --multicore {threads} {input} 2> {log}"
# ==========================================================================================
rule sortbam_se:
    input:
        PATHOUT+"04_mapped_n_deduped/{sample}_trimmed_bismark_bt2.deduplicated.bam"
    output:
        PATHOUT+"04a_sorted/{sample}_trimmed_bismark_bt2.deduplicated.sorted.bam"
    shell:
        "nice -"+str(NICE)+" samtools sort {input} -o {output}"

rule sortbam_pe:
    input:
        PATHOUT+"04_mapped_n_deduped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated.bam"
    output:
        PATHOUT+"04a_sorted/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated.sorted.bam",
    shell:
        "nice -"+str(NICE)+" samtools sort {input} -o {output}"
# ==========================================================================================

rule bismark_se_deduplication:
    input:
        PATHOUT+"04_mapped_n_deduped/{sample}_trimmed_bismark_bt2.bam"
    output:
        PATHOUT+"04_mapped_n_deduped/{sample}_trimmed_bismark_bt2.deduplicated.bam",
        PATHOUT+"04_mapped_n_deduped/{sample}_trimmed_bismark_bt2.deduplication_report.txt"
    params:
        bam="--bam ",
        sampath="--samtools_path "+SAMTOOLS
    log:
        PATHOUT+"04_mapped_n_deduped/{sample}_deduplication.log"
    message: """-----------   Deduplicating single-end read alignments ---------------------- """
    shell:
        """nice -"+str(NICE)+" {DEDUPLICATE_BISMARK} {params} {input} 2> {log} """

rule bismark_pe_deduplication:
    input:
        PATHOUT+"04_mapped_n_deduped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.bam"
    output:
        PATHOUT+"04_mapped_n_deduped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplicated.bam",
        PATHOUT+"04_mapped_n_deduped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.deduplication_report.txt"
    params:
        bam="--bam ",
        sampath="--samtools_path "+SAMTOOLS,
        paired="--paired "
    log:
        PATHOUT+"04_mapped_n_deduped/{sample}_deduplication.log"
    message: """-----------   Deduplicating paired-end read alignments ---------------------- """
    shell:
        """nice -"+str(NICE)+" {DEDUPLICATE_BISMARK} {params} {input} 2> {log} """

# ==========================================================================================

rule bismark_se:
    input:
       PATHOUT+"02_trimmed/{sample}_trimmed.fq.gz"
    output:
        PATHOUT+"04_mapped_n_deduped/{sample}_trimmed_bismark_bt2.bam",
        PATHOUT+"04_mapped_n_deduped/{sample}_trimmed_bismark_bt2.nucleotide_stats.txt",
        PATHOUT+"04_mapped_n_deduped/{sample}_trimmed_bismark_bt2_SE_report.txt"
    threads: 2
    params:
        N = "-N 1",
        L = "-L 20",
        genomeFolder = "--genome_folder " + GENOMEPATH,
        outdir = "--output_dir  "+PATHOUT+"04_mapped_n_deduped/",
        nucCov = "--nucleotide_coverage",
	    nonDir = NON_DIR_FLAG,             #--- THIS IS EMPTY IF DATA IS DIRECTIONAL (DEFAULT) OTHERWISE "non_directional". 
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(BOWTIE2) ,
        useBowtie2  = "--bowtie2 ",
        samtools    = "--samtools_path "+ os.path.dirname(SAMTOOLS),
        tempdir     = "--temp_dir "+PATHOUT
    log:
        PATHOUT+"04_mapped_n_deduped/{sample}_bismark_se_mapping.log"
    message: """-------------   Mapping single-end reads to genome {VERSION}. ------------- """
    shell:
        "nice -"+str(NICE)+" {BISMARK} {params} --multicore {threads} {input} 2> {log}"

#--------

rule bismark_pe:
    input:
        fin1 = PATHOUT+"02_trimmed/{sample}"+RCODE+"1_val_1.fq.gz",
        fin2 = PATHOUT+"02_trimmed/{sample}"+RCODE+"2_val_2.fq.gz"
    output:
        PATHOUT+"04_mapped_n_deduped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.bam",
        PATHOUT+"04_mapped_n_deduped/{sample}"+RCODE+"1_val_1_bismark_bt2_pe.nucleotide_stats.txt",
        PATHOUT+"04_mapped_n_deduped/{sample}"+RCODE+"1_val_1_bismark_bt2_PE_report.txt"
    threads: 2
    params:
        N = "-N 1",
        L = "-L 20",
        genomeFolder = "--genome_folder " + GENOMEPATH,
        outdir = "--output_dir  "+PATHOUT+"04_mapped_n_deduped/",
        nucCov = "--nucleotide_coverage",
        nonDir = NON_DIR_FLAG,             #--- THIS IS EMPTY IF DATA IS DIRECTIONAL (DEFAULT) OTHERWISE "non_directional".
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(BOWTIE2) ,
        useBowtie2  = "--bowtie2 ",
        samtools    = "--samtools_path "+ os.path.dirname(SAMTOOLS),
        tempdir     = "--temp_dir "+PATHOUT
    log:
        PATHOUT+"04_mapped_n_deduped/{sample}_bismark_pe_mapping.log"
    message: """-------------   Mapping paired-end reads to genome {VERSION}. ------------- """
    shell:
        "nice -"+str(NICE)+" {BISMARK} {params} --multicore {threads} -1 {input.fin1} -2 {input.fin2} 2> {log}"

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
        "nice -"+str(NICE)+" {BISMARK_GENOME_PREPARATION} {params} {input} 2> {log}"

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
    message: """ ------------  Quality checking trimmmed single-end data with Fastqc ------------- """
    shell:
        "nice -"+str(NICE)+" {FASTQC} {params.outdir} {input} 2> {log}"

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
