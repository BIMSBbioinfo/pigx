#!/usr/bin/env python3.5

#============================================================================================================
# SNAKEMAKE FILE WRITTEN BY THE AKALIN GROUP AT MDC, BERLIN, 2017
# Alexander Gosdschan, Katarzyna Wreczycka, Bren Osberg
# To process bisulfite sequencing data from raw fastq files to performing integrated bioinformatics analysis.
#============================================================================================================


# This is how now it should be done
#snakemake  --forceall --snakefile /home/kwreczy/repositories/pigx_bsseq//Snakefile_test.py -configfile=./config.py --config tablesheet=/home/kwreczy/repositories/pigx_bsseq/test_dataset/TableSheet_test.csv

#------  include function definitions, and set os:

import os,csv,json
include   : "./scripts/func_defs.py"
include   : "./scripts/functions_parsing.py"


#------  Parse table sheet and write config file:

# Check configuration parameters  
tablesheetpath = parse_config_args(config)

# Read a Table Sheet in the CSV format
TABLESHEET =  config["tablesheet"]
# TODO: for now instead of reading a big table I splitted it into few smaller ones
# TODO: split this table to separate files
# save it in eg tmp and then remove it
#f = open(TABLESHEET,'r')
#text = f.readlines()
general_paramspath = open("/home/kwreczy/repositories/pigx_bsseq/test_dataset/GENERALPARAMETERS.txt",'r').readlines() # TODO:
pathtable = "/home/kwreczy/repositories/pigx_bsseq/test_dataset/SAMPLES.txt" # TODO:
progspath= "/home/kwreczy/repositories/pigx_bsseq/test_dataset/PROGS.json" # TODO:

# Load general parameters
gen_params=parseGeneralParams2dict(general_paramspath)
# Load parameters specific to samples
sample_params = parseTable2dict( pathtable, skip=1)
# Load paths to tools
with open(progspath) as data_file:    
    progs=json.load(data_file)

# Create a config file  
config=gen_params
config.update(sample_params)
config.update(progs)

# Check if given directories exist, if they dont create them
for x in [config["PATHIN"], config["PATHOUT"], config["LOG"]]:
  if not os.path.exists(x):
    os.makedirs(x)

# Save the config file
with open(config["PATHOUT"]+"config.json", 'w') as outfile:
 json.dump(config, outfile)

print(config)
print(config["PATHOUT"]+"config.json")

#---------------------------------     DEFINE PATHS AND FILE NAMES:  ----------------------------------

PATHIN          = config["PATHIN"]         #--- location of the data files to be imported
PATHOUT         = config["PATHOUT"]        #--- where to send the output
GTOOLBOX        = config["GTOOLBOX"]       #--- where the programs are stored to carry out the necessary operations
GENOMEPATH      = config["GENOMEPATH"]     #--- where the reference genome being mapped to is stored

LOGS            = config["LOG"]           #--- subfolder name for the logs of some programs

VERSION         = config["GENOME_VERSION"]  #--- version of the genome being mapped to.

CHROM_INFO      = config["CHROM_INFO"]     #--- details of the reference genome (length, etc.) haploid chroms have been removed.
NUMTHREADS      = config["NUMTHREADS"]



#-------------------------------      DEFINE PROGRAMS TO BE EXECUTED: ---------------------------------

FASTQC                         =  GTOOLBOX+config["PROGS"]["FASTQC"]            #--- self-explanatory program names.
TRIMGALORE                     =  GTOOLBOX+config["PROGS"]["TRIMGALORE"]
CUTADAPT                       =  GTOOLBOX+config["PROGS"]["CUTADAPT"]
BISMARK_GENOME_PREPARATION     =  GTOOLBOX+config["PROGS"]["BISMARK_GENOME_PREPARATION"]
BISMARK                        =  GTOOLBOX+config["PROGS"]["BISMARK"]
BOWTIE2                        =  GTOOLBOX+config["PROGS"]["BOWTIE2"]
DEDUPLICATE_BISMARK            =  GTOOLBOX+config["PROGS"]["DEDUPLICATE_BISMARK"]
BISMARK_METHYLATION_EXTRACTOR  =  GTOOLBOX+config["PROGS"]["BISMARK_METHYLATION_EXTRACTOR"]
BISMARK2REPORT                 =  GTOOLBOX+config["PROGS"]["BISMARK2REPORT"]

SAMTOOLS                       =  GTOOLBOX+config["PROGS"]["SAMTOOLS"]


#---------------------------     LIST THE OUTPUT DIRECTORIED AND SUBDIRECTORIED TO BE PRODUCED     ------------------------------

DIR_sorted=PATHOUT+'06_sorted/'
DIR_mapped_n_deduped=PATHOUT+'04_mapped_n_deduped/'
DIR_posttrim_QC=PATHOUT+'03_posttrim_QC/'
DIR_trimmed=PATHOUT+'02_trimmed/'
DIR_rawqc=PATHOUT+'01_rawqc/'

for x in [DIR_rawqc, DIR_trimmed, DIR_posttrim_QC, DIR_mapped_n_deduped, DIR_sorted]:
  if not os.path.exists(x):
    os.makedirs(x)


# #---------------------------     LIST THE OUTPUT FILES TO BE PRODUCED     ------------------------------



## OUTPUTFILE contain files named nnot by fastq files but for sampleid!!!! (besides fastqc.)

OUTPUT_FILES = [
      #               ======  rule 01 raw QC    =========
      [ expand (list_files(DIR_rawqc, getFilenames(config["SAMPLES"][x]['fastq']), "_fastqc.html")  ) for x in config["SAMPLES"].keys()  ],
      #----RULE 2 IS ALWAYS EXECUTED, TRIMMING IS A PREREQUISITE FOR SUBSEQUENT RULES ----
      #               ======  rule 03 posttrim_QC_ ======
      #[ expand ( list_files_posttrim_QC(DIR_posttrim_QC, config["units"][x],".html")  ) for x in config["samples"].keys()  ],
      #--- fastQC output files are not needed downstream and need to be called explicitly.
      #               ====rule 02 trimgalore ======
      #[ expand ( list_files_TG( DIR_trimmed, config["units"][x])  ) for x in config["samples"].keys()  ],
      [ expand (list_files_TG(DIR_trimmed, config["SAMPLES"][x]['fastq'], x  )) for x in config["SAMPLES"].keys()  ]
      #               ====rule 04 Mapping ======
      #[ expand ( list_files_bismark( DIR_mapped_n_deduped, config["units"][x] )  ) for x in config["samples"].keys()  ],
      #               ====rule 05 Deduplication ======
      #[ expand ( list_files_dedupe( DIR_mapped_n_deduped, config["units"][x], config["samples"][x], RCODE )  ) for x in config["samples"].keys()  ],
      #               ====rule 04a Sorting ======
      #[ expand (list_files_sortbam( DIR_sorted, config["units"][x], config["samples"][x] )  ) for x in config["samples"].keys()  ],
      #               ====rule 04a Indexing ======
      #[ expand (list_files_indexbam( DIR_sorted, config["units"][x], config["samples"][x] )  ) for x in config["samples"].keys()  ]

      ]

print('---------------OUTPUT_FILES----------------')
print(OUTPUT_FILES)


rule final:
  input:
    OUTPUT_FILES


# # ==========================================================================================
# 
# # #TODO: remove it. I created it only to force sortindex rule be ran be snakemake.
# # rule dummy_bai:
# #     input:
# #         DIR_sorted+"{sample}_trimmed_bismark_bt2.deduplicated.sorted.bam.bai"
# # 
# # 
# # rule sortindex:
# #     input:
# #         DIR_sorted+"{sample}_trimmed_bismark_bt2.deduplicated.sorted.bam"
# #     output:
# #         "{DIR_sorted}/{sample}_trimmed_bismark_bt2.deduplicated.sorted.bam.bai"
# #     threads: 2
# #     shell: 
# #         "{SAMTOOLS} index {input}"    
# # 
# # 
# # # ==========================================================================================
# # 
# # rule sortbam:
# #     input:
# #         DIR_mapped_n_deduped+"{sample}_bismark_bt2.deduplicated.bam"
# #     output:
# #         "{DIR_sorted}/{sample}_bismark_bt2.deduplicated.sorted.bam"
# #     threads: 2
# #     shell: 
# #         "{SAMTOOLS} sort --threads {threads} {input} -o {output}"    
# 
# 
# # ==========================================================================================
# 
# #--------
# # rule bismark_pe_deduplication:
# #     input:
# #         DIR_mapped_n_deduped+"{sample}"+RCODE[0] +reads[0]+"_bismark_bt2_pe.bam",
# #         DIR_mapped_n_deduped+"{sample}"+RCODE[0] +reads[0]+"_bismark_bt2_PE_report.txt",
# #         DIR_mapped_n_deduped+"{sample}"+RCODE[1] +reads[1]+"_bismark_bt2_pe.bam",
# #         DIR_mapped_n_deduped+"{sample}"+RCODE[1] +reads[1]+"_bismark_bt2_PE_report.txt"
# #     output:
# #         DIR_mapped_n_deduped+"{sample}"+RCODE[0] +reads[0]+"_bismark_bt2_pe.deduplicated.bam",
# #         DIR_mapped_n_deduped+"{sample}"+RCODE[1] +reads[1]+"_bismark_bt2_pe.deduplicated.bam", ##TODO: add report
# #     log:
# #         DIR_mapped_n_deduped+"{sample}_deduplication.log"
# #     message: """-----------   Deduplicating paired-end read alignments ---------------------- """
# #     shell:
# #         "{samtools} rmdup {input} {output} 2> {log}"
# #         
# #         
# # rule bismark_se_deduplication:
# #     input:
# #         DIR_mapped_n_deduped+"{sample}_trimmed_bismark_bt2.bam"
# #     output:
# #         "{DIR_mapped_n_deduped}/{sample}_trimmed_bismark_bt2.deduplicated.bam",
# #         "{DIR_mapped_n_deduped}/{sample}_trimmed_bismark_bt2.deduplication_report.txt"
# #     params:
# #         bam="--bam ",
# #         sampath="--samtools_path "+SAMTOOLS
# #     log:
# #         "{DIR_mapped_n_deduped}/{sample}_deduplication.log"
# #     message: """-----------   Deduplicating single-end read alignments ---------------------- """
# #     shell:
# #         "{DEDUPLICATE_BISMARK} {params} {input} 2> {log} "
# 
# 
# # ==========================================================================================
# # 
# # rule bismark_pe:
# #     input:
# #         fin1=DIR_trimmed+"{sample}"+RCODE[0]+"_val_1.fq.gz",
# #         fin2=DIR_trimmed+"{sample}"+RCODE[1]+"_val_2.fq.gz"
# #     output:
# #         DIR_mapped_n_deduped+"{sample}"+RCODE[0]+"_val_1_bismark_bt2_pe.bam",
# #         DIR_mapped_n_deduped+"{sample}"+RCODE[0]+"_val_1_bismark_bt2_PE_report.txt"
# #     threads: 2
# #     params:
# #         extra        = config.get("bismark_args", ""),
# #         useBowtie2   = "--bowtie2 ",
# #         outdir       = "--output_dir  "+ DIR_mapped_n_deduped,
# #         pathToBowtie = "--path_to_bowtie "+ os.path.dirname(BOWTIE2),
# #         samtools     = "--samtools_path "+ os.path.dirname(SAMTOOLS),
# #         tempdir      = "--temp_dir " + PATHOUT,
# #     log:
# #         DIR_mapped_n_deduped+"{sample}_bismark_pe_mapping.log"
# #     message: """-------------   Mapping paired-end reads to genome {VERSION}. ------------- """
# #     shell:
# #         "{BISMARK} {params} --multicore {threads} -1 {input.fin1} -2 {input.fin2} 2> {log}"
# # 
# # 
# # rule bismark_se:
# #     input:
# #         DIR_trimmed+"{sample}_trimmed.fq.gz"
# #     output:
# #         "{DIR_mapped_n_deduped}/{sample}.bismark_bt2.bam",
# #         "{DIR_mapped_n_deduped}/{sample}_trimmed_bismark_bt2_SE_report.txt"
# #     threads: 2
# #     params:
# #         extra       = config.get("bismark_args", ""),
# #         outdir      = "--output_dir " + DIR_mapped_n_deduped,
# #         tempdir     = "--temp_dir " + PATHOUT,
# #         sampath     ="--samtools_path " + os.path.dirname(SAMTOOLS),
# #         useBowtie2  = "--bowtie2 ",
# #         bowtie2path = "--path_to_bowtie " + os.path.dirname(BOWTIE2)
# #     log:
# #         "{DIR_mapped_n_deduped}/{sample}_bismark_se_mapping.log"
# #     message: """-------------   Mapping single-end reads to genome {VERSION}. ------------- """
# #     shell:
# #         "{BISMARK} {params} --multicore {threads} {input} 2> {log}"


#        unit.pe1_trimmed.bismark_bt2.bam
#        mv unit.pe1_trimmed.bismark_bt2.bam {sample}.bismark_bt2.bam
        
# # 
# # 
# # 
# # # ==========================================================================================
# # # post-trimming quality control
# # 
# # rule fastqc_after_trimming_se:
# #     input:
# #        DIR_trimmed+"{sample}_trimmed.fq.gz",
# #     output:
# #     	  DIR_posttrim_QC+"{sample}_trimmed_fastqc.html",
# #     	  DIR_posttrim_QC+"{sample}_trimmed_fastqc.zip"
# #     params:
# #        fastqc_args = config.get("fastqc_args", ""),
# #        outdir = "--outdir "+DIR_posttrim_QC
# #     log:
# #    	   DIR_posttrim_QC+"{sample}_trimmed_fastqc.log"
# #     message: """ ------------  Quality checking trimmmed single-end data with Fastqc ------------- """
# #     shell:
# #        "{FASTQC} {params.outdir} {input} 2> {log}"
# # #--------
# # rule fastqc_after_trimming_pe:
# #     input:
# #         DIR_trimmed+"{sample}"+RCODE[0]+"_val_1.fq.gz",
# #         DIR_trimmed+"{sample}"+RCODE[1]+"_val_2.fq.gz"
# #     output:
# #     	DIR_posttrim_QC+"{sample}"+RCODE[0]+"_val_1_fastqc.html",
# #     	DIR_posttrim_QC+"{sample}"+RCODE[0]+"_val_1_fastqc.zip",
# #     	DIR_posttrim_QC+"{sample}"+RCODE[1]+"_val_2_fastqc.zip",
# #         DIR_posttrim_QC+"{sample}"+RCODE[1]+"_val_2_fastqc.html"
# #     params:
# #         fastqc_args = config.get("fastqc_args", ""),
# #         outdir = "--outdir "+DIR_posttrim_QC
# #     log:
# #    	    DIR_posttrim_QC+"{sample}_trimmed_fastqc.log"
# #     message: """ ------------  Quality checking trimmmed paired-end data with Fastqc ------------- """
# #     shell:
# #         "{FASTQC} {params.outdir} {input} 2> {log}"
# #         
# # # #
# # # # ==========================================================================================
# # # # trim the reads
# # # 
# 
# rule trimgalore_se:
#    input:
#        PATHIN+"{sample}"+INEXT
#    output:
#        DIR_trimmed+"{sample}_trimmed.fq.gz" #---- this ALWAYS outputs .fq.qz format.
#    params:
#        extra        = config.get("trim_galore_args", ""),
#        outdir       = "--output_dir "+  DIR_trimmed,
#        phred        = "--phred33",
#        gz           = "--gzip",
#        cutadapt     = "--path_to_cutadapt "+ CUTADAPT,
#    log:
#        DIR_trimmed+"{sample}.trimgalore.log"
#    message:
#        " ---------  Trimming raw single-end read data using {TRIMGALORE} -------  "
#    shell:
#        "{TRIMGALORE} {params} {input} -o {DIR_trimmed} 2> {log}"
# 
# 
# rule trimgalore_pe:
#     input:
#         PATHIN + "{sample}" + RCODE[0] + INEXT,
#         PATHIN + "{sample}" + RCODE[1] + INEXT
#     output:
#         DIR_trimmed+"{sample}" + RCODE[0] + "_val_1.fq.gz", #---- this ALWAYS outputs .fq.qz format.
#         DIR_trimmed+"{sample}" + RCODE[1] + "_val_2.fq.gz",
#     params:
#         extra          = config.get("trim_galore_args", ""),
#         outdir         = "--output_dir " + DIR_trimmed,
#         gz             = "--gzip",
#         cutadapt       = "--path_to_cutadapt "+ CUTADAPT,
#         paired         = "--paired"
#     log:
#         DIR_trimmed+"{sample}.trimgalore.log"
#     message:
#         " ---------  Trimming raw paired-end read data using {TRIMGALORE} -------  "
#     shell:
#         "echo {output}; {TRIMGALORE} {params} {input} -o {DIR_trimmed} 2> {log}"
# 

# 
#                    
# # # a=[getFilenames(config["SAMPLES"][x]['fastq']) for x in config["SAMPLES"].keys()]
# # print("aaaaaaaaaaaaaaaaaaaaaaaaaa")
# # # print(a)
# # a=expand(DIR_trimmed+'{name}{sufftrimgalore}.{inext}',
# #                    name=['sampleB.pe1', 'sampleB.pe2'],
# #                    sufftrimgalore=['_val_1','_val_2'],
# #                    inext="fq.gz")
# # print(a)
# 
# 

def get_trimgalore_input(wc):
   samps = config['SAMPLES'][wc.sample]['fastq']

   if type(samps) is str:
        samps = [samps]

   if(len(samps)==2):
     return( [os.path.join(PATHIN, samps[0]), os.path.join(PATHIN, samps[1])] )
   else:
     return( [os.path.join(PATHIN, samps[0])] )


rule trimgalore_pe:
    input:
        infile = get_trimgalore_input
    output:
        output1 = DIR_trimmed+"{sample}_val_1.fq.gz",
        output2 = DIR_trimmed+"{sample}_val_2.fq.gz",
    params:
        extra          = config.get("trim_galore_args", ""),
        outdir         = "--output_dir " + DIR_trimmed,
        gz             = "--gzip",
        cutadapt       = "--path_to_cutadapt "+ CUTADAPT,
    log:
        log=DIR_trimmed+"{sample}.trimgalore.log"
    message:
        " ---------  Trimming raw paired-end read data using {TRIMGALORE} -------  "
    run:
        cmd = " ".join(
        [TRIMGALORE,
         " --paired",
         input.infile[0], input.infile[1],
         params.extra, params.outdir, params.gz, params.cutadapt,
         "-o "+ DIR_trimmed,
         " 2> ", log.log
         ])

        if len(input.infile) == 2:
          # Here a hack to generate files named like samples and not like units, even though
          # trimgalore will produce files named using units
          cmd = cmd + "; touch "+output.output1 + "; touch "+output.output2
          shell(cmd)
        else:
          print("Dupa Jasio")


rule trimgalore_se:
    input:
        infile = get_trimgalore_input
    output:
        output = DIR_trimmed+"{sample}_trimmed.fq.gz" #here it doesnt matter what is here, I need it only for wildcard
    params:
        extra          = config.get("trim_galore_args", ""),
        outdir         = "--output_dir " + DIR_trimmed,
        gz             = "--gzip",
        cutadapt       = "--path_to_cutadapt "+ CUTADAPT,
    log:
        log=DIR_trimmed+"{sample}.trimgalore.log"
    message:
        " ---------  Trimming raw single-end read data using {TRIMGALORE} -------  "
    run:
        args = " ".join(
        [params.extra,
         params.outdir,
         params.gz,
         params.cutadapt,
         "-o "+ DIR_trimmed
         ])
        cmd=TRIMGALORE+args + " "+input.infile[0]+" 2> "+ log.log
        cmd = cmd + "; touch "+output.output
        shell(cmd)


def get_fastq_input(wc):
   filename = PATHIN+wc.name+".fq.gz"
   return(filename)
    

rule fastqc_raw:
    input:
        infile = get_fastq_input
    output:
        DIR_rawqc + "{name}_fastqc.html",
        DIR_rawqc + "{name}_fastqc.zip"
    params:
        fastqc_args = config.get("fastqc_args", ""),
        outdir = "--outdir "+ DIR_rawqc    # usually pass params as strings instead of wildcards.

    log:
        log=DIR_rawqc + "{name}_fastqc.log"
    message: """ ----------  Quality checking raw read data with {FASTQC}.  --------------   """
    run:
        cmd=FASTQC + " "+ params.fastqc_args + " "+ params.outdir + " "+  input.infile + " "+ " 2> "+ " "+log.log
        shell(cmd)


#print("------------cmd----------")
#print(cmd)
#"touch sampleid for fastqc or rename"
#"for bowtie vdran renames outfile from unit.bam to sample.bam"
