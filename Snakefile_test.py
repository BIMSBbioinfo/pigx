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
#include   : "bismark_rule.py"


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
      [ expand (list_files(DIR_rawqc, getFilenames(config["SAMPLES"][x]['fastq']), x,"_fastqc.html")  ) for x in config["SAMPLES"].keys()  ],
      #----RULE 2 IS ALWAYS EXECUTED, TRIMMING IS A PREREQUISITE FOR SUBSEQUENT RULES ----
      #               ======  rule 03 posttrim_QC_ ======
      #[ expand ( list_files_posttrim_QC(DIR_posttrim_QC, config["SAMPLES"][x]['fastq'],x, ".html")  ) for x in config["SAMPLES"].keys()  ],
      #--- fastQC output files are not needed downstream and need to be called explicitly.
      #               ====rule 02 trimgalore ======
      [ expand (list_files_TG(DIR_trimmed, getFilenames(config["SAMPLES"][x]['fastq']), x  )) for x in config["SAMPLES"].keys()  ],
      #               ====rule 04 Mapping ======
      [ expand ( list_files_bismark( DIR_mapped_n_deduped, getFilenames(config["SAMPLES"][x]['fastq']), x  )) for x in config["SAMPLES"].keys()  ]
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

def get_fq_after_trimming(wc):

   print('---------------wc------get_fq_after_trimming-----2222')
   print(list(wc))

   samps = config['SAMPLES'][wc.sample]['fastq_name']

   if type(samps) is str:
        samps = [samps]

   if(len(samps)==2):
     d=[os.path.join(DIR_trimmed, samps[0]+'_val_1.fq.gz'), os.path.join(DIR_trimmed, samps[1]+'_val_2.fq.gz')]
   else:
     d=[os.path.join(DIR_trimmed, samps[0]+'_trimmed.fq.gz')]
   print('d')
   print(d)
   return(d)

rule bismark_se:
    input:
        infile = get_fq_after_trimming
        #infile = lambda wc:  [DIR_trimmed + config['SAMPLES'][wc.sample]['fastq_name'][0]+'_trimmed.fq.gz']
    output:
        o=DIR_mapped_n_deduped+"{sample}.bam"
    threads: 2
    params:
        program     = BISMARK,
        extra       = config.get("bismark_args", ""),
        outdir      = "--output_dir " + DIR_mapped_n_deduped,
        tempdir     = "--temp_dir " + PATHOUT,
        sampath     ="--samtools_path " + os.path.dirname(SAMTOOLS),
        bowtie2path = "--path_to_bowtie " + os.path.dirname(BOWTIE2),
        useBowtie2  = "--bowtie2 ",
        genome_folder=GENOMEPATH
    log:
        log=DIR_mapped_n_deduped+"{sample}_bismark_mapping.log"
    message: """-------------   Mapping single-end reads to genome {VERSION}. ------------- """
    run:
        print('--------bismark-------infile---------')
        print(input.infile)

        cmd1 = " ".join(
        [BISMARK,
         params.extra,
         params.outdir,
         params.tempdir,
         params.sampath,
         params.bowtie2path,
         params.useBowtie2,
         params.genome_folder,
         '--multicore', str(threads)
         ])

        if len(input.infile)==1:
          cmd2 = " ".join(
            [" ", input.infile[0],
             '2>',log.log
             ])
        elif len(input.infile)==2:
          cmd2 = " ".join(
            [" ", " -1 ", input.infile[0],
             " -2 ", input.infile[1],
             '2>',log.log
             ])

        #real_list_files_bismark(input.infile)
        command = cmd1 + cmd2 + "; touch " +  output.o
        print("----------commad bismark--------")
        print(command)
        shell(command)

#        unit.pe1_trimmed.bismark_bt2.bam
#        mv unit.pe1_trimmed.bismark_bt2.bam {sample}.bismark_bt2.bam
        
# # 
# # 
# # 
# # # ==========================================================================================
# # # post-trimming quality control
# # 

######### this si wrong, change it!!!!


# rule fastqc_after_trimming_se:
#     input:
#        get_fastqc_after_trimming_input
#     output:
#     	  DIR_posttrim_QC+"{sample}_trimmed_fastqc.html",
#     	  DIR_posttrim_QC+"{sample}_trimmed_fastqc.zip"
#     params:
#        fastqc_args = config.get("fastqc_args", ""),
#        outdir = "--outdir "+DIR_posttrim_QC
#     log:
#    	   DIR_posttrim_QC+"{sample}_trimmed_fastqc.log"
#     message: """ ------------  Quality checking trimmmed single-end data with Fastqc ------------- """
#     run:
#         cmd=FASTQC + " "+ params.fastqc_args + " "+ params.outdir + " "+  input.infile[0] + " "+ " 2> "+ " "+log.log
#         # Here a hack to generate files named like samples and not like units, even though
#         # trimgalore will produce files named using units
#         cmd = cmd + "; for x in {output.o1} {output.o2}; do touch $x; done"
#         # Execute
#         shell(cmd)
       
#--------
# rule fastqc_after_trimming_pe:
#     input:
#         infile=get_fastqc_after_trimming_input
#     output:
#     	o1=DIR_posttrim_QC+"{sample}"+"_val_1_fastqc.html",
#     	o2=DIR_posttrim_QC+"{sample}"+"_val_1_fastqc.zip",
#     	o3=DIR_posttrim_QC+"{sample}"+"_val_2_fastqc.zip",
#         o4=DIR_posttrim_QC+"{sample}"+"_val_2_fastqc.html"
#     params:
#         fastqc_args = config.get("fastqc_args", ""),
#         outdir = "--outdir "+DIR_posttrim_QC
#     log:
#    	    DIR_posttrim_QC+"{sample}_trimmed_fastqc.log"
#     message: """ ------------  Quality checking trimmmed paired-end data with Fastqc ------------- """
#     #shell:
#     #    "{FASTQC} {params.outdir} {input} 2> {log}"
#     run:
#         cmd=FASTQC + " "+ params.fastqc_args + " "+ params.outdir + " "+  input.infile[0] + " "+ " 2> "+ " "+log.log
#         cmd=cmd+" ; "+ FASTQC + " "+ params.fastqc_args + " "+ params.outdir + " "+  input.infile[1] + " "+ " 2> "+ " "+log.log
#         # Here a hack to generate files named like samples and not like units, even though
#         # trimgalore will produce files named using units
#         cmd = cmd + "; for x in {output.o1} {output.o2} {output.o3} {output.o4}; do touch $x; done"
#         # Execute
#         shell(cmd)

#
# #
# # ==========================================================================================
# # trim the reads
# 
# def get_trimgalore_input(wc):
# 
#    print('---------------wc------get_trimgalore_input-----1111')
#    print(list(wc))
#    samps = config['SAMPLES'][wc.sample]['fastq']
# 
#    if type(samps) is str:
#         samps = [samps]
# 
#    if(len(samps)==2):
#      return( [os.path.join(PATHIN, samps[0]), os.path.join(PATHIN, samps[1])] )
#    else:
#      return( [os.path.join(PATHIN, samps[0])] )


rule trimgalore_pe:
    input:
        #infile = get_trimgalore_input
        #infile = lambda wc: [os.path.join(PATHIN, config['SAMPLES'][wc.sample]['fastq_name'][0]), os.path.join(PATHIN, config['SAMPLES'][wc.sample]['fastq_name'][1])] if len(config['SAMPLES'][wc.sample]['fastq_name'])==2 else [os.path.join(PATHIN, config['SAMPLES'][wc.sample]['fastq_name'][0])]
        infile = lambda wc: [os.path.join(PATHIN, config['SAMPLES'][wc.sample]['fastq_name'][0]), os.path.join(PATHIN, config['SAMPLES'][wc.sample]['fastq_name'][1])]
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

        print("---------------------print trim galore rule")
        if(len(input.infile)==0):
          print("AAAAAAAAAAA")
          print(output)
          print(intput)

        cmd = " ".join(
        [TRIMGALORE,
         " --paired",
         input.infile[0], input.infile[1],
         params.extra, params.outdir, params.gz, params.cutadapt,
         "-o "+ DIR_trimmed,
         " 2> ", log.log
         ])

        # Here a hack to generate files named like samples and not like units, even though
        # trimgalore will produce files named using units
        cmd = cmd + "; touch "+output.output1 + "; touch "+output.output2
        shell(cmd)




rule trimgalore_se:
    input:
        #infile = get_trimgalore_input
        #infile = lambda wc: os.path.join(PATHIN, config['SAMPLES'][wc.sample]['fastq_name'][0])
        #infile = lambda wc: [wc.sample]
        infile=PATHIN+'{sample}'+'.fq.gz'
        #infile = lambda wc: [PATHIN +config['SAMPLES'][wc.sample]['fastq_name'][0], os.path.join(PATHIN, config['SAMPLES'][wc.sample]['fastq_name'][1])] if len(config['SAMPLES'][wc.sample]['fastq_name'])==2 else [PATHIN+config['SAMPLES'][wc.sample]['fastq_name'][0]+".fq.gz"]
    output:
        output = DIR_trimmed+"{sample}_trimmed.fq.gz"
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
        cmd=TRIMGALORE+args + " "+input.infile+" 2> "+ log.log
        cmd = cmd + "; touch "+output.output
        shell(cmd)


# def get_fastq_input(wc):
#    print('---------------wc------fastqc-----000000000000000')
#    # wc.name is e.g. single.pe1
#    filename = PATHIN+wc.name+".fq.gz" #TODO: not all files should have to have .fq.gz suffix
#    return(filename)


rule fastqc_raw:
    input:
        #infile = get_fastq_input
        #infile = lambda wc: PATHIN+wc.name+".fq.gz" 
        infile=PATHIN+"{name}"+".fq.gz"
    output:
        DIR_rawqc + "{name}_fastqc.html",
        DIR_rawqc + "{name}_fastqc.zip"
    params:
        fastqc_args = config.get("fastqc_args", ""),
        outdir = "--outdir "+ DIR_rawqc

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
