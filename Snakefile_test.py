#!/usr/bin/env python3.5

#============================================================================================================
# SNAKEMAKE FILE WRITTEN BY THE AKALIN GROUP AT MDC, BERLIN, 2017
# Alexander Gosdschan, Katarzyna Wreczycka, Bren Osberg
# To process bisulfite sequencing data from raw fastq files to performing integrated bioinformatics analysis.
#============================================================================================================

# snakemake  --forceall --snakefile /home/kwreczy/repositories/pigx_bsseq//Snakefile_test.py --cores 2 --config tablesheet=/home/kwreczy/repositories/pigx_bsseq/test_dataset/TableSheet_test.csv gtoolbox=/home/kwreczy/tmp/ggtoolbox/ in=/home/kwreczy/repositories/pigx_bsseq/test_dataset/in/ out=~/tmp/my_output/ genome_folder=~/tmp/ce10/ log=/home/kwreczy/logs/ chrominfo=/home/kwreczy/repositories/pigx_bsseq/test_dataset/chromInfo.txt bismark_args="  --bowtie2 -N 1 -L 2 --genome_folder ~/tmp/ce10/"
# snakemake  --forceall --snakefile /home/kwreczy/repositories/pigx_bsseq//Snakefile_test.py --config configfile=my_output~/tmp/my_output/config.json
# snakemake  --forceall --snakefile /home/kwreczy/repositories/pigx_bsseq//Snakefile_test.py --cores 2 --config tablesheet=/home/kwreczy/repositories/pigx_bsseq/test_dataset/TableSheet_SRA.csv gtoolbox=/home/kwreczy/tmp/ggtoolbox/ in=/home/kwreczy/repositories/pigx_bsseq/test_dataset/in/ out=~/tmp/my_output/ genome_folder=~/tmp/ce10/ log=/home/kwreczy/logs/ chrominfo=/home/kwreczy/repositories/pigx_bsseq/test_dataset/chromInfo.txt bismark_args="  --bowtie2 -N 1 -L 2 --genome_folder ~/tmp/ce10/"


#------ set config file, include function definitions, and set os:
import os,csv,json
include   : "./scripts/func_defs.py"
include   : "./scripts/functions_parsing.py"
include   : "./rules/SRA2fastq/SRA2fastq_functions.py"


# Check configuration parameters  
config = parseConfig(config)

# Read a Table Sheet in the CSV format
TABLESHEET =  config["tablesheet"]
if(TABLESHEET is not None):
  
  # Parse Sheet Table provided by the user
  (rows, sample_ids, list_units) = parseTable(TABLESHEET)
  firstcol = [r[0] for r in rows]
  
  # Add samples, units and treatment information from a sheet table to config
  config["samples"] = dict(zip(sample_ids, sample_ids))
  config["units"] = dict(zip(sample_ids, list_units))
  config["treatment"] = dict(zip(sample_ids, [r[4] for r in rows][1:]))
  
  # All tools have to be accessible from the given 'progs' directory
  config["progs"] = {"FASTQC"          : "fastqc",
                 "TRIMGALORE"      : "trim_galore",
                 "CUTADAPT"        : "cutadapt",
                 "BOWTIE2"         : "bowtie2" ,
                 "BISMARK"         : "bismark",
                 "DEDUPLICATE_BISMARK" : "deduplicate_bismark",
                 "BISMARK_GENOME_PREPARATION"     : "bismark_genome_preparation",
                 "BISMARK_METHYLATION_EXTRACTOR"  : "bismark_methylation_extractor",
                 "BISMARK2REPORT"                 : "bismark2report",
                 "SAMTOOLS"                       : "samtools"}
  
  # Save config file
  # Config file will be updated if there is at least one SRA id.
  with open(config['out']+"config.json", 'w') as outfile:
    json.dump(config, outfile)



######### BEGIN the SRA part: Download fastq files based on their SRA ids


# Check if in first column of a csv table there are SRA ids
acc = ['PRJ', #Study accession
       #'SAMN', #Sample accession
       'SRS', #Secondary sample accession 
       'SRX', #Experiment accession
       'SRR'] #Run accession
SRA2download = [] # rows with SRA ids
SRA2download_indx = [] # indecies of rows with SRA ids
firscol = [x[0] for x in config["units"].values()]
for i in range(len(firstcol)):
  if ( (firstcol[i][:3] in acc ) or ( firstcol[i][:4]=='SAMN') ) and ( check_if_fastq(firstcol[i]) is True ):
    SRA2download.append(firstcol[i])
    SRA2download_indx.append(i)

    
if len(SRA2download)!=0:
  
  # Fastq files suppose to be in 'in' directory.
  # If SRA ids were provided, fastq file will be found
  # in the ENA database and SRA ids will be replaced in the
  # table as fastq files names.

 # Get indecies of rows with SRA ids
  sra2down_header_indx=[0]
  sra2down_header_indx.extend(SRA2download_indx)
  sra2down_header = [rows[i] for i in sra2down_header_indx]  
  # Get indecies of rows with no SRA ids
  nosra2down = [rows[i] if i not in sra2down_header_indx else None for i in range(len(rows))]
  nosra2down = list(filter(None.__ne__, nosra2down)) #remove None elements from 'nosra2down' list

  dict_SRA_ftp = get_dict_SRA_ftp(SRA2download, config['in'])
  sra_new_rows = create_new_rows(sra2down_header, dict_SRA_ftp)
  new_sheet_rows = [rows[0]] # header
  new_sheet_rows.extend(nosra2down) # non-SRA rows
  new_sheet_rows.extend(sra_new_rows) # SRA rows

  # Save updated Table Sheet to th output directory
  TABLESHEET_NEW = config['out'] +'TableSheet_config.csv'
  with open(TABLESHEET_NEW, 'w', newline='') as f: # TODO: remove it?
    writer = csv.writer(f)
    writer.writerows(new_sheet_rows)

  # Update config dictionary
  sample_ids_new = [x[2] for x in new_sheet_rows[1:]]
  list_units_new = [  list(filter(None,[x[0],x[1]])) for x in new_sheet_rows[1:]  ]
  config["samples"] = dict(zip(sample_ids_new, sample_ids_new))
  config["units"] = dict(zip(sample_ids_new, list_units_new))
  config["treatment"] = dict(zip(sample_ids_new, [r[4] for r in new_sheet_rows][1:]))

  ftplinks = sum(list(dict_SRA_ftp.values()), [])
  sraids = [filter_filename_no_ext( x ) for x in ftplinks]
  config["ftp_sra"] = dict(zip(sraids, ftplinks))
  
  # Update config file and save it in the output directory
  with open(config['out']+"config.json", 'w') as outfile: # TODO: remove it?
     json.dump(config, outfile)
  
  # Download fastq files here based on their SRA ids
  include: "rules/SRA2fastq/Snakefile"


######### END the SRA part


# Extract and remove extensions from list of units
samplesids = config['units'].keys() # dicts in python are not ordered
#There has to be nicer way than this
#Because right now, every time when there is at least 1 SRA ids in the sheet
# only .fastq.gz suffixes can be used
(list_units, inext) = remove_ext_from_units([config['units'][x] for x in samplesids])
new_units = dict(zip(samplesids, list_units))
config["inext"] = inext # Todo and what if file is not gzipped?
for x in samplesids:
  config['units'][x] = new_units[x]
  

#---------------------------------     DEFINE PATHS AND FILE NAMES:  ----------------------------------

PATHIN          = config["in"]         #--- location of the data files to be imported
PATHOUT         = config["out"]        #--- where to send the output
GTOOLBOX        = config["gtoolbox"]       #--- where the programs are stored to carry out the necessary operations
LOGS            = config["log"]           #--- subfolder name for the logs of some programs

VERSION         = config["genome_version"]  #--- version of the genome being mapped to.

CHROM_INFO      = config["chrominfo"]     #--- details of the reference genome (length, etc.) haploid chroms have been removed.
NUMTHREADS      = config["numthreads"]


INEXT=config["inext"] #TODO: it can be more flexible
RCODE = [".pe1", ".pe2"] #TODO: would be awesome if it could be more flexible, bbecause right now its impossible to mix it with SRA
#RCODE = ["_1", "_2"] #TODO: this si for SRA IDS
#RCODE_pref = {'sampleB':["pe1.", "pe2."], 'sampleB1':["pe1."]} 
#RCODE_suff = {'sampleB':["", "pe2."], 'sampleB1':[""]}


#-------------------------------      DEFINE PROGRAMS TO BE EXECUTED: ---------------------------------

FASTQC                         =  GTOOLBOX+config["progs"]["FASTQC"]            #--- self-explanatory program names.
TRIMGALORE                     =  GTOOLBOX+config["progs"]["TRIMGALORE"]
CUTADAPT                       =  GTOOLBOX+config["progs"]["CUTADAPT"]
BISMARK_GENOME_PREPARATION     =  GTOOLBOX+config["progs"]["BISMARK_GENOME_PREPARATION"]
BISMARK                        =  GTOOLBOX+config["progs"]["BISMARK"]
BOWTIE2                        =  GTOOLBOX+config["progs"]["BOWTIE2"]
DEDUPLICATE_BISMARK            =  GTOOLBOX+config["progs"]["DEDUPLICATE_BISMARK"]
BISMARK_METHYLATION_EXTRACTOR  =  GTOOLBOX+config["progs"]["BISMARK_METHYLATION_EXTRACTOR"]
BISMARK2REPORT                 =  GTOOLBOX+config["progs"]["BISMARK2REPORT"]

SAMTOOLS                       =  GTOOLBOX+config["progs"]["SAMTOOLS"] 


#---------------------------     LIST THE OUTPUT DIRECTORIED AND SUBDIRECTORIED TO BE PRODUCED     ------------------------------

for x in [config["in"], config["out"], config["log"]]:
  if not os.path.exists(x):
    os.makedirs(x)  
  
DIR_sorted=PATHOUT+'06_sorted/'
DIR_mapped_n_deduped=PATHOUT+'04_mapped_n_deduped/'
DIR_posttrim_QC=PATHOUT+'03_posttrim_QC/'
DIR_trimmed=PATHOUT+'02_trimmed/'
DIR_rawqc=PATHOUT+'01_rawqc/'

for x in [DIR_rawqc, DIR_trimmed, DIR_posttrim_QC, DIR_mapped_n_deduped, DIR_sorted]:
  if not os.path.exists(x):
    os.makedirs(x)


#---------------------------     LIST THE OUTPUT FILES TO BE PRODUCED     ------------------------------


# TODO: acutally here should be instead of config["samples"].keys() then config["samples"].values()....


OUTPUT_FILES = [
      #               ======  rule 01 raw QC    =========
      #[ expand (list_files(DIR_rawqc, config["units"][x], "_fastqc.html")  ) for x in config["samples"].keys()  ],
      #----RULE 2 IS ALWAYS EXECUTED, TRIMMING IS A PREREQUISITE FOR SUBSEQUENT RULES ----
      #               ======  rule 03 posttrim_QC_ ======
      #[ expand ( list_files_posttrim_QC(DIR_posttrim_QC, config["units"][x],".html")  ) for x in config["samples"].keys()  ],
      #--- fastQC output files are not needed downstream and need to be called explicitly.
      #               ====rule 02 trimgalore ======
      #[ expand ( list_files_TG( DIR_trimmed, config["units"][x])  ) for x in config["samples"].keys()  ],
      #               ====rule 04 Mapping ======
      [ expand ( list_files_bismark( DIR_mapped_n_deduped, config["units"][x] )  ) for x in config["samples"].keys()  ],
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
    

# ==========================================================================================

# #TODO: remove it. I created it only to force sortindex rule be ran be snakemake.
# rule dummy_bai:
#     input:
#         DIR_sorted+"{sample}_trimmed_bismark_bt2.deduplicated.sorted.bam.bai"
# 
# 
# rule sortindex:
#     input:
#         DIR_sorted+"{sample}_trimmed_bismark_bt2.deduplicated.sorted.bam"
#     output:
#         "{DIR_sorted}/{sample}_trimmed_bismark_bt2.deduplicated.sorted.bam.bai"
#     threads: 2
#     shell: 
#         "{SAMTOOLS} index {input}"    
# 
# 
# # ==========================================================================================
# 
# rule sortbam:
#     input:
#         DIR_mapped_n_deduped+"{sample}_trimmed_bismark_bt2.deduplicated.bam"
#     output:
#         "{DIR_sorted}/{sample}_trimmed_bismark_bt2.deduplicated.sorted.bam"
#     threads: 2
#     shell: 
#         "{SAMTOOLS} sort --threads {threads} {input} -o {output}"    


# ==========================================================================================

#--------
# rule bismark_pe_deduplication:
#     input:
#         DIR_mapped_n_deduped+"{sample}"+RCODE[0] +reads[0]+"_bismark_bt2_pe.bam",
#         DIR_mapped_n_deduped+"{sample}"+RCODE[0] +reads[0]+"_bismark_bt2_PE_report.txt",
#         DIR_mapped_n_deduped+"{sample}"+RCODE[1] +reads[1]+"_bismark_bt2_pe.bam",
#         DIR_mapped_n_deduped+"{sample}"+RCODE[1] +reads[1]+"_bismark_bt2_PE_report.txt"
#     output:
#         DIR_mapped_n_deduped+"{sample}"+RCODE[0] +reads[0]+"_bismark_bt2_pe.deduplicated.bam",
#         DIR_mapped_n_deduped+"{sample}"+RCODE[1] +reads[1]+"_bismark_bt2_pe.deduplicated.bam", ##TODO: add report
#     log:
#         DIR_mapped_n_deduped+"{sample}_deduplication.log"
#     message: """-----------   Deduplicating paired-end read alignments ---------------------- """
#     shell:
#         "{samtools} rmdup {input} {output} 2> {log}"
#         
#         
# rule bismark_se_deduplication:
#     input:
#         DIR_mapped_n_deduped+"{sample}_trimmed_bismark_bt2.bam"
#     output:
#         "{DIR_mapped_n_deduped}/{sample}_trimmed_bismark_bt2.deduplicated.bam",
#         "{DIR_mapped_n_deduped}/{sample}_trimmed_bismark_bt2.deduplication_report.txt"
#     params:
#         bam="--bam ",
#         sampath="--samtools_path "+SAMTOOLS
#     log:
#         "{DIR_mapped_n_deduped}/{sample}_deduplication.log"
#     message: """-----------   Deduplicating single-end read alignments ---------------------- """
#     shell:
#         "{DEDUPLICATE_BISMARK} {params} {input} 2> {log} "


# ==========================================================================================

rule bismark_pe:
    input:
        fin1=DIR_trimmed+"{sample}"+RCODE[0]+"_val_1.fq.gz",
        fin2=DIR_trimmed+"{sample}"+RCODE[1]+"_val_2.fq.gz"
    output:
        DIR_mapped_n_deduped+"{sample}"+RCODE[0]+"_val_1_bismark_bt2_pe.bam",
        DIR_mapped_n_deduped+"{sample}"+RCODE[0]+"_val_1_bismark_bt2_PE_report.txt"
    threads: 2
    params:
        extra        = config.get("bismark_args", ""),
        useBowtie2   = "--bowtie2 ",
        outdir       = "--output_dir  "+ DIR_mapped_n_deduped,
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(BOWTIE2),
        samtools     = "--samtools_path "+ os.path.dirname(SAMTOOLS),
        tempdir      = "--temp_dir " + PATHOUT,
    log:
        DIR_mapped_n_deduped+"{sample}_bismark_pe_mapping.log"
    message: """-------------   Mapping paired-end reads to genome {VERSION}. ------------- """
    shell:
        "{BISMARK} {params} --multicore {threads} -1 {input.fin1} -2 {input.fin2} 2> {log}"


rule bismark_se:
    input:
        DIR_trimmed+"{sample}_trimmed.fq.gz"
    output:
        "{DIR_mapped_n_deduped}/{sample}_trimmed_bismark_bt2.bam",
        "{DIR_mapped_n_deduped}/{sample}_trimmed_bismark_bt2_SE_report.txt"
    threads: 2
    params:
        extra       = config.get("bismark_args", ""),
        outdir      = "--output_dir " + DIR_mapped_n_deduped,
        tempdir     = "--temp_dir " + PATHOUT,
        sampath     ="--samtools_path " + os.path.dirname(SAMTOOLS),
        useBowtie2  = "--bowtie2 ",
        bowtie2path = "--path_to_bowtie " + os.path.dirname(BOWTIE2)
    log:
        "{DIR_mapped_n_deduped}/{sample}_bismark_se_mapping.log"
    message: """-------------   Mapping single-end reads to genome {VERSION}. ------------- """
    shell:
        "{BISMARK} {params} --multicore {threads} {input} 2> {log}"



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
    message: """ ------------  Quality checking trimmmed single-end data with Fastqc ------------- """
    shell:
       "{FASTQC} {params.outdir} {input} 2> {log}"
#--------
rule fastqc_after_trimming_pe:
    input:
        DIR_trimmed+"{sample}"+RCODE[0]+"_val_1.fq.gz",
        DIR_trimmed+"{sample}"+RCODE[1]+"_val_2.fq.gz"
    output:
    	DIR_posttrim_QC+"{sample}"+RCODE[0]+"_val_1_fastqc.html",
    	DIR_posttrim_QC+"{sample}"+RCODE[0]+"_val_1_fastqc.zip",
    	DIR_posttrim_QC+"{sample}"+RCODE[1]+"_val_2_fastqc.zip",
        DIR_posttrim_QC+"{sample}"+RCODE[1]+"_val_2_fastqc.html"
    params:
        fastqc_args = config.get("fastqc_args", ""),
        outdir = "--outdir "+DIR_posttrim_QC
    log:
   	    DIR_posttrim_QC+"{sample}_trimmed_fastqc.log"
    message: """ ------------  Quality checking trimmmed paired-end data with Fastqc ------------- """
    shell:
        "{FASTQC} {params.outdir} {input} 2> {log}"
        
# # #
# # # ==========================================================================================
# # # trim the reads
# # 

rule trimgalore_se:
   input:
       PATHIN+"{sample}"+INEXT
   output:
       DIR_trimmed+"{sample}_trimmed.fq.gz" #---- this ALWAYS outputs .fq.qz format.
   params:
       extra        = config.get("trim_galore_args", ""),
       outdir       = "--output_dir "+  DIR_trimmed,
       phred        = "--phred33",
       gz           = "--gzip",
       cutadapt     = "--path_to_cutadapt "+ CUTADAPT,
   log:
       DIR_trimmed+"{sample}.trimgalore.log"
   message:
       " ---------  Trimming raw single-end read data using {TRIMGALORE} -------  "
   shell:
       "{TRIMGALORE} {params} {input} -o {DIR_trimmed} 2> {log}"


rule trimgalore_pe:
    input:
        PATHIN + "{sample}" + RCODE[0] + INEXT,
        PATHIN + "{sample}" + RCODE[1] + INEXT
    output:
        DIR_trimmed+"{sample}" + RCODE[0] + "_val_1.fq.gz", #---- this ALWAYS outputs .fq.qz format.
        DIR_trimmed+"{sample}" + RCODE[1] + "_val_2.fq.gz",
    params:
        extra          = config.get("trim_galore_args", ""),
        outdir         = "--output_dir " + DIR_trimmed,
        gz             = "--gzip",
        cutadapt       = "--path_to_cutadapt "+ CUTADAPT,
        paired         = "--paired"
    log:
        DIR_trimmed+"{sample}.trimgalore.log"
    message:
        " ---------  Trimming raw paired-end read data using {TRIMGALORE} -------  "
    shell:
        "echo {output}; {TRIMGALORE} {params} {input} -o {DIR_trimmed} 2> {log}"


# ==========================================================================================
# raw quality control 
#

rule fastqc_raw: #----only need one: covers BOTH PE and SE cases.
    input:
        PATHIN+"{sample}"+INEXT
    output:
        DIR_rawqc + "{sample}_fastqc.html",
        DIR_rawqc + "{sample}_fastqc.zip"
    params:
        fastqc_args = config.get("fastqc_args", ""),
        outdir = "--outdir "+ DIR_rawqc    # usually pass params as strings instead of wildcards.

    log:
        DIR_rawqc + "{sample}_fastqc.log"
    message: """ ----------  Quality checking raw read data with {FASTQC}.  --------------   """
    shell:
        "{FASTQC} {params.outdir}  {input} 2> {log}"

        
