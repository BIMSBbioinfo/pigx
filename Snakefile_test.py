# -*- coding:latin-1 -*-


## Note: I ran it like this:
#/home/kwreczy/.local/bin/snakemake --snakefile ./Snakefile_test.py --cores 8 --forceall --config tablesheet=/home/kwreczy/repositories/makeNGSnake/test/test_dataset/TableSheet_SRA.csv gtoolbox=./programs/ in=/home/kwreczy/repositories/makeNGSnake/test/test_dataset/ out=/data/akalin/kwreczy/my_output/ genome_folder=/data/akalin/Base/Genomes/ce10/ log=/data/akalin/kwreczy/logs/ chrominfo=/home/kwreczy/repositories/makeNGSnake/test/test_dataset/chromInfo.txt bismark_args=" --non-directional --PBAT"
##


import os,csv
from scripts.functions import *
from rules.SRA2fastq.SRA2fastq_functions import *

# Check configuration parameters  
config = parseConfig(config)

if not os.path.exists(config["in"]):
    os.makedirs(config["in"])
if not os.path.exists(config["out"]):
    os.makedirs(config["out"])
if not os.path.exists(config["log"]):
    os.makedirs(config["log"])
    
TABLESHEET =   config["tablesheet"]

# Parse Sheet Table provided by the user
(rows, sample_ids, list_units) = parseTable(TABLESHEET)
firstcol = [r[0] for r in rows]

# Add samples, units and treatment information from a sheet table to config
config["SAMPLES"] = dict(zip(sample_ids, sample_ids))
config["UNITS"] = dict(zip(sample_ids, list_units))
config["TREATMENT"] = dict(zip(sample_ids, [r[4] for r in rows][1:]))

# Save config file
# Config file will be updated if there is at least one SRA id.
config=config2JSON(config)
with open(config['PATHS']['PATHOUT']+"config.json", 'w') as outfile:
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
for i in range(len(firstcol)):
  if ( (firstcol[i][:3] in acc ) or ( firstcol[i][:4]=='SAMN') ) and ( check_if_fastq(firstcol[i]) is True ):
    SRA2download.append(firstcol[i])
    SRA2download_indx.append(i)

if len(SRA2download)!=0:
  
  # Fastq files suppose to be in PATHIN directory.
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

  dict_SRA_ftp = get_dict_SRA_ftp(SRA2download, config['PATHS']['PATHIN'])
  sra_new_rows = create_new_rows(sra2down_header, dict_SRA_ftp)
  new_sheet_rows = [rows[0]] # header
  new_sheet_rows.extend(nosra2down) # non-SRA rows
  new_sheet_rows.extend(sra_new_rows) # SRA rows

  # Save updated Table Sheet
  TABLESHEET_NEW = config['PATHS']['PATHOUT'] +'TableSheet_config.csv' # TODO
  with open(TABLESHEET_NEW, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(new_sheet_rows)

  # Update config dictionary
  sample_ids_new = [x[2] for x in new_sheet_rows[1:]]
  list_units_new = [  list(filter(None,[x[0],x[1]])) for x in new_sheet_rows[1:]  ]
  config["SAMPLES"] = dict(zip(sample_ids_new, sample_ids_new))
  config["UNITS"] = dict(zip(sample_ids_new, list_units_new))
  config["TREATMENT"] = dict(zip(sample_ids_new, [r[4] for r in new_sheet_rows][1:]))
  ftplinks = sum(list(dict_SRA_ftp.values()), [])
  sraids = [filter_filename_no_ext( x ) for x in ftplinks]
  config["FTP_SRA"] = dict(zip(sraids, ftplinks))
  # Update config file
  with open(config['PATHS']['PATHOUT']+"config.json", 'w') as outfile:
     json.dump(config, outfile)
  
  # Download fastq files here based on their SRA ids
  include: "rules/SRA2fastq/Snakefile"
  
else:

  rule create_dummy_log:
    output:
        "dummy.txt"
    run:
        os.system("touch dummy.txt")

######### END the SRA part



#---------------------------------     DEFINE PATHS AND FILE NAMES:  ----------------------------------

RCODE           = "TMP"                   #--- the string that denotes which read the file corresponds to (only relevent for paired-end experiments.)
PATHIN          = config["PATHS"]["PATHIN"]         #--- location of the data files to be imported
PATHOUT         = config["PATHS"]["PATHOUT"]        #--- where to send the output
GTOOLBOX        = config["PATHS"]["GTOOLBOX"]       #--- where the programs are stored to carry out the necessary operations
GENOMEPATH      = config["PATHS"]["GENOMEPATH"]     #--- where the reference genome being mapped to is stored
LOGS            = config["PATHS"]["LOGS"]           #--- subfolder name for the logs of some programs

INEXT           = "TMP"                   #--- input file extension; usually .fq.gz, but can also be .bz2 among other possibilities.
VERSION         = config["GENOMEDAT"]["VERSION"]        #--- version of the genome being mapped to.

CHROM_INFO      = config["GENOMEDAT"]["CHROM_INFO"]     #--- details of the reference genome (length, etc.) haploid chroms have been removed.
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


#---------------------------     LIST THE OUTPUT FILES TO BE PRODUCED     ------------------------------

rule final: 
  input: 
    "anotherdummy.txt"
        
rule stupid_rule:
    input:
        "dummy.txt"
    output:
        "anotherdummy.txt"
    params:
      extra=config["bismark_args"]
    shell: "echo {params.extra} > anotherdummy.txt "













        
