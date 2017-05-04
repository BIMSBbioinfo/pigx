

def parseConfig(config_args):
  #      
  if ( config_args=={} ):
    raise Exception( """Missing arguments indicating 'tablesheet', 'gtoolbox','in','out','genome_folder', 'log' and 'chrominfo'.
                      Use argument --config and then type paths to the arguments, e.g.:
                      snakemake -s Snakemake.py --config tablesheet=./TableSheet.csv gtoolbox=./programs/ in=./my_fastq_files/ out=./my_output/ genome_folder=./Genomes/hg19/ log=./logs/ chrominfo=./hg19.chrom.sizes.ucsc.txt numthreads=2""" )
  #      
  tablesheet = config_args.get("tablesheet")
  if (tablesheet is None):
     raise Exception( """Missing an argument indicating 'tablesheet'.
                      Use argument --config and then type paths to the arguments, e.g.:
                      snakemake -s Snakemake.py --config tablesheet=./TableSheet.csv gtoolbox=./programs/ in=./my_fastq_files/ out=./my_output/ genome_folder=./Genomes/hg19/ log=./logs/  chrominfo=./hg19.chrom.sizes.ucsc.txt numthreads=2""" )
  #
  gtoolbox = config_args.get("gtoolbox")
  pathin = config_args.get("in")
  pathout = config_args.get("out")
  logs = config_args.get("log")
  config_args["gtoolbox"] = gtoolbox if gtoolbox is not None else "/bin/"
  config_args["in"] = pathin if pathin is not None else "./"
  config_args["out"] = pathout if pathout is not None else "./"
  config_args["log"] = logs if logs is not None else "./"
  #
  version = config_args.get("genome_version")
  config_args["version"] = version if version is not None else ""
  numthreads = config_args.get("numthreads")
  config_args["numthreads"] = numthreads if numthreads is not None else 2
  directional = config_args.get("directional")
  config_args["directional"] = directional if directional is not None else True
  #
  return(config_args)


def config2JSON(conf):
  import json
  # filter(None, my_list) removes empty elements from the list my_list
  if conf['SAMPLES'] is None:
     raise Exception( """Config file has to have a key called 'SAMPLES'.""" )
  if conf['TREATMENT'] is None:
     data['TREATMENT']=""
  #  
  conf['PATHS'] = {"GTOOLBOX":conf["gtoolbox"],
                 "PATHIN":conf["in"],
                 "PATHOUT":conf["out"],
                 #"GENOMEPATH":conf["genome_folder"],
                 "LOGS":conf["log"]}
  conf['GENOMEDAT'] = {"CHROM_INFO":conf["chrominfo"],
                 "VERSION":conf["version"]}
  conf["PROGS"] = {"FASTQC"          : "fastqc",
                 "TRIMGALORE"      : "trim_galore",
                 "CUTADAPT"        : "cutadapt",
                 "BOWTIE2"         : "bowtie2" ,
                 "BISMARK"         : "bismark",
                 "DEDUPLICATE_BISMARK" : "deduplicate_bismark",
                 "BISMARK_GENOME_PREPARATION"     : "bismark_genome_preparation",
                 "BISMARK_METHYLATION_EXTRACTOR"  : "bismark_methylation_extractor",
                 "BISMARK2REPORT"                 : "bismark2report",
                 "SAMTOOLS"                       : "samtools"}
  conf["NICE"]="19"
  conf["RCODE"]=".read"
  conf["DIRECTIONAL"]=conf["directional"]
  conf["NUMTHREADS"]=conf['numthreads']
  conf["INEXT"]=".fq.gz"                                                  # TODO
  conf["FTP_SRA"]=conf["FTP_SRA"] if "FTP_SRA" in conf.keys() else None   # TODO
  #
  return(conf)                 
                

def parseTable( path_table ):
  import csv
  sreader = csv.reader(open(path_table), delimiter=',')
  rows = [row for row in sreader]
  header = rows[0]
  minimal_header = ['Read1', 'Read2', 'Sample id', 'Read type', 'Treatment']
  if( header[:5] != minimal_header):
     raise Exception( "First columns of the input table have to be " + ",".join(minimal_header)+"." )
  sample_ids = [x[2] for x in rows[1:]]
  if ( len(set(sample_ids)) != len(sample_ids) ):
    raise Exception( "Column 'Sample id' has non-unique values." )
  # filter(None, my_list) removes empty elements from the list my_list
  list_units = [  list(filter(None,[x[0],x[1]])) for x in rows[1:]  ]
  if (  True in [len(x)==0 for x in list_units] ):
    raise Exception( "Each sample has to have entry in at least of the columns 'Read1' or 'Read2'." )
  return( [rows, sample_ids, list_units] )


def check_if_fastq(string):
  fq = ["fasta", "fastq", "fq"]
  if string.find(fq[0]) == -1 or string.find(fq[1]) == -1 or string.find(fq[2]) == -1:
    return True
  else:
    return False

