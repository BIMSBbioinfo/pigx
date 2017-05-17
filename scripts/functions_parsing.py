
def parseConfig(config_args):
  #      
  if ( config_args=={} ):
    raise Exception( """Missing arguments indicating 'tablesheet', 'gtoolbox','in','out','genome_folder', 'log' and 'chrominfo' or 'config'.
                      Use argument --config and then type paths to the arguments, e.g.:
                      snakemake -s Snakemake.py --config tablesheet=./TableSheet.csv gtoolbox=./programs/ in=./my_fastq_files/ out=./my_output/ genome_folder=./Genomes/hg19/ log=./logs/ chrominfo=./hg19.chrom.sizes.ucsc.txt numthreads=2 genome_version=hg19""" )
  #      
  tablesheet = config_args.get("tablesheet")
  config = config_args.get("configfile")
  if (tablesheet is not None and config is not None):
     raise Exception( """Missing an argument indicating 'tablesheet' or 'config'.
                      Use argument --config and then type paths to the arguments, e.g.:
                      snakemake -s Snakemake.py --config tablesheet=./TableSheet.csv gtoolbox=./programs/ in=./my_fastq_files/ out=./my_output/ genome_folder=./Genomes/hg19/ log=./logs/  chrominfo=./hg19.chrom.sizes.ucsc.txt numthreads=2 genome_version=hg19
                      or
                      snakemake -s Snakemake.py --config config=myconfigfile.json""" )
  #
  if (config is not None):
    with open(config) as config_file:    
      config_args = json.load(config_file)
      return(config_args)
  #  
  gtoolbox = config_args.get("gtoolbox")
  pathin = config_args.get("in")
  pathout = config_args.get("out")
  logs = config_args.get("log")
  config_args["gtoolbox"] = gtoolbox if gtoolbox is not None else ""
  config_args["in"] = pathin if pathin is not None else "./"
  config_args["out"] = pathout if pathout is not None else "./"
  config_args["log"] = logs if logs is not None else "./"
  #
  version = config_args.get("genome_version")
  config_args["genome_version"] = version if version is not None else ""
  numthreads = config_args.get("numthreads")
  config_args["numthreads"] = numthreads if numthreads is not None else 2
  #
  return(config_args)

                

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
    
    
def splitext_fqgz(string):
  ext = string.split(".")[-2:]
  core = string.split(".")[:-2]
  ext = ".".join(ext)
  if check_if_fastq(ext)==False:
    print("Input files are not fastq files!!")
  core = ".".join(core)
  return (core, ext)
  
  
def uniq_inext(list_units_only_ext):
  ext = list(set(list_units_only_ext))
  if len(ext)!=1:
    print("Input fastq files have different suffixes!")
    exit()
  return(ext)
  
  
def remove_ext_from_units(list_of_lists):
  list_units_no_ext = []
  list_units_only_ext = []
  for x in list_of_lists:
    if len(x)==2:
      core = [ splitext_fqgz(x[0])[0], splitext_fqgz(x[1])[0] ]
      ext = [ splitext_fqgz(x[0])[1], splitext_fqgz(x[1])[1] ]
      list_units_no_ext.append(core)
      list_units_only_ext.extend(ext)
    if len(x)==1:
      core = [ splitext_fqgz(x[0])[0] ]
      ext =  [ splitext_fqgz(x[0])[1] ]
      list_units_no_ext.append(core)
      list_units_only_ext.extend(ext)

  #TODO:make it more flexible?
  ext = uniq_inext(list_units_only_ext)
    
  return ( list_units_no_ext, "."+ext[0]) 

    

