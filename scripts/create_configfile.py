#!/usr/bin/env python3.5

import os, sys, json


def parse_config_args(config_args):
  """
  Check if user provided a tablesheet
  """
  if ( config_args=={} ):
    raise Exception( """Missing argument indicating 'tablesheet'.
                      Use argument --config and then type paths to the arguments, e.g.:
                      snakemake -s Snakemake.py --config tablesheet=./TableSheet.csv""")
  #      
  tablesheet = config_args.get("tablesheet")
  if (tablesheet is None):
     raise Exception( """Missing argument indicating 'tablesheet'.
                      Use argument --config and then type paths to the arguments, e.g.:
                      snakemake -s Snakemake.py --config tablesheet=./TableSheet.csv""" )

  return(tablesheet)


def splitInputFile2separateFile(file):
  """ It splits file into 3 separate files: 1.txt, 2.txt, 3.txt
  using regex [[ .* ]].
  """
  import os
  cmd="awk '/\[\[.*\]\]/{g++} { print $0 > g"+'"'+'.txt"'+"}'" + " " + file
  os.system(cmd)
  

def parseGeneralParams2dict(file, skip=None):
  """lines are lines from the file with general parameters for snakemake, such as:
  PATHIN='./in'
  PATHOUT=''
  GENOMEPATH=''
  LOGS=''
  CHROM_INFO=''
  GENOME_VERSION=''
  bismark_args=''
  fastqc_args=''
  trim_galore_args=''
  It returns a dictionary in which keys are variables (e.g. PATHIN) 
  and values are given by a user.
  """
  text=open(file,'r').readlines()

  if skip is not None and skip>=1:
    text = text[skip:]
  
  # remove \n characters after arguments
  text = [x.rstrip() for x in text]
  # remove empty lines
  text = list(filter(None, text))
  # list of lists (variable, value)
  list_of_list =[x.replace("'", "").split("=") for x in text]
  # convert it to list of variables and list of values
  (keys,values) = list(map(list, zip(*list_of_list)))
  dict_params=dict(zip(keys, values))
  return(dict_params)
    
    
def getFilenames(mylist):
  if len(mylist)==2:
    return( [splitext_fqgz(mylist[0])[0], splitext_fqgz(mylist[1])[0]] )
  if len(mylist)==1:
    return( [splitext_fqgz(mylist[0])[0]] )
  else:
    raise Exception("Sth went wrong in getFilenames())")
    
    
def getExtension(mylist):
  if len(mylist)==2:
    return( [splitext_fqgz(mylist[0])[1], splitext_fqgz(mylist[1])[1]] )
  if len(mylist)==1:
    return( [splitext_fqgz(mylist[0])[1]] )
  else:
    raise Exception("Sth went wrong in getExtension())")   
    

def parseTable2dict( path_table, skip=None):
  """
  Parse csv table with information about samples, eg:
  
  Read1,Read2,SampleID,Read type,Treatment
  sampleB.pe1.fq.gz,sampleB.pe2.fq.gz,sampleB,WGBS,B,,
  pe1.single.fq.gz,,sampleB1,WGBS,B,,
  
  It returns a dictionary required for the config file.
  """
  import csv
  sreader = csv.reader(open(path_table), delimiter=',')
  rows = [row for row in sreader]
  
  # remove empty lines
  rows = list(filter(None, rows))
  
  if skip is not None and skip>=1:
    rows = [x for x in rows[skip:]]

  header = rows[0]
  minimal_header = ['Read1', 'Read2', 'SampleID', 'ReadType', 'Treatment']
  
  if( header[:5] != minimal_header):
     raise Exception( "First columns of the input table have to be " + ",".join(minimal_header)+"." )
  sample_ids = [x[2] for x in rows[1:]]
  if ( len(set(sample_ids)) != len(sample_ids) ):
    raise Exception( "Column 'Sample id' has non-unique values." )
  # filter(None, my_list) removes empty elements from the list my_list
  list_units = [  list(filter(None,[x[0],x[1]])) for x in rows[1:]  ]
  if (  True in [len(x)==0 for x in list_units] ):
    raise Exception( "Each sample has to have entry in at least of the columns 'Read1' or 'Read2'." )
  units = dict(zip(sample_ids, list_units))
  
  # Create a dictionary with all params, keys are samples ids
  outputdict={}
  for i in range(len(rows[1:])):
    row = rows[1:][i]
    sampleid_dict = {}
    for j in range(len(header[2:])):
      try :
        sampleid_dict[header[2:][j]] = row[2:][j]
      except IndexError:
        raise Exception( "Number of columns in row"+j+" doesn't match number of elements in header." )
      
    sampleid_dict.update(  { 'files' : units[row[2]] }  )
    sampleid_dict.update(  { 'fastq_name' : getFilenames(units[row[2]]) }  )
    sampleid_dict.update(  { 'fastq_ext' : getExtension(units[row[2]]) }  )
    outputdict[row[2]] = sampleid_dict
    
  outputdict1 = {}
  outputdict1['SAMPLES'] =  outputdict 
  return( outputdict1 )


def createConfigfile(tablesheet, outfile, *args):
  """
  Create a config file in the JSON format. Arguments:
  tablesheet indicates a file that contain general parameters,
             info. about samples, sample specific arguments,
             and differential methylation
  outfile is a path to config file
  *args    is an additional argument, a path to the JSON file 
           or a dictionary which content will be added to the config file
  """
  dir_tablesheet = os.path.dirname(tablesheet)
  args=args[0] # only 1 additional argument
  
  splitInputFile2separateFile(tablesheet)
  general_paramspath = "1.txt"
  pathtable = "2.txt"
  sample_params = "3.txt"
  
  # Load general parameters
  gen_params=parseGeneralParams2dict(general_paramspath, skip=1)

  # Load parameters specific to samples
  sample_params = parseTable2dict( pathtable, skip=1)
  sample_params = dict(sample_params)
  
  # Create a config file  
  config=gen_params
  config.update(sample_params)
  
  # Remove temporary files
  os.remove("1.txt")
  os.remove("2.txt")
  os.remove("3.txt")

  # Add additional args to the config file
  if isinstance(args,dict):
    # convert OrderedDict to dict
    args=dict(args)
    config.update( args )
    
  if isinstance(args, str):
    with open(args) as data_file:    
      progs=json.load(data_file)
      # convert OrderedDict to dict
      progs=dict(progs)
      config.update( progs )

  # Save the config file
  with open(outfile, 'w') as outfile:
    dumps = json.dumps(config, 
                        indent=4, sort_keys=True, 
                        separators=(",",": "), ensure_ascii=True)
    outfile.write(dumps)

  return(config)


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
  

def main(argv):
    """
    Without the main sentinel, the code would be executed even if the script were imported as a module.
    """
    if len(argv)>1:
      
      # input table sheet
      tablesheet = argv[1]
      # path to config file
      outfile = argv[2]
      
      # Save config file
      if len(argv)<=2:
        config = createConfigfile(tablesheet, outfile)
      elif len(argv)==3:
        # if there is 3rd argument that indicate paths to tools
        # in JSON format include it to the config file
        progsjson = argv[3]
        config = createConfigfile(tablesheet, outfile, progsjson)


if __name__ == "__main__":
    main(sys.argv)  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

