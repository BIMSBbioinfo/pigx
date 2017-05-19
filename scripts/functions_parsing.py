
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


def parseGeneralParams2dict(lines):
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
  It returns a dictionary in which keys are variables (e.g. PATHIN) and values are e.g. ('./in')
  """
  text=lines
  # remove empty lines
  text = list(filter(None, text))
  # remove \n and \t characters after arguments
  text = [x.rstrip('\r\n')for x in text]
  # list of tuples (variable, value)
  list_of_tuples = [(x.split("=")[0], x.split("=")[1].replace("'", "")) for x in text[1:]]
  # convert it to list of variables and list of values
  (keys,values) = list(map(list, zip(*list_of_tuples)))
  dict_params=dict(zip(keys, values))
  return(dict_params)
                

def parseTable2dict( path_table, skip=None):
  """
  Parse csv table with information about samples, eg:
  
  Read1,Read2,SampleID,Read type,Treatment
  sampleB.pe1.fq.gz,sampleB.pe2.fq.gz,sampleB,WGBS,B,,
  pe1.single.fq.gz,,sampleB1,WGBS,B,,
  
  It returns lines fo the file, sample ids (3rd column), 
  dictionary where keys are samples ids and values units (Read1 and Read1 columns)
  and other columns
  """
  import csv
  sreader = csv.reader(open(path_table), delimiter=',')
  rows = [row for row in sreader]
  
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
      
    sampleid_dict.update(  { 'fastq' : units[row[2]] }  )
    #sampleid_dict.update(  { 'inext' : units[row[2]] }  )
    #sampleid_dict.update(  { 'units' : units[row[2]] }  )
    outputdict[row[2]] = sampleid_dict
    
  outputdict1 = {}
  outputdict1['SAMPLES'] =  outputdict 
  return( outputdict1 )


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


def getFilenames(mylist):
  if len(mylist)==2:
    return( [splitext_fqgz(mylist[0])[0], splitext_fqgz(mylist[1])[0]] )
  if len(mylist)==1:
    return( [splitext_fqgz(mylist[0])[0]] )
  else:
    raise Exception("Sth went wrong in getFilenames())")

#  TODO: propably can be just removed   
# def remove_ext_from_units(list_of_lists):
#   list_units_no_ext = []
#   list_units_only_ext = []
#   for x in list_of_lists:
#     if len(x)==2:
#       core = [ splitext_fqgz(x[0])[0], splitext_fqgz(x[1])[0] ]
#       ext = [ splitext_fqgz(x[0])[1], splitext_fqgz(x[1])[1] ]
#       list_units_no_ext.append(core)
#       list_units_only_ext.extend(ext)
#     if len(x)==1:
#       core = [ splitext_fqgz(x[0])[0] ]
#       ext =  [ splitext_fqgz(x[0])[1] ]
#       list_units_no_ext.append(core)
#       list_units_only_ext.extend(ext)
# 
#   #TODO:make it more flexible?
#   ext = uniq_inext(list_units_only_ext)
#     
#   return ( list_units_no_ext, "."+ext[0]) 
# 
#     
# 
