
def wget_file(file_path, file_url):
    import os
    """
    Download a file named 'filename' from a url 'file_url' to 'file_path'  
    """
    if(os.path.isfile(file_path)):
        print('Using locally cached version of a file found here:\n' + file_path)
    else:
        os.system("wget -O " + file_path + " " + file_url)
        
        
def url_ENA(id):
    return "http://www.ebi.ac.uk/ena/data/warehouse/filereport\?accession\="+id+"\&result=read_run\&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,library_strategy,library_selection,fastq_ftp"

                       
def SRA2fastq(exp_id, output_path, library_strategy=None):
    """function that extracts URL/ftp links to fastq files of runs from project/experiment/run id from the
       ENA database and returns a dict with key=name_of_file and value = link or list of links
    """    
    wget_file(output_path + exp_id+"_ENA.txt", url_ENA(exp_id))
    #
    lines = open(output_path + exp_id + "_ENA.txt", 'r').readlines()
    header = lines[0].rstrip().split("\t")
    lines = [line.rstrip().split("\t") for line in lines[1:]]
    #
    if library_strategy is not None:
      libstra = [lines[10] for line in lines]
      which.lib = libstra.index(library_strategy)
      lines = [lines[i] for i in which.lib]
    #
    d = {} # {"name_of_srr":"ftp link to srr"}
    for line in lines:
      run_acc = line[5]
      try:
        fq_link = line[12].split(";")
      except (TypeError, IndexError):
        print("No fastq file in the ENA database for " + run_acc+". This entry will be skipped.")
        continue
      d[run_acc] = fq_link
    return(d)


def filter_filename_no_ext(x):
    return( x.split("/")[-1].split(".")[0] )  


def filter_filename_with_ext(x):
    return( x.split("/")[-1] )   


def merge_two_dicts(x, y):
  """Given two dicts, merge them into a new dict as a shallow copy."""
  z = x.copy()
  z.update(y)
  return z


def check_if_fastq(string):
  fq = ["fasta", "fastq", "fq"]
  if string.find(fq[0]) == -1 or string.find(fq[1]) == -1 or string.find(fq[2]) == -1:
    return True
  else:
    return False



# Example
#a = get_dict_SRA_ftp(["SRR2075688"], "/data/akalin/kwreczy/my_fastq_files/")
def get_dict_SRA_ftp(sra_ids, path2ENAreport):
  """
  sra_ids a list of strings indicatinf SRA ids
  path2ENAreport a path to a directory where ENA reports will be downloaded
  """
  dict_fastq_links = {}
  for sra_id in sra_ids:
    # sra id shouldn't have fastq.gz extension
    if check_if_fastq(sra_id)==False:
      raise Exception( "SRA id shouldn't have an fastq.gz extension" )
    dicti = SRA2fastq(sra_id, path2ENAreport)
    dict_fastq_links = merge_two_dicts(dict_fastq_links, dicti)
  return(dict_fastq_links)  
  
    
def create_new_rows(rows1, dicti):
  """
  rows1 a list of list indicating input table sheet
  dicti 
  """
  new_rows = []
  for r in range(1,len(rows1)):
    row = rows1[r]
    if_fq = check_if_fastq(row[0]) # if there is SRA id then its in the first column.
    if if_fq is False:
      new_rows.append(row)
      continue
    #
    for k in dicti.keys():
      sra_id = k
      new_row = row
      #
      if len(dicti[sra_id]) == 2:
        # paired-end
        new_row[0] = filter_filename_with_ext( dicti[sra_id ][0] )
        new_row[1] = filter_filename_with_ext( dicti[sra_id ][1] )
        new_row[2] = filter_filename_no_ext( dicti[sra_id ][0] ).split("_")[0]
      elif len(dicti[sra_id]) == 1:
        # single-end
        new_row[0] = filter_filename_with_ext( dicti[sra_id ][0] )
        new_row[2] = filter_filename_no_ext( dicti[sra_id ][0] )
      else:
        print("Sth went wrong.")
      new_rows.append(new_row[:])  # [:] has to be here http://stackoverflow.com/questions/5280799/list-append-changing-all-elements-to-the-appended-item
  return(new_rows)
    
    
    
