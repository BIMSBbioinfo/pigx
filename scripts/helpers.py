import os

def makelink( actual_file, link_location, linkname):
  if not os.path.isfile(actual_file):
    print("ERROR: linking to non-existent file: "+actual_file ); status = -1
  elif not os.path.isdir(link_location):
    print("ERROR: PATHOUT and/or subdirectory does not exist for linking: "+link_location); status = -1
  else:
    status = os.system( "ln -sn "+actual_file+"  "+link_location+linkname+" 2>/dev/null" )
    #-- create a link only if non previously exists. suppress warnings (silent) if it's already there.
  return status


def fq_suffix(filename):
  return any(filename.endswith(ext) for ext in [".fq", ".fastq", ".fasta"])


def is_zipped(filename):
  return any(filename.endswith(ext) for ext in [".gz", ".bz2"])


def splitext_fqgz(string):
  if is_zipped(string):
    string, zipext = os.path.splitext(string)
  else:
    zipext = ""

  if fq_suffix(string):
    base, ext = os.path.splitext(string)
    return (base, ext + zipext)
  else:
    print("ERROR: Input files are not fastq files!")

