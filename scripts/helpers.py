import os

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

