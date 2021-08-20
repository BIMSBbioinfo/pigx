import os
import argparse
import sys
import subprocess

def renaming(dir, bash_script, naming_table):
    os.chdir(dir)
    res = subprocess.call(bash_script)

def merging(dir):
    if os.path.isdir(dir):
        for file1 in os.listdir(dir):
            for file2 in os.listdir(dir):
                if file1 == file2:
                    new = file1.split(".", 2)[0] + "_merged.fastq.gz"
                    zcat = subprocess.call("zcat ", file1, " ", file2, "> ", new)
    else:
        raise Exception("directory not found")

def sorting():
    pass 

def main():
    pass

if __name__ == '__main__':
    main()