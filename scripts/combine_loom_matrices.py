import os
import loompy
import sys
import numpy
import pandas
import argparse

# ------------------------------------------------------------------ #
# get list of files from a txt file that contains space separated list of file paths
def get_filepaths (file):
    filepaths = []
    with open(file, "r") as f:
        for line in f:
            filepaths.extend(line.split())
    return(filepaths)

# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert STAR mtx to loom')
    parser.add_argument('--input_files',  action="store", dest="input_files", nargs='+')
    parser.add_argument('--output_file', action="store", dest="output_file")

    args = parser.parse_args()

    # -------------------------------------------------------------- #
    input_files     = args.input_files
    output_filepath = args.output_file


    loompy.combine(input_files, output_filepath, key = 'gene_id')
