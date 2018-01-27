import os
import loompy
import sys
import numpy
import pandas

def get_filepaths (file):
    filepaths = []
    with open(file, "r") as f:
        for line in f:
            filepaths.extend(line.split())
    return(filepaths)

if __name__ == '__main__':
    
    input_file = sys.argv[1]
    output_filepath = sys.argv[2]
    
    loomfiles = get_filepaths(input_file)
    
    loompy.combine(loomfiles, output_filepath, key = 'Genes')
