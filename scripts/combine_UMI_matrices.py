import os
import sys
import loompy

if __name__ == '__main__':
    
    loomfiles = sys.argv[1:-1]
    output_filepath = sys.argv[-1]
    loompy.combine(loomfiles, output_filepath, key = 'Genes')
