import sys
import loompy
import pandas


# convert a UMI matrix file (plain txt) into a loom file
def convertToLoom (input_file, output_filepath): 
    df = pandas.read_table(input_file, index_col=0, nrows = 10) #todo remove nrows after finished testing
    matrix = pandas.DataFrame.as_matrix(df)
    row_attrs = { "Genes": list(df.index) }
    col_attrs = { "Cells": list(df.columns) }
    loompy.create(output_filepath, matrix, row_attrs, col_attrs)

if __name__ == '__main__':
    
    input_file = sys.argv[1]
    output_filepath = sys.argv[2] 
    
    print("converting matrix file ", input_file ,"into .loom format")
    convertToLoom(input_file, output_filepath)    
