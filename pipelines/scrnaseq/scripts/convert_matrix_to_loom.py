import sys
import loompy
import pandas
import re
import numpy 

#from BCBio import GFF

#find the complete list of gene ids from the gtf file
def get_all_gene_ids (gtf_file):
    gene_ids = []
    with open(gtf_file) as f:
        i = 0
        for line in f:
            m = re.search('gene_id "(.+?)"', line)
            if m != None:
                gene_ids.append(m.group(1))
    return(set(gene_ids))
            
# parse the gene x cell matrix file and return a matrix object, row attributes along with column attributes (necessary to construct a loom file)
def parse_umi_matrix_file (matrix_file):
    df = pandas.read_table(input_file, index_col=0)
    matrix = pandas.DataFrame.as_matrix(df)
    row_names = list(df.index)
    col_names = list(df.columns)
    return(matrix, row_names, col_names)

#given a matrix, row_names, and col_names, find which row names don't exist in gene_ids 
#and fill the matrix with rows of zeros for each missing gene
def fill_missing_genes_into_matrix(gene_ids, matrix, row_names, col_names):
    missing = list(set(gene_ids) - set(row_names))
    print("Found ",len(missing),"genes missing in the given matrix with",len(row_names),"genes")
    # append the matrix with rows of zeros 
    new_rows = numpy.zeros((len(missing), len(col_names)), dtype = int)
    print("appending a matrix of zeros of dimensions",new_rows.shape)
    matrix = numpy.vstack((matrix, new_rows))
    row_names = row_names + missing
    return(matrix, row_names)

if __name__ == '__main__':
    
    sample_id = sys.argv[1]
    input_file = sys.argv[2]
    gtf_file = sys.argv[3]
    output_filepath = sys.argv[4]
                    
    print("Parsing umi matrix file",input_file)
    umi_matrix, row_names, col_names = parse_umi_matrix_file(input_file)
    print("umi matrix shape",umi_matrix.shape)
    print("rows:",len(row_names))
    print("cols:",len(col_names))
    
    print("Parsing gene ids from gtf file",gtf_file)
    gene_ids = get_all_gene_ids(gtf_file)
    print("Found",len(gene_ids),"genes in the gtf file")
    
    umi_matrix, row_names = fill_missing_genes_into_matrix(gene_ids, umi_matrix, row_names, col_names)
    
    print("Dimensions of the new matrix",umi_matrix.shape)
    
    print("Converting matrix of shape",umi_matrix.shape,"into file",output_filepath)
    row_attrs = { "Genes": row_names}
    col_attrs = { "cell_id": ['_'.join([sample_id, cn]) for cn in col_names] } #prepend sample ids to the cell ids
    loompy.create(output_filepath, umi_matrix, row_attrs, col_attrs)
    
    