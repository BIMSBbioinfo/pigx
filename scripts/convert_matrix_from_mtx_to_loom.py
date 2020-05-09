import sys
import loompy
import pandas as pd
import re
import numpy
import os
from scipy.io import mmread
import argparse

# ------------------------------------------------------------------ #
#find the complete list of gene ids from the gtf file
def get_gtf_gene_ids (gtf_file):
    gene_ids = []
    with open(gtf_file) as f:
        i = 0
        for line in f:
            m = re.search('gene_id "(.+?)"', line)
            if m != None:
                gene_ids.append(m.group(1))
    return(set(gene_ids))

# ------------------------------------------------------------------ #
#given a matrix, row_names, and col_names, find which row names don't exist in gene_ids
#and fill the matrix with rows of zeros for each missing gene
def fill_missing_genes_into_matrix(gene_ids, matrix, row_names, col_names):
    missing = list(set(gene_ids) - set(row_names))
    # append the matrix with rows of zeros
    if(len(missing) > 0):
        new_rows  = numpy.zeros((len(missing), len(col_names)), dtype = int)
        matrix    = numpy.vstack((matrix, new_rows))
        row_names = row_names.append(pd.Series(missing))

    row_names = pd.DataFrame({'gene_id' : row_names})
    return(matrix, row_names)

# ------------------------------------------------------------------ #
def dict_to_array(d):
    for i in d.keys():
        d[i] = numpy.asarray(d[i])
    return(d)
# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert STAR mtx to loom')
    # sample name
    parser.add_argument('--sample_id',         action="store", dest='sample_id')
    # location of the star solo output
    parser.add_argument('--input_dir',         action="store", dest="input_dir")
    # location of the gtf file
    parser.add_argument('--gtf_file',          action="store", dest="gtf_file")
    # names of star solo output directories
    parser.add_argument('--star_output_types_keys', action="store",
    dest="star_output_types_keys", nargs='+')
    # internal names for star solo output - exists because of velocity
    parser.add_argument('--star_output_types_vals', action="store", dest="star_output_types_vals", nargs='+')
    parser.add_argument('--output_file',       action="store", dest="output_file")
    parser.add_argument('--sample_sheet_file', action="store", dest="sample_sheet_file")
    parser.add_argument('--path_script', action="store", dest="PATH_SCRIPT")

    args = parser.parse_args()
    # -------------------------------------------------------------- #
    sample_id              = args.sample_id
    basepath               = args.input_dir
    gtf_file               = args.gtf_file
    star_output_types_keys = args.star_output_types_keys
    star_output_types_vals = args.star_output_types_vals
    output_file            = args.output_file
    sample_sheet_file      = args.sample_sheet_file

    # -------------------------------------------------------------- #
    # implements custom mmread function - enables multi column mm files
    sys.path.append(args.PATH_SCRIPT)
    import matrix_market_IO


    # -------------------------------------------------------------- #
    print("Parsing gene ids from gtf file",gtf_file)
    gene_ids = get_gtf_gene_ids(gtf_file)

    # -------------------------------------------------------------- #
    dge_dict = {}
    print("Reading input files ...")
    for index in range(len(star_output_types_keys)):

        star_output_type  = star_output_types_keys[index]
        star_output_value = star_output_types_vals[index]

        # reads in multi column mtx for velocity
        # column_to_read == 2 corresponds to spliced
        # column_to_read == 3 corresponds to unspliced
        column_to_read = 2
        if(star_output_value == 'Unspliced'):
            column_to_read = 3

        print(star_output_value)
        path_input    = os.path.join(basepath, star_output_type, 'raw')

        print('Features ...')
        path_genes    = os.path.join(path_input,'features.tsv')
        genes         = pd.read_csv(path_genes, sep='\t', header=None)
        genes.columns = ['gene_id','gene_id2']
        genes         = genes.drop(columns = ['gene_id2'])

        print('Barcode ...')
        path_barcode       = os.path.join(path_input,'barcodes.tsv')
        barcode            = pd.read_csv(path_barcode, sep='\t', header=None)
        barcode.columns    = ['cell_id']
        barcode['sample_name'] = sample_id
        barcode['index']   = range(barcode.shape[0])
        barcode['cell_id'] = barcode['sample_name'] + '_' + barcode['cell_id']

        sample_sheet = pd.read_csv(sample_sheet_file)
        sample_sheet = sample_sheet.drop(columns=['reads','barcode'])
        sample_sheet = sample_sheet.drop_duplicates()
        barcode      = barcode.merge(sample_sheet, on='sample_name')
        barcode      = barcode.sort_values(by='index')
        barcode      = barcode.drop(columns=['index'])


        # -------------------------------------------------------------- #
        # fills the GENE matrix with zeros for missing genes
        print('Matrix ...')
        matrix = matrix_market_IO.mmread_pigx(os.path.join(path_input, 'matrix.mtx'), column_to_read = column_to_read)
        # row_names_gene needs to be a named pandas object
        matrix_gene_umi, row_names_gene = fill_missing_genes_into_matrix(gene_ids, matrix.toarray(), genes['gene_id'], barcode['cell_id'])

        dge_dict[star_output_value] = {'matrix' : matrix_gene_umi, 'barcode' : barcode, 'genes' : row_names_gene}

    # ADD: test to check that the rownames and the column names correspond
    # between diferent matrices

    # -------------------------------------------------------------- #
    print('Creating loompy file')
    # extracts the exon count matrix
    matrix_list = {'' : dge_dict[star_output_types_vals[0]]['matrix']}

    # extracts the exon + intron count matrix
    # extracts the velocyto matrices
    for i in range(1,len(star_output_types_vals)):
        matrix_list[star_output_types_vals[i]] = dge_dict[star_output_types_vals[i]]['matrix']

    col_attrs = dge_dict[star_output_types_vals[0]]['barcode'].to_dict("list")
    col_attrs = dict_to_array(col_attrs)

    row_attrs = dge_dict[star_output_types_vals[0]]['genes'].to_dict("list")
    row_attrs = dict_to_array(row_attrs)

    loompy.create(output_file, matrix_list, row_attrs, col_attrs)
