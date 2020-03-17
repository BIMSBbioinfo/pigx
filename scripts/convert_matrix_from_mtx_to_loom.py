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
    parser.add_argument('--sample_id',         action="store", dest='sample_id')
    parser.add_argument('--input_dir',         action="store", dest="input_dir")
    parser.add_argument('--gtf_file',          action="store", dest="gtf_file")
    parser.add_argument('--star_output_types', action="store", dest="star_output_types", nargs='+')
    parser.add_argument('--output_file',       action="store", dest="output_file")
    parser.add_argument('--sample_sheet_file', action="store", dest="sample_sheet_file")

    args = parser.parse_args()
    # -------------------------------------------------------------- #
    sample_id         = args.sample_id
    basepath          = args.input_dir
    gtf_file          = args.gtf_file
    star_output_types = args.star_output_types
    output_file       = args.output_file
    sample_sheet_file = args.sample_sheet_file

    # -------------------------------------------------------------- #
    print("Parsing gene ids from gtf file",gtf_file)
    gene_ids = get_gtf_gene_ids(gtf_file)

    # -------------------------------------------------------------- #
    dge_list = {}
    print("Reading input files ...")
    for star_output_type in star_output_types:
        print(star_output_type)
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
        matrix_gene = mmread(os.path.join(path_input, 'matrix.mtx'))

        # row_names_gene needs to be a named pandas object
        matrix_gene_umi, row_names_gene = fill_missing_genes_into_matrix(gene_ids, matrix_gene.toarray(), genes['gene_id'], barcode['cell_id'])

        dge_list[star_output_type] = {'matrix' : matrix_gene, 'barcode' : barcode, 'genes' : row_names_gene}

    # ADD: test to check that the rownames and the column names correspond
    # between diferent matrices

    # -------------------------------------------------------------- #
    print('Creating loompy file')
    matrix_list = {'' : dge_list['Gene']['matrix']}

    # checks whether the Velocyto and GeneFull outputs
    # are in the STAR output folder
    if 'Velocyto' in set(star_output_types):
        matrix_list['unspliced'] = dge_list['Velocyto']['matrix']

    if 'GeneFull' in set(star_output_types):
        matrix_list['total'] = dge_list['GeneFull']['matrix']

    col_attrs = dge_list['Gene']['barcode'].to_dict("list")
    col_attrs = dict_to_array(col_attrs)

    row_attrs = dge_list['Gene']['genes'].to_dict("list")
    row_attrs = dict_to_array(row_attrs)

    loompy.create(output_file, matrix_list, row_attrs, col_attrs)
