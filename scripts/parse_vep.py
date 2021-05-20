import re
import os
import csv
import sys

def parse_vep(vep_dir, input, output):
    with open(input, 'r') as file:
        with open(output, 'w') as output:
            header = None
            for line in file:
                if line.startswith("#Uploaded_variation"):
                    header = line[1:].strip().split("\t")
                    break

            reader = csv.reader(file, delimiter="\t")
            writer = csv.writer(output, lineterminator='\n')

            # make new column SYMBOL for gene names 
            gene_column = []
            header.append("SYMBOL")
            gene_column.append(header)

            extra_column = header.index("Extra")

            for row in reader:
                # get gene name 
                extra = row[extra_column].split(";")
                for element in extra:
                    if re.search(r"SYMBOL=", element):
                        gene = element.split("=")[1:][0]
                        row.append(gene)
                        gene_column.append(row)

            writer.writerows(gene_column)

if __name__ == '__main__':
    vep_dir = sys.argv[1]
    input = sys.argv[2]
    output = sys.argv[3]
    parse_vep(vep_dir, input, output)
