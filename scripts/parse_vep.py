import re
import os
import csv
import sys

def parse_vep(vep_dir, input, output):
    header_line = None

    with open(input, 'r') as file:
        if header_line:
            raise Exception("file has not the correct format or does not has to be parsed anymore")

        with open(os.path.join(vep_dir, "tmp.txt"), 'w+') as clean_file:
            for num, line in enumerate(file, 1):
                if re.match("#Uploaded_variation", line):
                    header_line = num
                    clean_file.write(line[1:])
                    break
            for num, line in enumerate(file, header_line):
                clean_file.write(line)

    with open(os.path.join(vep_dir,'tmp.txt'), 'r') as file:
        with open(output, 'w') as output:
            reader = csv.reader(file, delimiter="\t")
            writer = csv.writer(output, lineterminator='\n')

            # make new column SYMBOL for gene names 
            # fixme: This is redundant with the making of the header line above, but works for now
            gene_column = []
            header = next(reader)
            header.append("SYMBOL")
            gene_column.append(header)

            for row in reader:  # fixme: I'm not sure if the following construct is ugly or not
                # get gene name 
                extra = row[13].split(";")
                for element in extra:
                    if re.search(r"SYMBOL=", element):
                        gene = element.split("=")[1:][0]
                        row.append(gene)
                        gene_column.append(row)

            writer.writerows(gene_column)
    # TODO: add step to clean "tmp.txt" 

if __name__ == '__main__':
    vep_dir = sys.argv[1]
    input = sys.argv[2]
    output = sys.argv[3]
    parse_vep(vep_dir, input, output)