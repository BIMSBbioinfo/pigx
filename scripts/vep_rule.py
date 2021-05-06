
import os
import re


# 1. run ensemble vep
rule vep:
    input: 
        vcf_file = os.path.join(VARIANTS_DIR, '{sample}_snv.vcf')
        species = "sars_cov_2" # TODO get from settings file?
    output: os.path.join(VEP_DIR, '{sample}_vep_sarscov2.txt') # TODO add VEP_DIR to snakefile.py globals
    log: 
    shell: 
      """
      vep --verbose --offline --DB_VERSION 101 --appris --biotype --buffer_size 5000 --check_existing\
      --distance 5000 --mane --protein --species {input.species} --symbol --transcript_version --tsl\
      --input_file {input.vcf_file} --output_file {output}
      """

# 2. parse the output so it can be used as input for the rmarkdown report
rule parse_vep_output: 
    input: os.path.join(VEP_DIR, '{sample}_vep_sarscov2.txt')
    output: os.path.join(VEP_DIR, '{sample}_vep_sarscov2_reportinput.txt')
    run: 
        header_line = None
        
        if not header_line:
          raise Exception("file has a format or does not has to be parsed anymore")
        
        # remove unnecessary lines 
        with open(input, 'r') as file:
            with open(output, 'w') as output_file:
                reader = csv.reader(file, delimiter="\t")
                writer = csv.writer(output_file, lineterminator='\n')
        
                # make new column SYMBOL for gene names
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