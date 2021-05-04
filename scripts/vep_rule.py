
import os
import sys
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
            with open(os.path.join(VEP_DIR, "../../tmp.txt"), 'w+') as clean_file:
        
                for num, line in enumerate(file,header_line):
                    if re.match("#Uploaded_variation", line):
                        header_line = num
                        clean_file.write(line)
                        break
        # TODO parse the info part of the vfc to be seperate columns