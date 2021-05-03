
# 2021/04/28
# Vic-Fabienne Schumann
# dependencies: ensemble-vep

# this script contains the functions and commands to run the guix ensemble-vep for sars-cov-2 and parse it's output
# so it can directly used for the variant reports

# check if database is downloaded and located in ~/.vep
# TODO: or does this will be handled by snakemake?
# NOTE: Input file should be given with the full path, otherwise I can't guarantee that it works.^^
# NOTE: parse_vep_cli_output() only works if the vep_report exists already. 
#       But in theory the other functions are executed earlier anyways

# !!!!
# prepare vep data base by downloading the sars-cov-2 assembly and annotation: 
# bash: "wget ftp://ftp.ensemblgenomes.org/pub/viruses/variation/indexed_vep_cache/sars_cov_2_vep_101_ASM985889v3.tar.gz"
# bash: "mv sars_cov_2*.tar.gz ~/.vep && cd ~/.vep && tar xf sars_cov_2*.tar.gz"

import os
import sys
import re
import subprocess

def check_vep_dir(vep_report_txt): 
    
    vep_dir = os.path.dirname(vep_report_txt)
    if not os.path.isdir(vep_dir): 
        
        # print("mkdir -v {}".format(vep_dir))
       subprocess.call("mkdir -v {}".format(vep_dir), shell=True)

def run_vep(vep_report, snv_file):
    
     if not os.path.isfile(vep_report):
        
        # print("vep --verbose --offline --DB_VERSION 101 --appris --biotype --buffer_size 5000 --check_existing "
        #                 "--distance 5000 --mane --protein --species sars_cov_2 --symbol --transcript_version --tsl "
        #                 "--input_file {} --output_file {}".format(snv_file,vep_report))
        
        subprocess.call("vep --verbose --offline --DB_VERSION 101 --appris --biotype --buffer_size 5000 --check_existing "
                        "--distance 5000 --mane --protein --species sars_cov_2 --symbol --transcript_version --tsl "
                        "--input_file {} --output_file {}".format(snv_file,vep_report_txt))
    

def prepare_vep_output(vep_report_txt):
    
    ''' checks if the first line of the vep output is the correct headerline, will parse the file
    accordingly if not '''
    
    filename = os.path.basename(vep_report_txt)
    
    with open(vep_report_txt, 'r') as vep_report:
        # check if it starts with "Uploaded_variation" in the first line 
        if not re.match("^Uploaded_variation", vep_report.readline()):
            print("some parsing to do")
            parse_vep_cli_output(vep_report_txt)
            
            # remove old file and rename new file with that name
            # print(f"rm -v {vep_report_txt} && mv -v tmp.txt {filename}")
            subprocess.call(f"rm {vep_report_txt} && mv tmp.txt {filename}")
        else:
            print("fine to go")
            
def parse_vep_cli_output(vep_report):
    
    ''' takes the txt output of ensemble vep cli as input and removes all info lines and lines starting with # so it can be 
    read in as a table with correct header by R '''
    
    vep_dir = os.path.dirname(os.path.abspath(vep_report))
    print(f"operating in dir: {vep_dir}")
    header_line = None
    
    with open(vep_report, 'r') as file:
        with open(os.path.join(vep_dir,"tmp.txt"), 'w+') as clean_file:
            for num, line in enumerate(file,1):
                if re.match("#Uploaded_variation",line):
                    header_line = num
                    clean_file.write(line[1:])
                    break
            if header_line is not None:
                for num, line in enumerate(file,header_line):
                    clean_file.write(line)  
            else:
                raise Exception("file has a format or does not has to be parsed anymore")

def main():
    
    vep_report_txt = sys.argv[1]
    snv_file = sys.argv[2]
    
    # check if there is a vep dir in the sample dir, or make one
    check_vep_dir(vep_report_txt)
     
    # check if there is a vep file, if not run ensemble vep
    run_vep(vep_report_txt,snv_file)
   
    # read file and check if it looks ok or not
    prepare_vep_output(vep_report_txt)
    
if __name__ == '__main__':
  main()

