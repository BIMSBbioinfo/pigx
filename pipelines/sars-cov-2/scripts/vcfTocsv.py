# 21/04/04
# Vic-Fabienne Schumann
# input: vcf file


import os
import re
import csv
import sys

def vcfTocsv(vcffile):
    
    ''' This functions take lofreq (!) output vcfs as input, parses the information given in it's "info" column
     and makes a csv file out of it so it can be loaded into an R file easily.'''
    
    filename = os.path.basename(vcffile)
    filename_parsed = os.path.splitext(filename)[0]
    dirname = os.path.dirname(vcffile)
    
    filebody = filename.split(".")[0]
    csvfilename = f"{filebody}.csv"
    
    os.chdir(dirname)
    
    with open(csvfilename, "w", newline='') as csvfile: # because i get the "has to have a write option error"
   
        fieldnames = ["Chromosome", "Pos", "Ref", "Var","DP", "AF", "SB", "DP4"] # DP4 is missing but don't know right now how to deal with it
    
        writer = csv.DictWriter(csvfile, # TODO: make this dynamic
                            fieldnames=fieldnames, 
                            delimiter=",")
        writer.writeheader()
    
        with open(filename,"r") as vcf:

            for line in vcf:
                if re.match("^NC_",line):
                    snv_info1 = re.split(r'\t+',line)
                    Chrom = snv_info1[0]
                    Pos = snv_info1[1]
                    Ref = snv_info1[3]
                    Var = snv_info1[4]

                    snv_info2 = re.split(r';',snv_info1[7])
                    for element in snv_info2:
                        if re.search(r'DP=',element):
                            DP = element.split("=")[1] # can this be concat. to fewer lines?
                        if re.search(r'AF=',element):
                            AF = element.split("=")[1] # this part feels to be a bit redundant
                        if re.search(r'SB=',element):
                            SB = element.split("=")[1]
                        if re.search(r'DP4=',element):
                            DP4 = element.split("=")[1]

                    writer.writerow({'Chromosome': Chrom,
                                     'Pos': Pos, 
                                     'Ref': Ref,
                                     'Var': Var,
                                     'DP': DP, 
                                     'AF': AF,
                                     'SB': SB,
                                     'DP4': DP4.split("\n")[0]}) # has exclamations in the output, maybe get rid of them, don't know how right now, also not super important right now
        
if __name__ == '__main__':
    vcffile = sys.argv[1]
    vcfTocsv(vcffile)

