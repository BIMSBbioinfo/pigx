# get a sample sheet with the basic information from a list of files 
# this version makes a file where the sample name will match the read name
# this will only create the basic skeleton with sample names till time, (do everything else on a spreadsheet for now)

import os
import csv
import sys
import re

def build_csv(csvfile, sample_dir):
    with open(csvfile, "w", newline='') as csvfile:
        fieldnames = ["name", "reads", "reads2", "date", "location_name", "coordinates_lat", "coordinates_long"]
        writer = csv.DictWriter(csvfile,
                                fieldnames=fieldnames,
                                delimiter=",")
        writer.writeheader()
        
        for file in os.listdir(sample_dir):
            if "_R" in file: 
                print("files are paired end")
                name = file.split("_R")[0]
                read1 = [file for file in os.listdir(sample_dir) if re.match(f"{name}_R1]*", file)][0]
                read2 = [file for file in os.listdir(sample_dir) if re.match(f"{name}_R2]*", file)][0] # can be empty in case of single end
                # TODO:
                # date = ""
                writer.writerow({'name': name,
                                 'reads': read1,
                                 'reads2': read2})
            else: 
                raise Exception("The files seem to be not paired-end reads. This script can't deal with that for now."
                                "You could change the read names to match the pattern of '*_R1.fastq.gz'.")
            
# this is here because WIP 
csvfile = "/home/vfs/PycharmProjects/Akalinlab_pathogenomics/pigx_sarscov2_ww/tests/auto_sample_sheet.csv"
sample_dir = "/mnt/bimsb_local/data/210818_ww_berlin2021_merged_fake"
build_csv(csvfile, sample_dir)
