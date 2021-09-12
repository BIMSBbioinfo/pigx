import os
import csv 

def make_table_with_cov( coverage_dir, output_file ):
    
    if os.path.isdir(coverage_dir):
        for file in os.listdir(coverage_dir):
            coverage = 0
            if "_coverage.csv" in file:
                samplename = file.split("_cov")[0]
                
                with open(os.path.join(coverage_dir, file), "r") as cov_csv:
                    rdr = csv.reader(cov_csv, delimiter = "\t")
                    for row in rdr: 
                        if not "#" in row[0] and "NC" in row[0]:
                            coverage = row[5]
                
                fields = [samplename, coverage]
                
                with open( output_file, "a") as output: 
                    wrtr = csv.writer(output)
                    wrtr.writerow(fields)
                
#coverage_dir = "/mnt/beast/pathogenomics/pigx_ww/210822_ww_berlin_febjune_publication/output/coverage"
coverage_dir = "/mnt/bimsb_local/pigx_sarscov2_ww/210531_wastewaterARTICEW109_pigxsarscov2ww_default/output/coverage"
#output_file = "/home/vfs/PycharmProjects/Akalinlab_pathogenomics/pigx_sarscov2_ww/tests/coverage_gathering.csv"
output_file = "/home/vfs/PycharmProjects/Akalinlab_pathogenomics/pigx_sarscov2_ww/tests/coverage_gathering_oldARTIC_analysis.csv"

make_table_with_cov(coverage_dir, output_file)
