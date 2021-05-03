## 21/04/20 - writing a sample space script 

import subprocess
import os

datadir = "/local/Projects/AAkalin_pathogenomics/data/wastewater_all"
sampledir = "/home/vschuma/pathogenomics/vpipe/akalinlab_pathogenomics/samples/wastewater_all"

#########
for file in datadir:
    samplename = os.path.basename(file).split("_R")[0]
    filename_clean = file.split("_001")
    filename_clean = ''.join(filename_clean)

    # make dir
    dir = ''.join([samplename, '/', 'raw_data'])
    subprocess.call('mkdir -vp {}'.format(dir), shell=True)

    # copy files from local to sample dirs
    clean_file_sample_dir = ''.join([dir, '/', os.path.basename(filename_clean)])
    subprocess.call('cp -v {} {}'.format(file, clean_file_sample_dir), shell=True)
