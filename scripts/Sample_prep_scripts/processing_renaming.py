# this is mainly to document the process of the renaming and merging of the sample files 

from Renaming_by_table import *

sample_dir_0707 = "/mnt/bimsb_local/data/210808_ww_berlin2021_0707_fake"
sample_dir_0721 = "/mnt/bimsb_local/data/210808_ww_berlin2021_0721_fake"
sample_dir_merged = "/mnt/bimsb_local/data/210818_ww_berlin2021_merged_fake"
bash_script = "/home/vfs/PycharmProjects/Akalinlab_pathogenomics/pigx_sarscov2_ww/scripts/Sample_prep_scripts/merging.sh"
naming_table = "/home/vfs/PycharmProjects/Akalinlab_pathogenomics/pigx_sarscov2_ww/misc_pub_related/sample_names_mayjune.csv"

renaming(sample_dir_0707, bash_script, naming_table)
renaming(sample_dir_0721, bash_script, naming_table)

# manually removed all not-ww and all not-berlin files 
# sample_dir_0707: 44 files
# sample_dir_0721: 24 files (missing the april samples)

# merging(sample_dir_0707, sample_dir_0721, sample_dir_merged)

# manually renamed  P939_320422_20-22hmix to P939_210422_20-22hmix because I think that was a excel typo 

# sample_dir_april_data = "/mnt/bimsb_local/data/ww_berlin_complete_fake"
# merge_by_date(sample_dir_0707, sample_dir_april_data, sample_dir_merged)

