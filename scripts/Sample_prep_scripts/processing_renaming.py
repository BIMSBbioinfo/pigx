# this is mainly to document the process of the renaming and merging of the sample files 

from Renaming_by_table import *

sample_dir_0707 = "/mnt/bimsb_local/data/210818_ww_berlin2021_0707"
sample_dir_0721 = "/mnt/bimsb_local/data/210808_ww_berlin2021_0721"
sample_dir_merged = "/mnt/bimsb_local/data/210818_ww_berlin2021_all_merged"
bash_script = "/home/vfs/PycharmProjects/Akalinlab_pathogenomics/pigx_sarscov2_ww/scripts/Sample_prep_scripts/merging.sh"
naming_table = "/home/vfs/PycharmProjects/Akalinlab_pathogenomics/pigx_sarscov2_ww/misc_pub_related/sample_names_mayjune.csv"

# renaming(sample_dir_0707, bash_script, naming_table)
# renaming(sample_dir_0721, bash_script, naming_table)

# manually removed all not-ww and all not-berlin files 
# sample_dir_0707: 44 files (fake) 42 files (real)
# sample_dir_0721: 24 files (missing the april samples) (fake and real)

# merging(sample_dir_0707, sample_dir_0721, sample_dir_merged)
# --> 22 files, 12 samples

# manually renamed  P939_320422_20-22hmix to P939_210422_20-22hmix because I think that was a excel typo 

sample_dir_april_data = "/mnt/bimsb_local/data/ww_berlin_complete"
merge_by_date(sample_dir_0707, sample_dir_april_data, sample_dir_merged)

