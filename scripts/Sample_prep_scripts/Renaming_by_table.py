import os
import sys
import subprocess


def renaming(sample_dir, bash_script, naming_table):
    """ takes a csv table with 3 cols and renames sample_names based in it's ID to it's sample name accordingly """

    if os.path.isdir(sample_dir):
        if len(os.listdir(sample_dir)) == 0:
            raise Exception("Directory is empty")
        else:   
            os.chdir(sample_dir)
            print("Entering ", sample_dir)
            res = subprocess.call([bash_script, naming_table])
            print(res)
    else:
        raise Exception("Directory does not exist")


def merging(sample_dir_1, sample_dir_2, target_dir):
    """ works with gz-zipped files that have the same name but are in different directories.
    Files with the same name will be merged, compressed and placed in a new target directory"""

    if os.path.isdir(sample_dir_1) and os.path.isdir(sample_dir_2):
        for file1 in os.listdir(sample_dir_1): 
            for file2 in os.listdir(sample_dir_2):
                if file1 == file2:
                    new = os.path.join(target_dir, ("merged_" + file1))
                    zcat_command = f'cat {os.path.join(sample_dir_1, file1)} {os.path.join(sample_dir_2, file2)} > {new}'
                    res = os.system(zcat_command)
                    print(zcat_command)
    else:
        raise Exception("directory not found")

def name_split(str):
    """ depending from which dir, splits sample name so the date/time stemp only is left """
    namebody = ""
    if str.startswith('P939'):
        namebody = str.split("39_")[1].split("mix_R")[0]
    elif str.startswith('WW'):
        namebody = str.split("W_")[1].split("_R")[0]
        
    return namebody

def merge_by_date(sample_dir_1, sample_dir_2, target_dir):
    """ works with gz-zipped files that have the same name but are in different directories.
    Files with the same date/time stemp (or whatever namebody they have) will be merged, compressed 
    and placed in a new target directory"""

    if os.path.isdir(sample_dir_1) and os.path.isdir(sample_dir_2):
        for file1 in os.listdir(sample_dir_1):
            namebody_1 = name_split(file1)
            suffix_1 = file1.split("_R")[1]
            
            for file2 in os.listdir(sample_dir_2):
                namebody_2 = name_split(file2)
                suffix_2 = file2.split("_R")[1]
                
                if namebody_1 == namebody_2 and suffix_1 == suffix_2:
                    new = os.path.join(target_dir, ("merged_"+file1))
                    zcat_command = f'cat {os.path.join(sample_dir_1, file1)} {os.path.join(sample_dir_2, file2)} > {new}'
                    res = os.system(zcat_command)
                    print(zcat_command)
    else:
        raise Exception("directory not found")


# not functional as it is right now, use the single functions
def main():
    sample_dir = sys.argv[0] # todo should be list of dirs
    bash_script = sys.argv[1]
    naming_table = sys.argv[2]

    # todo loop through list of input dirs
    renaming(dir, bash_script, naming_table) 
    merging(dir)


if __name__ == '__main__':
    # main()
    print("main is currently not functional , pls use the single functions")