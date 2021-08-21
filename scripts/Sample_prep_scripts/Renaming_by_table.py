import os
import sys
import subprocess


def renaming(sample_dir, bash_script, naming_table):
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
    if os.path.isdir(sample_dir_1) && os.path.isdir(sample_dir_2):
        for file1 in os.listdir(sample_dir_1): 
            for file2 in os.listdir(sample_dir_2):
                if file1 == file2:
                    new = os.path.join(target_dir, file1.split(".", 2)[0], "_merged.fastq.gz")
                    zcat_command = f'zcat {file1} {file2} > {new}'
                    res = subprocess.call(zcat_command)
                    print(zcat_command)
    else:
        raise Exception("directory not found")


def main():
    sample_dir = sys.argv[0]
    bash_script = sys.argv[1]
    naming_table = sys.argv[2]

    renaming(dir, bash_script, naming_table)
    merging(dir)


if __name__ == '__main__':
    main()
