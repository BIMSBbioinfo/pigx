import os
import sys
import subprocess


def renaming(sample_dir, bash_script, naming_table):
    os.chdir(sample_dir)
    print("Entering ", sample_dir)
    res = subprocess.call([bash_script, naming_table])
    print(res)


def merging(sample_dir):
    if os.path.isdir(sample_dir):
        for i in os.listdir(sample_dir):
            file1 = i.split("_R")[0]
            for j in os.listdir(sample_dir):
                file2 = i.split("_R")[0]
                if file1 + "_R1.fastq.gz" == file2:
                    new = file1.split(".", 2)[0] + "_merged.fastq.gz"
                    zcat_command = f'zcat {file1} {file2} > {new}'
                    res = subprocess.call(zcat_command)
                    print(zcat)
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
