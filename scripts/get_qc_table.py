import sys
import re
import csv


def parsing_cov(coverage_file):
    with open(coverage_file, "r") as cov_file:
        reader = csv.reader(cov_file, delimiter="\t")

        for row in reader:
            for element in row:
                if not re.match("#", element):
                    aligned_reads = row[3]
                    total_cov = row[5]
                    mean_depth = row[6]
        return (aligned_reads, total_cov, mean_depth)


def parsing_amplicon(amplicon_file):
    with open(amplicon_file, "r") as amp_file:
        reader = csv.reader(amp_file, delimiter="\t")

        num_covered_amplicons = 0
        covered_amplicon = 0
        drop_out_amplicons = []
        full_drop_out = []
        half_covered = []

        for row in reader:

            amplicon = row[3].split("_")[1]

            if amplicon != str(covered_amplicon):
                covered_amplicon = amplicon
                num_covered_amplicons += 1

            if int(row[6]) == 0:
                drop_out_amplicons.append(amplicon)

        if len(drop_out_amplicons) != 0:
            for ampli in drop_out_amplicons:
                cnt = drop_out_amplicons.count(ampli)
                if cnt == 2:
                    full_drop_out.append(ampli)
                if cnt == 1:
                    half_covered.append(ampli)

    full_drop_out = list(set(full_drop_out))
    fully_covered_amplicons = (int(num_covered_amplicons) - len(drop_out_amplicons)/2)

    return (fully_covered_amplicons, half_covered, full_drop_out)


def make_report_input_csv(coverage_file, amplicon_file, output_file):
    coverage_info = parsing_cov(coverage_file)
    amplicon_info = parsing_amplicon(amplicon_file)
    print(coverage_info)
    print(amplicon_info)

    header = [
        "Total number of amplicons fully covered",
        "Amplicons partially covered",
        "Drop-out amplicons",
        "Total number aligned reads",
        "Coverage",
        "Mean depth",
    ]
    fields = [
        amplicon_info[0],
        amplicon_info[1],
        amplicon_info[2],
        coverage_info[0],
        coverage_info[1],
        coverage_info[2],
    ]

    with open(output_file, "a") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerow(fields)

if __name__ == '__main__':
    coverage_file = sys.argv[1]
    amplicon_file = sys.argv[2]
    output = sys.argv[3]
    make_report_input_csv(coverage_file, amplicon_file, output)
