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

def parsing_mut_loc_coverage(muts_loc_coverage):
    ''' takes the output of samtools bedcov with 5 columns and parses how many locations are covered,
    how much and sums this up'''
    with open(muts_loc_coverage, "r") as muts_loc_coverage:
        reader = csv.reader(muts_loc_coverage, delimiter="\t")

        num_covered_mutations_locations = 0
        covered_mutation_locations = 0
        drop_out_mutation_locations = []

        for row in reader:

            mut_location = row[3].split("_")[1]
            coverage_value = 0

            if mut_location != str(covered_mutation_locations):
                # count total amount of locations to check
                num_covered_mutations_locations += 1
                # track coverage per base
                coverage_value = int(row[4]) / 10  # reported sum per base, we use 10 base range
                # track muations with no coverage
                if coverage_value == 0:
                    drop_out_mutation_locations.append(mut_location)

    # get number of actual covered locations by subtracting 0 covered locations from total checked locations
    covered_mutation_locations = (int(num_covered_mutations_locations) - len(drop_out_mutation_locations))

    return (num_covered_mutations_locations, covered_mutation_locations, drop_out_mutation_locations)

def make_report_input_csv(genome_coverage_file, muts_loc_coverage_file, output_file):
    genome_coverage_info = parsing_cov(genome_coverage_file)
    mut_loc_cov_info = parsing_mut_loc_coverage(muts_loc_coverage_file)
    print(genome_coverage_info)
    print(mut_loc_cov_info)

    header = [
        "Total number of tracked mutations",
        "Total number of mutations covered",
        "Number of mutations not covered",
        "Total number aligned reads",
        "Percentage ref.genome covered",
        "Mean depth ref.genome coverage",
    ]
    fields = [
        mut_loc_cov_info[0],
        mut_loc_cov_info[1],
        mut_loc_cov_info[2],
        genome_coverage_info[0],
        genome_coverage_info[1],
        genome_coverage_info[2],
    ]

    with open(output_file, "a") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerow(fields)

if __name__ == '__main__':
    genome_coverage_file = sys.argv[1]
    muts_loc_coverage_file = sys.argv[2]
    output = sys.argv[3]
    make_report_input_csv(genome_coverage_file, muts_loc_coverage_file, output)
