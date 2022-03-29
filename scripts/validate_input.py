import os
import csv
import yaml
import argparse
from glob import glob
import itertools

def read_sample_sheet(path):
    with open(path, 'r') as fp:
        rows =  [row for row in csv.reader(fp, delimiter=',')]
        header = rows[0]; rows = rows[1:]
        sample_sheet = [dict(zip(header, row)) for row in rows]
    return sample_sheet

def read_config_file(path):
    with open(path, 'rt') as infile:
        config = yaml.load(infile)
    return config

def validate_memory_restrictions(config):
    for rule in config['execution']['rules']:
        if 'memory' in config['execution']['rules'][rule]:
            value = config['execution']['rules'][rule]['memory']
            if not type(value) == int and not value.isdigit():
                raise Exception("ERROR: memory limits must be expressed as a plain number of megabytes.  Got '{}' in '{}'.".format(value, rule))

def validate_config(config):
    # Check that all locations exist
    for loc in config['locations']:
        if (not loc == 'output-dir') and (not (os.path.isdir(config['locations'][loc]) or os.path.isfile(config['locations'][loc]))):
            raise Exception("ERROR: The following necessary directory/file does not exist: '{}' ({})".format(config['locations'][loc], loc))
        if not loc == 'output-dir':
            if config['locations'][loc].endswith(".gz") or config['locations'][loc].endswith(".bz2") or config['locations'][loc].endswith(".xz"):
                raise Exception("ERROR: The {} file '{}' is referenced in its compressed form like it was downloaded. However, the tools of this workflow expects the reference as plain FASTA/GTF files that are directly readable. Please unpack these files and update the settings to the respective new filename. This does _not_ hold for the FASTQ.gz files of the samples which shall remain compressed.".format(loc,config['locations'][loc]))

    validate_memory_restrictions(config)

    sample_sheet = read_sample_sheet(config['locations']['sample-sheet'])
    
    # Check if the required fields are found in the sample sheet
    # also add the list of covariates from the DE analyses
    covariates = [config['DEanalyses'][x]['covariates'] for x in config['DEanalyses'].keys()]
    # cleanup and get the set of unique covariates
    covariates = [y.strip() for y  in itertools.chain(*[x.split(',') for x in covariates])]
    # remove empty strings
    covariates = [x for x in covariates if x]
    required_fields = set(['name', 'reads', 'reads2', 'sample_type'] + covariates)
    not_found = required_fields.difference(set(sample_sheet[0].keys()))
    if len(not_found) > 0:
        raise Exception("ERROR: Required field(s) {} could not be found in the sample sheet file '{}'".format(not_found, config['locations']['sample-sheet']))

    # Check that requested analyses make sense
    if 'DEanalyses' in config:
        for analysis in config['DEanalyses']:
            for group in config['DEanalyses'][analysis]['case_sample_groups'] .split(',') + config['DEanalyses'][analysis]['control_sample_groups'].split(','):
                group = group.strip() #remove any leading/trailing whitespaces in the sample group names
                if not any(row['sample_type'] == group for row in sample_sheet):
                    raise Exception("ERROR: no samples in sample sheet have sample type '{}', specified in analysis {}.".format(group, analysis))

    # Check that reads files exist; sample names are unique to each row; 
    samples = {}        
    for row in sample_sheet:
        sample = row['name']
        if sample in samples:
            raise Exception('ERROR: name "{}" is not unique. Replace it with a unique name in the sample_sheet.'.format(sample))
        else:
            samples[sample] = 1
        filenames = [row['reads'], row['reads2']] if row['reads2'] else [row['reads']]
        # don't allow paths in the sample sheet, only allow the base name of the file
        if sum([os.path.basename(x) != x for x in filenames]) > 0:
            raise Exception(" ".join(["ERROR: read file names in the sample sheet must be basenames", 
            "not paths to the files. See sample ", sample, "in the sample sheet"]))
        for filename in filenames:
            fullpath = os.path.join(config['locations']['reads-dir'], filename)
            if not os.path.isfile(fullpath):
                filenameFlankedWithWhitespace = filename.startswith(' ')  or filename.endswith(' ') or filename.startswith('\t') or filename.endswith('\t')
                if filenameFlankedWithWhitespace:
                    raise Exception("ERROR: missing reads file: '{}', likely caused by blanks flanking the filename, please correct.".format(fullpath))
                else:
                    raise Exception("ERROR: missing reads file: '{}'".format(fullpath))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config-file', required=True, help='Path of configuration file [settings.yaml]')
    parser.add_argument('-s', '--sample-sheet-file', required=True, help='Path of sample sheet [sample_sheet.csv]')
    args = parser.parse_args()

    config = read_config_file(args.config_file)
    config['locations']['sample-sheet'] = args.sample_sheet_file
    validate_config(config)
