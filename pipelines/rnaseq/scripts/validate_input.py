import os
import csv
import yaml
import argparse
from glob import glob

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

def validate_config(config):
    # Check that all locations exist
    for loc in config['locations']:
        if (not loc == 'output-dir') and (not (os.path.isdir(config['locations'][loc]) or os.path.isfile(config['locations'][loc]))):
            raise Exception("ERROR: The following necessary directory/file does not exist: {} ({})".format(config['locations'][loc], loc))

    sample_sheet = read_sample_sheet(config['locations']['sample-sheet'])
    
    # Check if the required fields are found in the sample sheet
    required_fields = set(['name', 'reads', 'reads2', 'sample_type'])
    not_found = required_fields.difference(set(sample_sheet[0].keys()))
    if len(not_found) > 0:
        raise Exception("ERROR: Required field(s) {} could not be found in the sample sheet file '{}'".format(not_found, config['locations']['sample-sheet']))

    # Check that requested analyses make sense
    if 'DEanalyses' in config:
        for analysis in config['DEanalyses']:
            for group in config['DEanalyses'][analysis]['case_sample_groups'] .split(',') + config['DEanalyses'][analysis]['control_sample_groups'].split(','):
                group = group.strip() #remove any leading/trailing whitespaces in the sample group names
                if not any(row['sample_type'] == group for row in sample_sheet):
                    raise Exception('ERROR: no samples in sample sheet have sample type {}, specified in analysis {}.'.format(group, analysis))

    # Check that reads files exist; sample names are unique to each row; 
    samples = {}        

    for row in sample_sheet:
        sample = row['name']
        if sample in samples:
            raise Exception('ERROR: name "{}" is not unique. Replace it with a unique name in the sample_sheet.'.format(sample))
        else:
            samples[sample] = 1
                   
        filenames = [row['reads'], row['reads2']] if row['reads2'] else [row['reads']]
        for filename in filenames:
            fullpath = os.path.join(config['locations']['reads-dir'], filename)
            if not os.path.isfile(fullpath):
                raise Exception('ERROR: missing reads file: {}'.format(fullpath))

                    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config-file', required=True, help='Path of configuration file [settings.yaml]')
    parser.add_argument('-s', '--sample-sheet-file', required=True, help='Path of sample sheet [sample_sheet.csv]')
    args = parser.parse_args()

    config = read_config_file(args.config_file)
    config['locations']['sample-sheet'] = args.sample_sheet_file
    validate_config(config)