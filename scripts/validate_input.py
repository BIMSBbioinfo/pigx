import os
import csv
import yaml
import argparse
from glob import glob


def read_config_file(path):
    with open(path, 'rt') as infile:
        config = yaml.load(infile)
    return config

def validate_config(config):
    # --------------------------------------------------------------------------------- #
    # Check that all locations exist
    for loc in config['locations']:
        if (not loc == 'output-dir') and config['locations'][loc] and (not (os.path.isdir(config['locations'][loc]) or os.path.isfile(config['locations'][loc]))):
            raise Exception("ERROR: The following necessary directory/file does not exist: {} ({})".format(config['locations'][loc], loc))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config-file', required=True, help='Path of configuration file [settings.yaml]')
    parser.add_argument('-s', '--sample-sheet-file', required=True, help='Path of sample sheet [sample_sheet.csv]')
    args = parser.parse_args()

    config = read_config_file(args.config_file)
    config['locations']['sample-sheet'] = args.sample_sheet_file
    validate_config(config)
