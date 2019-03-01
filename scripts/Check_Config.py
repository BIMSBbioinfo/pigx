import itertools
import os
import yaml
import sys
from itertools import chain
import xlrd
import csv

# ---------------------------------------------------------------------------- #
# config dict of dicts
# sample_sheet_dict - list of dicts
# sample_sheet_column_names - list of allowed column names for sample sheet
def validate_config(config, structure_variables):

    sample_sheet_dict = read_SAMPLE_SHEET(config)

    message = ''
    message = check_sample_sheet(sample_sheet_dict, config,  structure_variables, message)
    if len(message) > 0:
        message = 'ERROR: Sample Sheet is not properly formated:\n' + message
        sys.exit(message)

    message = check_settings(sample_sheet_dict, config, structure_variables , message)
    if len(message) > 0:
        message = 'ERROR: Config file is not properly formated:\n' + message
        sys.exit(message)

    return(0)


# ---------------------------------------------------------------------------- #
# checks the proper structure of the sample_sheet file
# ---------------------------------------------------------------------------- #
def check_sample_sheet(sample_sheet_dict, config, structure_variables, message):

    sample_sheet_column_names = set(structure_variables['SAMPLE_SHEET_COLUMN_NAMES'])
    if len(sample_sheet_dict) > 0:
        not_found = sample_sheet_column_names.difference(set(sample_sheet_dict[0].keys()))
        if(len(not_found) > 0):
            message = message + "\t required columns " + str(not_found) + " were not found in sample sheet: " + config['locations']['sample-sheet'] + "\n"
    else:
        message = message + "\t" + "There are no input samples!\n"

    return message


# ---------------------------------------------------------------------------- #
# checks the proper structure of the settings file
# ---------------------------------------------------------------------------- #
def check_settings(sample_sheet_dict, config, structure_variables, message):

    # ---------------------------------------------------------------------------- #
    # checks for proper top level config categories
    message = check_params(config, structure_variables['SETTING_SUBSECTIONS'], message)

    # ---------------------------------------------------------------------------- #
    # sets obligatory genome files
    # checks for index or genome specification - obligatory
    locations_dict = config['locations']
    for obligatory_file in structure_variables['OBLIGATORY_FILES']:
        if locations_dict[obligatory_file] == None:
            message = message + "\t" + obligatory_file + " is not specified\n"


    # checks whether the locations files exist if they are specified
    for location in sorted(list(set(locations_dict.keys()))):
        message = check_file_exists(locations_dict, location, message)


    # ---------------------------------------------------------------------------- #
    # checks whether the fasta header contains whitespaces
    message = check_fasta_header(locations_dict['genome-file'], message)

    # ---------------------------------------------------------------------------- #
    # checks whether the general section is properly formatted
    message = check_general_section(config['general'], structure_variables, message)

    # ---------------------------------------------------------------------------- #
    # checks for ChIP and Cont specifications
    peak_calling_desciptors = ['ChIP', 'Cont']
    if 'peak_calling' in set(config.keys()):
        if len(config['peak_calling'].keys()) > 0:
            for samp in config['peak_calling'].keys():
                if config['peak_calling'][samp]['ChIP'] == None:
                    message = message + '\t' + samp + ": " + "ChIP not specified\n"

            # if sample_sheet_dict['peak_calling'][samp]['Cont'] == None:
#                 message = message + '\t' + samp + ": " + "Cont not specified\n"

    # checks for correspondence between peak calling and samples
        if(len(sample_sheet_dict) > 0 and len(config['peak_calling']) > 0):
            samples = [sample_dict['SampleName'] for sample_dict in sample_sheet_dict]
            keys = list(config['peak_calling'].keys())
            peaks = [[config['peak_calling'][i]['ChIP'],
                      config['peak_calling'][i]['Cont']] for i in keys]
            peaks = flatten(peaks)
            peaks = list(filter(None, peaks))
            samples_diff = (set(peaks) - set(samples))
            if len(samples_diff) > 0:
                message = message + "\tsome peak calling samples are not specified\n"


    # ------------------------------------------------------------------------ #
    # checks whether the idr samples correspond to the peaks_calling samples
    if 'idr' in set(config.keys()):
        if(len(config['peak_calling']) > 0 and len(config['idr']) > 0):
            peaks_calling = set(config['peak_calling'].keys())
            for i in config['idr'].keys():
                peaks_idr = set([config['idr'][i][j] for j in config['idr'][i].keys()])
                if len(peaks_idr - peaks_calling) > 0:
                    message = message + "\tIDR: " + i + " Contains samples not in peak calling\n"

    # ------------------------------------------------------------------------ #
    # checks for proper feature combination
    # This check is temporary. Once Check_sample_sheet_dict is updated, can be removed.
    if 'feature_combination' in set(config.keys()):
        if len(config['feature_combination']) > 0:
            feature_keys = config['feature_combination'].keys()
            samps = []
            if 'idr' in set(config.keys()):
                samps = samps + list(config['idr'].keys())

            if 'peak_calling' in set(config.keys()):
                samps = samps + list(config['peak_calling'].keys())

            samps = set(samps)

            for key in feature_keys:
                key_diff = len(set(config['feature_combination'][key])  - samps)
                if(key_diff > 0):
                    message = message + "\tfeature_combination contains unknown peak files"

    # ---------------------------------------------------------------------------- #
    # checks for correspondence between Spike-in in sample_sheet_dict and settings
    spike_in_indicator = False
    for sample in sample_sheet_dict:
        if 'Spike-in' in set(sample.keys()):
            if sample['Spike-in'].upper() == 'YES':
                spike_in_indicator = True

    # checks whether the spikein-file field is defined
    if spike_in_indicator:
        if not config['locations']['spikein-file']:
            message = message + 'spike-in normalisation is requested, but: '
            message = message + 'spike-in genome fasta file is not defined. Please set locations: spikein-file'            
        else:
            # checks whether the spikein-file fasta is a valid file
            message = check_file_exists(config['locations'],'spikein-file',message)
            # checks whether the spikein fasta file contains header with spaces
            message = check_fasta_header(locations_dict['spikein-file'], message)

    # ------------------------------------------------------------------------ #
    # checks whether the designated files exist
    message = check_sample_exists(sample_sheet_dict, config, message)


    return(message)


# ---------------------------------------------------------------------------- #
# checks whether the supplied file or directory exists in the settings file
def check_file_exists(locations_dict, file_name, message=''):

    if not locations_dict[file_name] == None:
       dirfile = locations_dict[file_name]
       dir_ind = os.path.isfile(dirfile) or os.path.isdir(dirfile)
       if not dir_ind:
           message = message + "\t" + file_name + " is not a valid file\n"

    return(message)

# ---------------------------------------------------------------------------- #
def check_sample_exists(sample_sheet_dict, config, message=''):

    # checks whether the fastq path is specified
    locations_dict = config['locations']
    if not locations_dict['input-dir']:
        message = message + "\tfastq input directory is not specified\n"

    elif not os.path.isdir(locations_dict['input-dir']):
        message = message + "\tfastq input directory does not exist\n"

    else:
        input_dir = locations_dict['input-dir']
        for sample_dict in sample_sheet_dict:
            files = []
            if sample_dict['library_type'] == 'single':
                files = [sample_dict['Read']]

            elif sample_dict['library_type'] == 'paired':
                files = [sample_dict['Read'], sample_dict['Read2']]

                for file in files:
                    if not os.path.isfile(os.path.join(input_dir, file)):
                        message = message + '\t' + file + ": file does not exist\n"
    return(message)

# ---------------------------------------------------------------------------- #
# given a dict, checks for existence of params
def check_params(config_dict, params, message):
    params_diff = list(set(config_dict.keys()) - set(params))
    params_str = " ".join(params_diff)
    if len(params_diff) > 0:
        message = message + "config file contains unknown parameters:" + params_str + "\n"

    return(message)

# ---------------------------------------------------------------------------- #
# checks fasta header for spaces
def check_fasta_header(genome_file, message):
    import re
    import sys
    import magic as mg
    import gzip

    genome_file_type = mg.from_file(genome_file, mime=True)
    if genome_file_type.find('gzip') > 0:
        file = gzip.open(genome_file, 'rb')
    else:
        file = open(genome_file, "rb")
    
    # get ascii for '>', ' ' and '\t'
    gt = ord('>')
    ws = ord(' ')
    tab= ord('\t')
    bar = ord('|')

    for line in file:
        if line[0] == gt:
            if (ws in line) or (tab in line) or (bar in line) : 
                message = message + 'Genome fasta headers contain whitespaces.\n Please reformat the headers\n'
                return(message)

    return(message)


# ---------------------------------------------------------------------------- #
# given a string checks if it is a number
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# ---------------------------------------------------------------------------- #
# checks general section for common mistakes
def check_general_section(general_dict, structure_variables, message):

    # checks wether assembly is set
    if not general_dict['analysis']:
        message = message + "GENERAL: analysis type needs to be supplied\n"
    else:
        if not general_dict['analysis'].upper() in structure_variables['ANALYSIS_TYPES']:
            message = message + "GENERAL: supported analysis types are:" + ','.join(sort(list(structure_variables['ANALYSIS_TYPES']))) + "\n"

    if not general_dict['assembly']:
        message = message + "GENERAL: genome assembly is not given\n"

    # checks whether bigwig extend is a number
    if not (is_number(general_dict['params']['export_bigwig']['extend'])):
        message = message + "GENERAL PARAMS: export_bigwig: extend must be a number\n"

    # checks whether mapq is a number
    if not (is_number(general_dict['params']['bam_filter']['mapq'])):
        message = message + "GENERAL PARAMS: bam_filter: mapq must be a number\n"

    # checks whether idr-threshold is a number
    if not (is_number(general_dict['params']['idr']['idr-threshold'])):
        message = message + "GENERAL PARAMS: idr: idr-threshold must be a number\n"

    # checks whether custom scripts get a number
    for param in ['extract_signal', 'peak_statistics', 'width_params']:
        for key in general_dict['params'][param].keys():
            if not (is_number(general_dict['params'][param][key])):
                message = message + "GENERAL PARAMS: {}: {} must be a number\n".format(param,key)

    return(message)
