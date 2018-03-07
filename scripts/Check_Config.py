import itertools
import os
import yaml
import sys
from itertools import chain
import xlrd
import csv






# ---------------------------------------------------------------------------- #
def validate_config(config):



    sample_sheet_dict = read_sample_sheet(config['locations']['sample-sheet'])

    message = ''
    message = check_settings(sample_sheet_dict, config, message)
    message = check_sample_sheet(sample_sheet_dict, config, message)

    if len(message) > 0:
        message = 'ERROR: Config file is not properly formated:\n' + message
        sys.exit(message)

    return(0)



def read_sample_sheet(path):
    with open(path, 'r') as fp:
        sample_sheet = [row for row in csv.DictReader(fp, skipinitialspace=True)]
    return sample_sheet


# ---------------------------------------------------------------------------- #
# checks the proper structure of the sample_sheet file
# ---------------------------------------------------------------------------- #
def check_sample_sheet(sample_sheet_dict, config, message):

    sample_desciptors = set(['SampleName', 'Read', 'Read2'])
    
    if len(sample_sheet_dict) > 0:
        not_found = sample_desciptors.difference(set(sample_sheet_dict[0].keys()))
        if(len(not_found) > 0):
            message = message + "\t required columns " + str(not_found) + " were not found in sample sheet: " + config['locations']['sample-sheet'] + "\n"
    else:
        message = message + "\t" + "There are no input samples!\n"

    return message


# ---------------------------------------------------------------------------- #
# checks the proper structure of the settings file
# ---------------------------------------------------------------------------- #
def check_settings(sample_sheet_dict, config, message):


    # ---------------------------------------------------------------------------- #
    # checks for proper top level config categories
    params_list = ['locations', 'general', 'execution', 'tools',
                   'peak_calling', 'idr', 'hub', 'feature_combination']
    message = check_params(config, params_list, message)

    # ---------------------------------------------------------------------------- #
    # checks for index or genome specification
    locations_dict = config['locations']
    if (locations_dict['genome-file'] == None) and (locations_dict['index-dir'] == None):
        message = message + "\t" + "neither genome nor index are specified\n"

    # checks whether the locations files exist if they are specified
    for location in sorted(list(set(locations_dict.keys()))):
        message = check_file_exists(locations_dict, location, message)

    return(message)

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

    # ------------------------------------------------------------------------ #
    # checks whether the designated files exist
    message = check_sample_exists(sample_sheet_dict, config, message)

    # checks whether extend is a number
    # if not (is.number(conig['general']['params']['export_bigwig']['extend'])):
    #     message = message + "extend must be a number\n"

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
    else:
        input_dir = locations_dict['input-dir']
        if sample_sheet_dict:
            files = [[sample_dict['Read'],
                      sample_dict['Read2']]
                     for sample_dict in
                     sample_sheet_dict]
            files = list(itertools.chain(*files))
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
# given a list of lists, returns a flattened version
def flatten(l):
    out = []
    for item in l:
        if isinstance(item, (list, tuple)):
            out.extend(flatten(item))
        else:
            out.append(item)
    return out
