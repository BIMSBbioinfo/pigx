import itertools
import os
import yaml
import sys

# ---------------------------------------------------------------------------- #
def validate_config(settings_dict, sample_sheet_file):
    with open(sample_sheet_file, 'r') as stream:
        sample_sheet_dict = yaml.load(stream)

    message = ''
    message = check_settings(settings_dict, message)
    message = check_sample_sheet(sample_sheet_dict, settings_dict, message)

    if len(message) > 0:
        message = 'ERROR: Config file is not properly formated:\n' + message
        sys.exit(message)

    return(0)

# ---------------------------------------------------------------------------- #
# checks the proper structure of the settings file
# ---------------------------------------------------------------------------- #
def check_settings(settings_dict, message):

    # ---------------------------------------------------------------------------- #
    params_list  = ['locations','general','execution','tools']
    message = check_params(settings_dict, params_list, message)

    # ---------------------------------------------------------------------------- #
    # checks for index or genome specification
    locations_dict = settings_dict['locations']
    if (locations_dict['genome-file'] == None) and (locations_dict['index-dir'] == None):
        message = message + "\t" + "neither genome nor index are specified\n"

    # checks whether the locations files exist if they are specified
    for location in sorted(list(set(locations_dict.keys()))):
        message = check_file_exists(locations_dict, location, message)

    return(message)



# ---------------------------------------------------------------------------- #
# checks the proper structure of the sample_sheet file
# ---------------------------------------------------------------------------- #
def check_sample_sheet(sample_sheet_dict, settings_dict, message):

    from itertools import chain

    # ---------------------------------------------------------------------------- #
    # checks for proper top level config categories
    params_list  = ['samples','peak_calling','idr','hub','feature_combination']
    message = check_params(sample_sheet_dict, params_list, message)

    # ---------------------------------------------------------------------------- #
    # checks for sample specification
    sample_desciptors = set(['fastq','library'])
    if len(sample_sheet_dict['samples'].keys()) > 0:
        for samp in sample_sheet_dict['samples'].keys():
            if(len(set(sample_sheet_dict['samples'][samp]) - sample_desciptors) > 0):
                message = message + "\t" + "sample: " + samp + " contains unrecognized sample     descriptors\n"
            # checks whether the fastq file names are given as a list
            if(not isinstance(sample_sheet_dict['samples'][samp]['fastq'], list)):
                message = message + "\t" + "sample: " + samp + " fastq should be a yaml list\n"

    else:
        message = message + "\t" + "There are no input samples!\n"

    # ---------------------------------------------------------------------------- #
    # checks for ChIP and Cont specifications
    peak_calling_desciptors = ['ChIP','Cont']
    if 'peak_calling' in set(sample_sheet_dict.keys()):
        if len(sample_sheet_dict['peak_calling'].keys()) > 0:
            for samp in sample_sheet_dict['peak_calling'].keys():
                if sample_sheet_dict['peak_calling'][samp]['ChIP'] == None:
                    message = message + '\t' + samp + ": " + "ChIP not specified\n"

            # if sample_sheet_dict['peak_calling'][samp]['Cont'] == None:
#                 message = message + '\t' + samp + ": " + "Cont not specified\n"

    # checks for correspondence between peak calling and samples
        if(len(sample_sheet_dict['samples']) > 0 and len(sample_sheet_dict['peak_calling']) > 0):
            samples = list(sample_sheet_dict['samples'].keys())
            keys = list(sample_sheet_dict['peak_calling'].keys())
            peaks = [[sample_sheet_dict['peak_calling'][i]['ChIP'],
                      sample_sheet_dict['peak_calling'][i]['Cont']]  for i in keys]
            peaks = flatten(peaks)
            peaks = list(filter(None, peaks))
            samples_diff = (set(peaks) - set(samples))
            if len(samples_diff) > 0:
                message = message + "\tsome peak calling samples are not specified\n"


    # ------------------------------------------------------------------------ #
    # checks whether the idr samples correspond to the peaks_calling samples
    if 'idr' in set(sample_sheet_dict.keys()):
        if(len(sample_sheet_dict['peak_calling']) > 0 and len(sample_sheet_dict['idr']) > 0):
            peaks_calling = set(sample_sheet_dict['peak_calling'].keys())
            for i in sample_sheet_dict['idr'].keys():
                peaks_idr = set([sample_sheet_dict['idr'][i][j] for j in sample_sheet_dict['idr'][i].keys()])
                if len(peaks_idr - peaks_calling) > 0:
                    message = message + "\tIDR: " + i + " Contains samples not in peak calling\n"

    # ------------------------------------------------------------------------ #
    # checks for proper feature combination
    # This check is temporary. Once Check_sample_sheet_dict is updated, can be removed.
    if 'feature_combination' in set(sample_sheet_dict.keys()):
        if len(sample_sheet_dict['feature_combination']) > 0:
            feature_keys = sample_sheet_dict['feature_combination'].keys()
            samps = []
            if 'idr' in set(sample_sheet_dict.keys()):
                samps = samps + list(sample_sheet_dict['idr'].keys())

            if 'peak_calling' in set(sample_sheet_dict.keys()):
                samps = samps + list(sample_sheet_dict['peak_calling'].keys())

            samps = set(samps)

            for key in feature_keys:
                key_diff = len(set(sample_sheet_dict['feature_combination'][key])  - samps)
                if(key_diff > 0):
                    message = message + "\tfeature_combination contains unknown peak files"

    # ------------------------------------------------------------------------ #
    # checks whether the designated files exist
    message = check_sample_exists(sample_sheet_dict, settings_dict, message)

    # checks whether extend is a number
    # if not (is.number(config['params']['extend'])):
    #     message = message + "extend must be a number\n"

    return(message)


# ---------------------------------------------------------------------------- #
# checks whether the supplied file or directory exists in the settings file
def check_file_exists(locations_dict, file_name, message=''):
    import os
    from itertools import chain
    if not locations_dict[file_name] == None:
       dirfile = locations_dict[file_name]
       dir_ind = os.path.isfile(dirfile) or os.path.isdir(dirfile)
       if not dir_ind:
           message = message + "\t" + file_name + "is not a valid file\n"

    return(message)

# ---------------------------------------------------------------------------- #
def check_sample_exists(sample_sheet_dict, settings_dict, message=''):
    import os
    from itertools import chain

    # checks whether the fastq path is specified
    locations_dict = settings_dict['locations']
    if locations_dict['input-dir'] == None:
        message = message + "\tfastq input directory is not specified\n"
    else:
        input_dir = locations_dict['input-dir']
        if not sample_sheet_dict['samples'] == None:
            files = [sample_sheet_dict['samples'][s]['fastq'] for s in sample_sheet_dict['samples'].keys()]
            files = [([x] if isinstance(x,str) else x) for x in files]
            files = list(itertools.chain(*files))
            for file in files:
                if not os.path.isfile(os.path.join(input_dir, file)):
                    message = message + '\t'+ file + ": file does not exist\n"

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
