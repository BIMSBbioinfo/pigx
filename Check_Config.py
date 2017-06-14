import itertools
import os
# ---------------------------------------------------------------------------- #
# checks the proper structure of the config file
# ---------------------------------------------------------------------------- #
def check_config(config):

    from itertools import chain
    message = ''
    # checks for proper top level config categories
    params = ['genome','index','fastq','params','samples','peak_calling','idr','software']

    # ---------------------------------------------------------------------------- #
    params_diff = set(config.keys()) - set(params)
    if len(params_diff) > 0:
        message = message + "\t" + "config file contains unknown parameters\n"

    # ---------------------------------------------------------------------------- #
    # checks for index or genome specification
    if (config['genome']['fasta'] == None) and (config['index'] == None):
        message = message + "\t" + "neither genome nor index are specified\n"

    # ---------------------------------------------------------------------------- #
    # checks for sample specification
    sample_desciptors = set(['fastq','library'])
    if len(config['samples'].keys()) > 0:
        for samp in config['samples'].keys():
            if(len(set(config['samples'][samp]) - sample_desciptors) > 0):
                message = message + "\t" + "sample: " + samp + " contains unrecognized sample     descriptors\n"
            # checks whether the fastq file names are given as a list
            if(not isinstance(config['samples'][samp]['fastq'], list)):
                message = message + "\t" + "sample: " + samp + " fastq should be a yaml list\n"
            
    else:
        message = message + "\t" + "There are no input samples!\n"

    # ---------------------------------------------------------------------------- #
    # checks for ChIP and Cont specifications
    peak_calling_desciptors = ['ChIP','Cont']
    if len(config['peak_calling'].keys()) > 0:
        for samp in config['peak_calling'].keys():
            if config['peak_calling'][samp]['ChIP'] == None:
                message = message + '\t' + samp + ": " + "ChIP not specified\n"

            if config['peak_calling'][samp]['Cont'] == None:
                message = message + '\t' + samp + ": " + "Cont not specified\n"

    # checks for correspondence between peak calling and samples
    if(len(config['samples']) > 0 and len(config['peak_calling']) > 0):
        samples = list(config['samples'].keys())
        keys = list(config['peak_calling'].keys())
        peaks = [[config['peak_calling'][i]['ChIP'],
                  config['peak_calling'][i]['Cont']]  for i in keys]
        peaks = flatten(peaks)
        samples_diff = (set(peaks[0]) - set(samples))
        if len(samples_diff) > 0:
            message = message + "\tsome peak calling samples are not specified\n"


    # ---------------------------------------------------------------------------- #
    # checks whether the idr samples correspond to the peaks_calling samples
    if(len(config['peak_calling']) > 0 and len(config['idr']) > 0):
        peaks_calling = set(config['peak_calling'].keys())
        for i in config['idr'].keys():
            peaks_idr = set([config['idr'][i][j] for j in config['idr'][i].keys()])
            if len(peaks_idr - peaks_calling) > 0:
                message = message + "\tIDR: " + i + " Contains samples not in peak calling\n"


    # ---------------------------------------------------------------------------- #
    # checks whether the designated files exist
    message = check_file_exists(config, message)

    # checks whether extend is a number
    # if not (is.number(config['params']['extend'])):
    #     message = message + "extend must be a number\n"

    if len(message) > 0:
        message = 'ERROR: Config file is not properly formated:\n' + message
        print(message)
        return(1);

    return(0)

# ---------------------------------------------------------------------------- #
def check_file_exists(config, message=''):
    import os
    from itertools import chain
    import IPython
#     IPython.embed()

    # checks whether the genome fasta file exists if the file is specified
    if not config['genome']['fasta'] == None:
        if not os.path.isfile(config['genome']['fasta']):
            message = message + "\tgenome file is not a valid file\n"

    # checks whether the fastq path is specified
    if config['fastq'] == None:
        message = message + "\tfastq input directory is not specified\n"
    else:
        if not config['samples'] == None:
            files = [config['samples'][s]['fastq'] for s in config['samples'].keys()]
            files = [([x] if isinstance(x,str) else x) for x in files]
            files = list(itertools.chain(*files))
            for file in files:
                if not os.path.isfile(os.path.join(config['fastq'],file)):
                    message = message + '\t'+file + "file does not exist\n"

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



# ---------------------------------------------------------------------------- #
# checks the config file for validity
# config check
def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False
