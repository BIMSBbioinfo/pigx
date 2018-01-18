# ---------------------------------------------------------------------------- #
import subprocess
import re
import os
import sys

# ---------------------------------------------------------------------------- #
def get_fastq_input(wc):
    samps = SAMPLE_SHEET['samples'][wc.name]['fastq']

    if type(samps) is str:
        samps = [samps]

    infiles = [os.path.join(PATH_FASTQ, i) for i in samps]
    return(infiles)

# ---------------------------------------------------------------------------- #
def get_library_type(wc):
    lib = SAMPLE_SHEET['samples'][wc.name]['library']
    return(lib)

# ---------------------------------------------------------------------------- #
# given a config dictionary sets the default value
def set_default(param, default, dictionary):
	VALUE = ''
	if param in dictionary.keys() and dictionary[param] != None:
		VALUE = dictionary[param]
	else:
		VALUE = default

	return(VALUE)

# ---------------------------------------------------------------------------- #
# given a app name calls the help and parses the parameters
def get_app_params(app_name):
    app_return = subprocess.check_output(SOFTWARE[app_name]['executable'] +' '+ SOFTWARE[app_name]['help'], shell=True)
    app_return = str(app_return)
    vals = list(set(re.findall('(\-{1,2}[a-zA-Z][\w\-]*)\W' , app_return)))
    keys = [re.sub('^-+','',i) for i in vals]
    args = dict(zip(keys, vals))
    return(args)

# ---------------------------------------------------------------------------- #
# WRITE TESTS FOR THIS FUNCTION!!!
def join_params(app, app_params, params_set):
    app_name    = os.path.basename(app)
    params_all  = get_app_params(app_name)
    names       = set(params_all.keys())
    params_diff = set(params_set) - names
    if len(params_diff) > 0:
        print(app_name + 'contains unknown parameters:' + ", ".join(list(params_diff)))
        return(0)
    params_set_keys = set(params_set.keys())
    # check whether some of the arguments are not allowed
    if 'remove' in SOFTWARE[app_name].keys():
        params_set_keys = params_set_keys- set(SOFTWARE[app_name]['remove'])
    params_set_keys = list(params_set_keys)
    params = [params_all[i] +' '+ str(params_set[i]) for i in params_set_keys]
    params = " ".join(params)
    return(params)

# ---------------------------------------------------------------------------- #
# given a sample name extracts whether the sample is broadPeak or narrowPeak
def get_macs2_suffix(name, dictionary):
    sample = dictionary[name]
    suffix = 'narrowPeak'
    if 'macs2' in set(sample.keys()):
        if 'broad' in set(sample['macs2'].keys()):
                suffix = 'broadPeak'
    return(suffix)

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
def RunRscript(input, output, params, script):

    # if isinstance(input, list):
    #     input = dict(zip(input, input))

    print(input.items())
    params_dump = json.dumps(dict(params.items()),sort_keys=True,
                       separators=(",",":"), ensure_ascii=True)
    input_dump  = json.dumps(dict(input.items()),sort_keys=True,
                       separators=(",",":"), ensure_ascii=True)
    output_dump = json.dumps(dict(output.items()),sort_keys=True,
                       separators=(",",":"), ensure_ascii=True)

    cmd = " ".join(['nice -19',str(params.Rscript), os.path.join(SCRIPT_PATH, script),
                    "--basedir", SCRIPT_PATH,
                    "--input",  "{input_dump:q}",
                    "--output", "{output_dump:q}",
                    "--params", "{params_dump:q}"])

    shell(cmd)
