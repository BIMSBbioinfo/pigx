# ---------------------------------------------------------------------------- #
import subprocess
import re
import os
import sys


# ---------------------------------------------------------------------------- #
# given a fastq file strips of either of the following extensions:
# "fastq.gz","fastq","fq.gz","fq" 
# this is required for parsing trim galore output
def replace_fastq_ext(fqfile,replacement):
    p = re.compile('.f(ast)*q(\.gz)*')
    repl = p.sub(replacement,fqfile)
    return(repl)

# ---------------------------------------------------------------------------- #
# given a sample name returns fastq location(s)
def get_fastq_input(name):
    samps = lookup('SampleName', name, ['Read', 'Read2'])
    infiles = [os.path.join(PATH_FASTQ, i) for i in samps if i]
    return(infiles)

# ---------------------------------------------------------------------------- #
# given a sample name returns fastq location(s)
def get_trimmed_input(name):
    fqfiles = [file for file in lookup('SampleName',name,['Read','Read2']) if file]
    if len(fqfiles) == 2:
        trimmed_files = [os.path.join(PATH_TRIMMED,name, "{}_{}.fastq.gz".format(name,read)) for read in ["R1","R2"]]
    else:    
        trimmed_files = [os.path.join(PATH_TRIMMED,name, "{}_R.fastq.gz".format(name))]
    return(trimmed_files)

# ---------------------------------------------------------------------------- #
# given a sample name returns library type, depending on number of fastq files
def get_library_type(name):
    files = [file for file in lookup('SampleName', name, ['Read', 'Read2']) if file]
    if len(files) == 2:
        lib = "paired"
    else:
        lib = "single"
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

    if not sample == None:
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

    cmd = " ".join(['nice -19',str(params.Rscript),
                    SOFTWARE['Rscript']['args'],
                    os.path.join(SCRIPT_PATH, script),
                    "--basedir", SCRIPT_PATH,
                    "--input",  "{input_dump:q}",
                    "--output", "{output_dump:q}",
                    "--params", "{params_dump:q}"])

    shell(cmd)
