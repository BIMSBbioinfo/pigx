from snakemake import shell
import sys
import subprocess
import re

# ----------------------------------------------------------------------------- #
# java_head_difference is the reduction in heap allocation given to the java executible
# compared to SGE submission: if SGE submission requests 16G, java will be
# started with (16 - java_heap_difference)G
def java_tool(java, threads, mem, tempdir, tool_path, tool_name, java_heap_difference=5):

    # removes 1 g from java memory heap
    mem_size   = int(float(mem[:-1]))
    mem_suffix = mem[-1]
    if(mem_size <= java_heap_difference):
        mem_size = java_heap_difference - 1
    else:
        mem_size = mem_size - java_heap_difference

    mem_reduced = str(mem_size) + mem_suffix

    tool = ' '.join([
        java,
        '-XX:ParallelGCThreads=' + str(threads),
        '-Xmx'                   + str(mem_reduced),
        '-Xss'                   + str('16M'),
        '-Djava.io.tmpdir='      + str(tempdir),
        '-jar',
        str(tool_path),
        str(tool_name)
    ])
    return tool

# ----------------------------------------------------------------------------- # prints the command to STDERR and executes
def print_shell(command):
    print(command, file=sys.stderr)
    shell(command)
    

# ----------------------------------------------------------------------------- 
# extracts the adapter start and adapter length from the adapter hash
# sample_name is the sample name defined in sample sheet
# type : umi_barcode / cell_barcode
def adapter_params(sample_name, type):

    if not type in set(['cell_barcode','umi_barcode']):
        sys.exit('invalid barcode type')

    method = SAMPLE_SHEET.fetch_field(sample_name,'method')[0]
    adapter_params = ADAPTER_PARAMETERS[method][type]
    barcode_hash = {
        'start'  : adapter_params['base_min'],
        'length' : adapter_params['base_max'] - adapter_params['base_min'] + 1
    }
    return(barcode_hash)
# ----------------------------------------------------------------------------- # calculates the barcode length from the sample sheet
def get_adapter_size(name):
    cb_adapter   = adapter_params(name, 'cell_barcode')
    umi_adapter  = adapter_params(name, 'umi_barcode')
    adapter_size = cb_adapter['length'] + umi_adapter['length']
    return(adapter_size)


# ---------------------------------------------------------------------------- #
# given a app name calls the help and parses the parameters
def get_app_params(app_name):
    app_return = subprocess.check_output(SOFTWARE[app_name]['executable'] +' '+ SOFTWARE[app_name]['help'], shell=True)
    app_return = str(app_return)
    vals = list(set(re.findall('^(\-{1,2}[a-zA-Z][\w\-]*)\W' , app_return)))
    keys = [re.sub('^-+','',i) for i in vals]
    args = dict(zip(keys, vals))
    return(args)

# ---------------------------------------------------------------------------- #
# star requires a separate function for parsing parameters
# help doesn't include --
def get_star_params():
    app_return = subprocess.check_output(SOFTWARE['STAR']['executable'] +' '+ SOFTWARE['STAR']['help'], shell=True)
    app_return = str(app_return)
    keys = list(set(re.findall('\\\\n([a-zA-Z]+)\W' , app_return)))
    vals = ['--' + i for i in keys]
    args = dict(zip(keys, vals))
    return(args)


# ---------------------------------------------------------------------------- #
# WRITE TESTS FOR THIS FUNCTION!!!
def join_params(app, app_params, params_set):
    app_name    = os.path.basename(app)
    if app == 'STAR':
        params_all = get_star_params()
    else:
        params_all  = get_app_params(app_name)
    
    # checks whether any parameters are defined
    if params_set == None:
        params = ""
    else:
        names       = set(params_all.keys())
        params_diff = set(params_set) - names
        if len(params_diff) > 0:
            message = app_name + 'contains unknown parameters: ' + ", ".join(list(params_diff))
            sys.exit(message)
        
        params_set_keys = set(params_set.keys())
        # check whether some of the arguments are not allowed
        if 'remove' in SOFTWARE[app_name].keys():
            params_set_keys = params_set_keys- set(SOFTWARE[app_name]['remove'])
        params_set_keys = list(params_set_keys)
        params = [params_all[i] +' '+ str(params_set[i]) for i in params_set_keys]
        params = " ".join(params)
        
    return(params)
    