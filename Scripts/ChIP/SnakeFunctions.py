
import subprocess
import re
import os
import sys

# ---------------------------------------------------------------------------- #
# given a config dictionary sets the default value
def set_default(param, default, config):
	VALUE = ''
	if param in config.keys() and config[param] != None:
		VALUE = config[param]
	else:
		VALUE = default

	return(VALUE)

# ---------------------------------------------------------------------------- #
# checks the config file for validity
# config check
def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False


def check_config(config):

	message = ''

	# checks for proper top level config categories
	params = ['genome','genome_fasta','index','fastq','params','samples','peak_calling','idr','software']
	params_diff = set(config.keys()) - set(params)
	if len(params_diff) > 0:
		message = message + "config file contains unknown parameters"

	# checks for index or genome specification
	if (config['genome'] == None) and (config['index'] == None):
		message = message + "neither genome nor index are specified\n"

	# checks for correspondence between peak calling and samples
	samples = list(config['samples'].keys())
	peaks = [config['peak_calling'][i].values() for i in list(config['peak_calling'].keys())]
	samples_diff = (set(peaks[0]) - set(samples))
	if len(samples_diff) > 0:
		message = message + "some peak calling samples are not specified\n"

	# checks whether extend is a number
	# if not (is.number(config['params']['extend'])):
	#     message = message + "extend must be a number\n"

	if len(message) > 0:
		print(message)
		return(1);

	return(0)



# ---------------------------------------------------------------------------- #
# given a app name calls the help and parses the parameters
def get_app_params(app, app_name, app_params):
    import subprocess
    import re
    app_return = subprocess.check_output(app +' '+ app_params[app_name]['help'], shell=True)
    app_return = str(app_return)
    vals = list(set(re.findall('(\-{1,2}[a-zA-Z][\w\-]*)\W' , app_return)))
    keys = [re.sub('^-+','',i) for i in vals]
    args = dict(zip(keys, vals))
    return(args)

# ---------------------------------------------------------------------------- #
def join_params(app, app_params, config):
    import subprocess
    import os
    import re
    app_name = os.path.basename(app)
    params_all = get_app_params(app, app_name, app_params)
    names  = set(params_all.keys())
    params_set = config[app_name]
    params_diff = set(params_set) - names
    if len(params_diff) > 0:
        print(app_name + 'contains unknown parameters:' + ", ".join(list(params_diff)))
        return(0)
    params_set_keys = set(params_set.keys())
    # check whether some of the arguments are disallower
    if 'remove' in app_params[app_name].keys():
        params_set_keys = params_set_keys- set(app_params[app_name]['remove'])
    params_set_keys = list(params_set_keys)
    params = [params_all[i] +' '+ str(params_set[i]) for i in params_set_keys]
    params = " ".join(params)
    return(params)




