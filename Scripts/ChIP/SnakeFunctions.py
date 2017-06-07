
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
# given a app name calls the help and parses the parameters
def get_app_params(app, app_name, app_params):
    app_return = subprocess.check_output(app +' '+ app_params[app_name]['help'], shell=True)
    app_return = str(app_return)
    vals = list(set(re.findall('(\-{1,2}[a-zA-Z][\w\-]*)\W' , app_return)))
    keys = [re.sub('^-+','',i) for i in vals]
    args = dict(zip(keys, vals))
    return(args)

# ---------------------------------------------------------------------------- #
def join_params(app, app_params, params_set):
    app_name = os.path.basename(app)
    params_all = get_app_params(app, app_name, app_params)
    names  = set(params_all.keys())
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

