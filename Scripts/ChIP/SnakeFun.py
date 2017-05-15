# ---------------------------------------------------------------------------- #
# given a config dictionary sets the default value
def set_default(param, value, config):
	if param in config.keys() and config[param] != None:
		VALUE = config[param]
	else:
		VALUE = value

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
	params = ['genome','index','params','samples','peak_calling','idr','fastq','software']
	params_diff = set(config.keys()) - set(params)
	if len(params_diff) > 0:
		message = message + "config file contains dissalowed categories\n"

	# checks for correspondence between peak calling and samples
	samples = list(config['samples'].keys())
	peaks = [config['peak_calling'][i].values() for i in list(config['peak_calling'].keys())]
	samples_diff = (set(peaks[0]) - set(samples))
	if len(samples_diff) > 0:
		message = message + "some peak calling samples are not specified\n"

	# checks for index or genome specification
	if len(config['genome']) == 0 and len(config['index']):
		message = message + "neither genome nor index are specified\n"

	# checks whether extend is a number
	# if not (is.number(config['params']['extend'])):
	#     message = message + "extend must be a number\n"

	if len(message) > 0:
		print(message)
		return(1);

	return(0)
