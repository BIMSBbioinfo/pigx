# ---------------------------------------------------------------------------- #
"""
Running pytests:
python -m pytest test_Config.py

"""

import os
import sys
sys.path.append('..')
import yaml
import pytest
import collections
from Check_Config import *

# ---------------------------------------------------------------------------- #
config_path = '../'
config_replace = yaml.load(open('config_replace.yaml'))
config_replace_sample_sheet = config_replace['sample_sheet']
config_replace_settings     = config_replace['settings']

# ---------------------------------------------------------------------------- #
# Recursive function that updates values in hash d with given values from hash u
# u has to correspond to a complete path to the values in d
def update(d, u):
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            r = update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d

# ---------------------------------------------------------------------------- #
# loader functions
def load_config(filename, config_path):
    yaml_file = os.path.join(config_path, (filename + '.yaml'))
    return yaml.load(open(yaml_file))

@pytest.fixture
def load_sample_sheet():
    return load_config('sample_sheet', config_path)

# test settings
@pytest.fixture
def load_settings():
   return load_config('settings', config_path)

# ---------------------------------------------------------------------------- #
@pytest.fixture(params = config_replace_settings.keys())
def load_settings_replace(request):
    settings_corrupt = update(load_settings(), config_replace_settings[request.param])
    return settings_corrupt


def test_settings_fail(load_settings_replace):
    sample_sheet = load_sample_sheet()
    with pytest.raises(SystemExit):
        check_proper_settings_configuration(load_settings_replace, sample_sheet)


# ---------------------------------------------------------------------------- #
# test sample_sheet
# test default values
def test_sample_sheet_default():
    settings     = load_settings()
    sample_sheet = load_sample_sheet()
    assert check_proper_settings_configuration(settings,sample_sheet) == 0

# constructs a config corruptor function
@pytest.fixture(params = config_replace_sample_sheet.keys())
def load_sample_sheet_replace(request):
    sample_sheet_corrupt = update(load_sample_sheet(), config_replace_sample_sheet[request.param])
    return sample_sheet_corrupt

def test_sample_sheet_fail(load_sample_sheet_replace):
    settings = load_settings()
    with pytest.raises(SystemExit):
        check_proper_settings_configuration(settings, load_sample_sheet_replace)


# loads config corruptors and tests whether the config - check works
