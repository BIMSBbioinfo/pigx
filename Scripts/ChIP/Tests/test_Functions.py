import os
import sys
sys.path.append('..')
import yaml
import pytest
import collections
from SnakeFunctions import *

# ---------------------------------------------------------------------------- #
# Test all

config = yaml.load(open('../config.yaml'))
software = yaml.load(open('../software.yaml'))

get_app_params

@pytest.fixture
def load_config_all():
    return yaml.load(open('../config.yaml'))

# test config: properly formatted config
def test_config(load_config_all):
    assert check_config(load_config_all) == 0

# loads config corruptors and tests whether the config - check works
config_replace = yaml.load(open('config_replace.yaml'))
@pytest.fixture(params = config_replace.keys())
def load_config_replace(request):
    config = yaml.load(open('../config.yaml'))
    return update(config, config_replace[request.param])

def test_config_fail(load_config_replace):
    assert check_config(load_config_replace) == 1