import os
import sys
sys.path.append('..')
import yaml
import pytest
import collections
from Check_Config import *

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
# Test all

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


