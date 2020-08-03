"""
Test functionality of leedsim.classes.py module without loading the whole TLEEDMAP

Created: 2020-01-12
Author: Michele Riva

This module is just a wrapper that imports all test submodules and can be
imported by tleedmap.test.test_all to run tests of the class submodule
"""

import sys
import os 

tests_path = os.path.realpath(os.path.dirname(__file__))
base_path = os.path.realpath(os.path.join(tests_path, '..'))
for path in [tests_path, base_path]:
    if path not in sys.path:
        sys.path.append(path)

import pytest
import numpy as np

import classes
from test_realspace import *
