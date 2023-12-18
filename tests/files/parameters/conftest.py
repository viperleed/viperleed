"""Test configuration for viperleed.tleedmlib.files.parameters.

Created on 2023-10-18

@author: Michele Riva (@michele-riva)
"""

import shutil

from pytest_cases import fixture

from viperleed.tleedmlib.files.parameters._read import read

from .case_parameters import CasesParametersFile


@fixture(scope='session', name='read_parameters')
def factory_read_parameters(tmp_path_factory):
    """Return an Rparams read from file and the path to it."""
    def _move_to_temp_and_read(args):
        _, fpath = args
        tmp_dir = tmp_path_factory.mktemp(basename='parameters', numbered=True)
        new_fpath = tmp_dir / 'PARAMETERS'
        shutil.copy2(fpath, new_fpath)
        rpars = read(new_fpath)
        return new_fpath, rpars
    return _move_to_temp_and_read


@fixture
def read_one_param_file(data_path, read_parameters):
    """Read one example PARAMETERS file."""
    args = CasesParametersFile().case_stop(data_path)
    return read_parameters(args)
