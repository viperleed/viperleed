"""Tests for module _read of viperleed.tleedmlib.files.parameters.

Created on 2023-10-17

@author: Michele Riva (@michele-riva)
"""

import logging
import shutil

import pytest
from pytest_cases import fixture, parametrize_with_cases

from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tleedmlib.files.parameters.errors import (
    ParameterNotRecognizedError, ParameterHasNoValueError
    )
from viperleed.tleedmlib.files.parameters._read import read, update
from viperleed.tleedmlib.files.parameters._reader import (
    ParametersReader, RawLineParametersReader
    )

from ...helpers import exclude_tags, CaseTag
from .case_parameters import CasesParametersFile


@fixture(scope='class', name='read_parameters')
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


_NOT_RAISING = {'cases': CasesParametersFile,
                'filter': exclude_tags(CaseTag.RAISES)}


class TestReadSuccessful:
    """Collection of tests for successfully reading from a PARAMETERS file."""

    @parametrize_with_cases('args', **_NOT_RAISING)
    def test_reads_rpars(self, args, read_parameters):
        """Check that reading of file succeeds."""
        _, rpars = read_parameters(args)
        assert isinstance(rpars, Rparams)

    @parametrize_with_cases('args', **_NOT_RAISING)
    def test_read_not_empty(self, args, read_parameters):
        """Check that reading of file succeeds."""
        _, rpars = read_parameters(args)
        info, _ = args
        if 'PARAMETERS_empty' in info.parameters.param_path:
            pytest.skip(reason='File is empty')
        assert rpars.readParams

    @parametrize_with_cases('args', **_NOT_RAISING)
    def test_read_expected(self, args, read_parameters):
        """Check that reading of file succeeds."""
        info, fpath = args
        _, rpars = read_parameters(args)
        assert set(rpars.readParams.keys()) == set(info.parameters.expected)


class TestReadFailing:
    """Collection of tests for reading from an improper PARAMETERS file."""

    def test_file_not_found(self):
        """Check complaints if no PARAMETERS file is found."""
        with pytest.raises(FileNotFoundError):
            read('__not_a_PARAMETERS_file')

    def test_missing_equals_logs(self, data_path, read_parameters,
                                 caplog, re_match):
        """Check that a PARAMETERS files without an '=' causes complaints."""
        args = CasesParametersFile().case_misses_equals(data_path)
        with caplog.at_level(logging.WARNING):
            _, rpars = read_parameters(args)
        # Notice the [\s\S] to match also newlines
        assert re_match(r'.*[\s\S]*line.*without.*=.*sign.*SKIPPED',
                        caplog.text)
        assert rpars.halt >= 2

    def test_missing_value_raises(self, data_path, read_parameters):
        """Check complaints when reading an unknown parameter from file."""
        args = CasesParametersFile().case_no_value_parameter(data_path)
        with pytest.raises(ParameterHasNoValueError):
            read_parameters(args)

    def test_unknown_parameter_raises(self, data_path, read_parameters):
        """Check complaints when reading an unknown parameter from file."""
        args = CasesParametersFile().case_unknown_parameter(data_path)
        with pytest.raises(ParameterNotRecognizedError):
            read_parameters(args)
