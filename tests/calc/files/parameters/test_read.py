"""Tests for modules _read/_reader of viperleed.calc.files.parameters."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-17'
__license__ = 'GPLv3+'

import logging

import pytest
from pytest_cases import fixture, parametrize_with_cases

from viperleed.calc.classes.rparams import Rparams
from viperleed.calc.files.parameters.errors import (
    ParameterNotRecognizedError, ParameterHasNoValueError
    )
from viperleed.calc.files.parameters._read import read, update
from viperleed.calc.files.parameters._reader import (
    ParametersReader, RawLineParametersReader
    )
from viperleed.calc.files.parameters._utils import Assignment

from ....helpers import exclude_tags, duplicate_all
from ...tags import CaseTag
from .case_parameters import CasesParametersFile


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
        info, _ = args
        _, rpars = read_parameters(args)
        assert set(rpars.readParams.keys()) == set(info.parameters.expected)

    def test_stop_commented_out(self, data_path, read_parameters, re_match):
        """Check that a 'STOP' is commented out when reading."""
        args = CasesParametersFile().case_stop(data_path)
        fpath, _ = read_parameters(args)
        with fpath.open('r', encoding='utf-8') as parameters_file:
            assert re_match(r'.*[\s\S]*[!%#]\s*STOP\s*! Disabled',
                            parameters_file.read())


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

    def test_typo(self, data_path, read_parameters):
        """Check suggestions when a typo is present."""
        args = CasesParametersFile().case_typo(data_path)
        with pytest.raises(ParameterNotRecognizedError) as exc:
            read_parameters(args)
        assert exc.match(".*ENERGEIS:\n.*Did you mean 'THEO_ENERGIES'?")


class TestUpdate:
    """Collection of tests for updating from a PARAMETERS file."""

    def test_file_not_found(self):
        """Check exceptions raised when no PARAMETERS file is found."""
        with pytest.raises(FileNotFoundError):
            update(None)

    def test_nothing_updated(self, read_one_param_file):
        """Check that nothing is changed if no PARAMETERS is unchanged."""
        fpath, rpars = read_one_param_file
        before, *_ = duplicate_all(rpars.readParams)
        update(rpars, update_from=fpath.parent)
        assert rpars.readParams == before

    def test_updated_irrelevant(self, read_one_param_file):
        """Check nothing changes if a non-watched parameter is changed."""
        fpath, rpars = read_one_param_file
        before, *_ = duplicate_all(rpars.readParams)
        with fpath.open('a', encoding='utf-8') as parameters_file:
            parameters_file.write('BULK_LIKE_BELOW = 0.25\n')
        update(rpars, update_from=fpath.parent)
        assert rpars.readParams == before

    def test_updated_irrelevant_typo(self, read_one_param_file):
        """Check no complaints when a wrong parameter is written."""
        fpath, rpars = read_one_param_file
        before, *_ = duplicate_all(rpars.readParams)
        with fpath.open('a', encoding='utf-8') as parameters_file:
            parameters_file.write('SURFACE_ABOVE = 0.25\n')
        update(rpars, update_from=fpath.parent)
        assert rpars.readParams == before

    def test_updated_relevant(self, read_one_param_file):
        """Check no stopping occurs if STOP=False is written."""
        fpath, rpars = read_one_param_file
        with fpath.open('a', encoding='utf-8') as parameters_file:
            parameters_file.write('SEARCH_CONVERGENCE dgen = 2519 1.5\n')
        update(rpars, update_from=fpath.parent)
        assert rpars.searchConvInit['dgen']['dec'] == 2519
        assert rpars.SEARCH_MAX_DGEN['dec'] == 2519
        assert rpars.SEARCH_MAX_DGEN_SCALING['dec'] == 1.5

    def test_updated_stop_false(self, read_one_param_file):
        """Check no stopping occurs if STOP=False is written."""
        fpath, rpars = read_one_param_file
        with fpath.open('a', encoding='utf-8') as parameters_file:
            parameters_file.write('STOP = false\n')
        update(rpars, update_from=fpath.parent)
        assert not rpars.STOP

    def test_updated_stop_true(self, read_one_param_file):
        """Check no stopping occurs if STOP=False is written."""
        fpath, rpars = read_one_param_file
        with fpath.open('a', encoding='utf-8') as parameters_file:
            parameters_file.write('STOP\n')
        update(rpars, update_from=fpath.parent)
        assert rpars.STOP


class TestReader:
    """Collection of simple tests for ParametersReader classes."""

    @fixture(name='path_to_params', scope='session')
    def fixture_path_to_params(self, data_path):
        """Return the path to a PARAMETERS file."""
        _, fpath = CasesParametersFile().case_stop(data_path)
        return fpath

    def test_reader(self, path_to_params):
        """Check the lines returned by a ParametersReader."""
        with ParametersReader(path_to_params) as reader:
            # pylint: disable=protected-access
            assert next(reader)[1] == Assignment('1-3', 'RUN')
            assert reader._current_line == 3
            assert next(reader)[1] == Assignment('50 700 3', 'THEO_ENERGIES')
            assert reader._current_line == 4

    def test_raw_reader(self, path_to_params):
        """Check the lines returned by a RawLineParametersReader."""
        expected_lines = (
            ('', '! ####### GLOBAL PARAMETERS #######\n'),
            ('', '\n'),
            ('RUN', 'RUN = 1-3\n'),
            ('THEO_ENERGIES', 'THEO_ENERGIES =   50 700 3\n'),
            )
        with RawLineParametersReader(path_to_params) as reader:
            # pylint: disable=protected-access
            for i, (expected, _read) in enumerate(zip(expected_lines, reader)):
                assert _read == expected
                assert reader._current_line == i + 1
