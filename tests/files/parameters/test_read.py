"""Tests for modules _read/_reader of viperleed.tleedmlib.files.parameters.

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
from viperleed.tleedmlib.files.parameters._utils import Assignment

from ...helpers import exclude_tags, duplicate_all, CaseTag
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
        info, _ = args
        _, rpars = read_parameters(args)
        assert set(rpars.readParams.keys()) == set(info.parameters.expected)

    def test_stop_commented_out(self, data_path, read_parameters):
        """Check that a 'STOP' is commented out when reading."""
        args = CasesParametersFile().case_stop(data_path)
        fpath, _ = read_parameters(args)
        with fpath.open('r', encoding='utf-8') as parameters_file:
            assert '! STOP' in parameters_file.read()


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


class TestUpdate:
    """Collection of tests for updating from a PARAMETERS file."""

    def test_file_not_found(self):
        """Check exceptions raised when no PARAMETERS file is found."""
        with pytest.raises(FileNotFoundError):
            update(None)

    @fixture(name='read_once')
    def fixture_read_once(self, data_path, read_parameters):
        """Read one PARAMETERS file."""
        args = CasesParametersFile().case_stop(data_path)
        return read_parameters(args)

    def test_nothing_updated(self, read_once):
        """Check that nothing is changed if no PARAMETERS is unchanged."""
        fpath, rpars = read_once
        before, *_ = duplicate_all(rpars.readParams)
        update(rpars, update_from=fpath.parent)
        assert rpars.readParams == before

    def test_updated_irrelevant(self, read_once):
        """Check nothing changes if a non-watched parameter is changed."""
        fpath, rpars = read_once
        before, *_ = duplicate_all(rpars.readParams)
        with fpath.open('a', encoding='utf-8') as parameters_file:
            parameters_file.write('BULK_LIKE_BELOW = 0.25\n')
        update(rpars, update_from=fpath.parent)
        assert rpars.readParams == before

    def test_updated_irrelevant_typo(self, read_once):
        """Check no complaints when a wrong parameter is written."""
        fpath, rpars = read_once
        before, *_ = duplicate_all(rpars.readParams)
        with fpath.open('a', encoding='utf-8') as parameters_file:
            parameters_file.write('SURFACE_ABOVE = 0.25\n')
        update(rpars, update_from=fpath.parent)
        assert rpars.readParams == before

    def test_updated_relevant(self, read_once):
        """Check no stopping occurs if STOP=False is written."""
        fpath, rpars = read_once
        with fpath.open('a', encoding='utf-8') as parameters_file:
            parameters_file.write('SEARCH_CONVERGENCE dgen = 2519 1.5\n')
        update(rpars, update_from=fpath.parent)
        assert rpars.searchConvInit['dgen']['dec'] == 2519
        assert rpars.SEARCH_MAX_DGEN['dec'] == 2519
        assert rpars.SEARCH_MAX_DGEN_SCALING['dec'] == 1.5

    def test_updated_stop_false(self, read_once):
        """Check no stopping occurs if STOP=False is written."""
        fpath, rpars = read_once
        with fpath.open('a', encoding='utf-8') as parameters_file:
            parameters_file.write('STOP = false\n')
        update(rpars, update_from=fpath.parent)
        assert not rpars.STOP

    def test_updated_stop_true(self, read_once):
        """Check no stopping occurs if STOP=False is written."""
        fpath, rpars = read_once
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
            assert next(reader)[1] == Assignment('1-3', 'RUN')
            assert next(reader)[1] == Assignment('50 700 3', 'THEO_ENERGIES')

    def test_raw_reader(self, path_to_params):
        """Check the lines returned by a RawLineParametersReader."""
        expected_lines = (
            ('', '! ####### GLOBAL PARAMETERS #######\n'),
            ('', '\n'),
            ('RUN', 'RUN = 1-3\n'),
            ('THEO_ENERGIES', 'THEO_ENERGIES =   50 700 3\n'),
            )
        with RawLineParametersReader(path_to_params) as reader:
            for expected, _read in zip(expected_lines, reader):
                assert _read == expected
