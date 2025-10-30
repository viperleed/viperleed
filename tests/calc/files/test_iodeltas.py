"""Tests for module iodeltas of viperleed.calc.files."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-07-08'
__license__ = 'GPLv3+'

from pathlib import Path

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.files.iodeltas import _fetch_auxbeams
from viperleed.calc.files.iodeltas import _read_file_with_newline
from viperleed.calc.files.iodeltas import collect_static_input_files
from viperleed.calc.files.iodeltas import write_delta_input_file
from viperleed.calc.lib.context import execute_in_dir


_MODULE = 'viperleed.calc.files.iodeltas'


@fixture(name='rpars')
def fixture_rpars(mocker):
    """Return a fake Rparams."""
    return mocker.MagicMock()


class TestCollectStaticInputFiles:
    """Tests for function collect_static_input_files."""

    @fixture(name='mocks')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details of the deltas function with mocks."""
        return {
            'delta_basic': mocker.patch(f'{_MODULE}.generateDeltaBasic'),
            'fetch_aux': mocker.patch(f'{_MODULE}._fetch_auxbeams'),
            'read': mocker.patch(f'{_MODULE}._read_file_with_newline'),
            }

    def test_implementation(self, mocks, mocker):
        """Check delegation to helper functions."""
        _, mock_rpars = mock_args = mocker.MagicMock(), mocker.MagicMock()
        result = collect_static_input_files(*mock_args)
        mocks['delta_basic'].assert_called_once_with(*mock_args)
        mocks['fetch_aux'].assert_called_once_with(mock_rpars)
        read_calls = [
            mocker.call(mocks['fetch_aux'].return_value),
            mocker.call('PHASESHIFTS'),
            ]
        assert mocks['read'].mock_calls == read_calls
        expect = (
            mocks['delta_basic'].return_value,
            mocks['read'].return_value,
            mocks['read'].return_value,
            )
        assert result == expect

    @parametrize(exc=(BaseException, Exception))
    def test_propagate_exceptions(self, exc, mocks):
        """Check that no exceptions are caught."""
        for mock in mocks.values():
            mock.side_effect = exc
            with pytest.raises(exc):
                collect_static_input_files('mock_slab', 'mock_rpars')


class TestFetchAuxbeams:
    """Tests for the _fetch_auxbeams helper."""

    @fixture(name='mocks')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details of the deltas function with mocks."""
        return {
            'is_file': mocker.patch('pathlib.Path.is_file', return_value=True),
            'copy': mocker.patch('shutil.copy2'),
            'writeAUXBEAMS': mocker.patch(
                'viperleed.calc.files.beams.writeAUXBEAMS',
                ),
            }

    @fixture(name='call')
    # pylint: disable-next=unused-argument  # Cannot .mark fixtures
    def fixture_call(self, rpars, mocks):
        """Call _fetch_auxbeams and check the return value."""
        def _call():
            result = _fetch_auxbeams(rpars)
            assert result == Path('AUXBEAMS')
        return _call

    def test_exists(self, call, mocks):
        """Check that nothing happens if AUXBEAMS is present in CWD."""
        call()
        not_called = (
            'copy',
            'writeAUXBEAMS',
            )
        for mock_name in not_called:
            mocks[mock_name].assert_not_called()

    def test_fetch_from_supp(self, call, rpars, mocks):
        """Check pulling AUXBEAMS from SUPP."""
        mocks['is_file'].return_value = False
        call()
        mocks['copy'].assert_called_once()
        # writeAUXBEAMS is called because is_file()
        # returns False twice for Path('AUXBEAMS')
        mocks['writeAUXBEAMS'].assert_called_once_with(ivbeams=rpars.ivbeams,
                                                       beamlist=rpars.beamlist)

    def test_fetch_from_supp_fails(self, call, rpars, mocks, caplog):
        """Check failing to copy AUXBEAMS from SUPP causes log warnings."""
        mocks['copy'].side_effect = OSError
        self.test_fetch_from_supp(call, rpars, mocks)
        expect_log = 'Failed to copy AUXBEAMS from SUPP'
        assert expect_log in caplog.text

    def test_not_found(self, call, rpars, mocks):
        """Check behavior when AUXBEAMS is not present in CWD or SUPP."""
        mocks['copy'].side_effect = FileNotFoundError
        self.test_fetch_from_supp(call, rpars, mocks)

    def test_write_anew_fails(self, call, rpars, mocks):
        """Check that failures in writeAUXBEAMS are propagated."""
        mocks['writeAUXBEAMS'].side_effect = Exception('writeAUXBEAMS failed')
        with pytest.raises(Exception, match='writeAUXBEAMS failed'):
            self.test_fetch_from_supp(call, rpars, mocks)


class TestReadFileWithNewline:
    """Tests for the _read_file_with_newline helper."""

    def test_cant_read_file(self, mocker, caplog):
        """Check that failures to read a necessary file are propagated."""
        mocker.patch('pathlib.Path.read_text', side_effect=OSError)
        file = 'read_fails'
        with pytest.raises(OSError):
            _read_file_with_newline(file)
        expect_log = f'Could not read {file} for delta'
        assert expect_log in caplog.text

    def test_file_has_newline(self, mocker):
        """Check that deltas reads existing static input files."""
        with_newline = 'contents\n'
        mocker.patch('pathlib.Path.read_text', return_value=with_newline)
        result = _read_file_with_newline('file')
        assert result == with_newline

    def test_file_no_newline(self, mocker):
        """Check that static input files always end with a newline."""
        no_newline = '''
Missing
newline
at the end'''
        mocker.patch('pathlib.Path.read_text', return_value=no_newline)
        result = _read_file_with_newline('file')
        expect = no_newline + '\n'
        assert result == expect


class TestWriteDeltaInputFile:
    """Tests for the write_delta_input_file function."""

    @fixture(name='mock_tasks')
    def fixture_mock_tasks(self, mocker):
        """Return fake compile and run tasks."""
        compile_tasks = (
            mocker.MagicMock(param='param 1. Has one runtaks'),
            mocker.MagicMock(param='param 2. Has multiple runtasks'),
            mocker.MagicMock(param='param 3. Has no runtasks'),
            )
        run_tasks = (
            mocker.MagicMock(comptask=compile_tasks[0],
                             deltaname='DEL_runtask_1',
                             din_short='comp 1, short 1'),
            mocker.MagicMock(comptask=compile_tasks[1],
                             deltaname='DEL_runtask_2',
                             din_short='comp 2, short 1'),
            mocker.MagicMock(comptask=compile_tasks[1],
                             deltaname='DEL_runtask_3',
                             din_short='comp 2, short 2'),
            )
        return compile_tasks, run_tasks

    def test_success(self, tmp_path, mock_tasks):
        """Check the successful writing of a delta-input file."""
        with execute_in_dir(tmp_path):
            write_delta_input_file(*mock_tasks)
        contents = (tmp_path/'delta-input').read_text(encoding='utf-8')
        expect_contents = '''\
# ABOUT THIS FILE:
# Input for the delta-calculations is collected here. The blocks of data are
# new 'PARAM' files, which are used to recompile the fortran code, and input
# for generation of specific DELTA files. Lines starting with '#' are comments
# on the function of the next block of data.
# In the DELTA file blocks, [AUXBEAMS] and [PHASESHIFTS] denote where the
# entire contents of the AUXBEAMS and PHASESHIFTS files should be inserted.

#### NEW 'PARAM' FILE: ####

param 1. Has one runtaks

#### INPUT for new DELTA file DEL_runtask_1: ####

comp 1, short 1

#### NEW 'PARAM' FILE: ####

param 2. Has multiple runtasks

#### INPUT for new DELTA file DEL_runtask_2: ####

comp 2, short 1

#### INPUT for new DELTA file DEL_runtask_3: ####

comp 2, short 2

#### NEW 'PARAM' FILE: ####

param 3. Has no runtasks
'''
        assert contents == expect_contents

    def test_write_fails(self, mock_tasks, mocker, caplog):
        """Check that failure to write the input file is tolerated."""
        input_file_name = 'delta-input'
        mocker.patch('pathlib.Path.write_text', side_effect=OSError)
        write_delta_input_file(*mock_tasks)
        expect_log = f'Failed to write file {input_file_name!r}'
        assert expect_log in caplog.text
