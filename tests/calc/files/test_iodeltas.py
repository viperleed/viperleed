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

from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.files.iodeltas import collect_static_input_files


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
        files_read = []
        def _mock_read_text(path, **kwargs):
            # pylint: disable-next=magic-value-comparison
            assert kwargs['encoding'] == 'utf-8'
            files_read.append(path)
            return 'fake text\n'
        mocker.patch('pathlib.Path.read_text', _mock_read_text)

        return {
            'copy': mocker.patch('shutil.copy2'),
            'is_file': mocker.patch('pathlib.Path.is_file', return_value=True),
            'files_read': files_read,
            'delta_basic': mocker.patch(f'{_MODULE}.generateDeltaBasic'),
            'writeAUXBEAMS': mocker.patch(
                'viperleed.calc.files.beams.writeAUXBEAMS',
                ),
            }

    def test_auxbeams_from_supp(self, rpars, mocks, mocker):
        """Check pulling AUXBEAMS from SUPP."""
        def _root_missing(path):
            return path.parent.name == DEFAULT_SUPP

        mocker.patch('pathlib.Path.is_file', _root_missing)
        collect_static_input_files('mock_slab', rpars)
        mocks['copy'].assert_called_once()
        # writeAUXBEAMS is called because is_file()
        # returns False twice for Path('AUXBEAMS')
        mocks['writeAUXBEAMS'].assert_called_once_with(ivbeams=rpars.ivbeams,
                                                       beamlist=rpars.beamlist)

    def test_auxbeams_from_supp_fails(self, rpars, mocks, mocker, caplog):
        """Check failing to copy AUXBEAMS from SUPP causes log warnings."""
        mocks['copy'].side_effect = OSError
        self.test_auxbeams_from_supp(rpars, mocks, mocker)
        expect_log = 'Failed to copy AUXBEAMS from SUPP'
        assert expect_log in caplog.text

    @pytest.mark.usefixtures('mocks')
    @parametrize(file=('AUXBEAMS', 'PHASESHIFTS'))
    def test_cant_read_needed_file(self, file, rpars, mocker, caplog):
        """Check that failures to read a necessary file are propagated."""
        old_read_text = Path.read_text
        def _fail_to_read_file(path, **kwargs):
            if path.name == file:
                raise OSError
            return old_read_text(path, **kwargs)
        mocker.patch('pathlib.Path.read_text', _fail_to_read_file)
        with pytest.raises(OSError):
            collect_static_input_files('mock_slab', rpars)
        expect_log = f'Could not read {file} for delta'
        assert expect_log in caplog.text

    def test_files_exist(self, rpars, mocks):
        """Check that deltas reads existing static input files."""
        collect_static_input_files('mock_slab', rpars)
        expect_files_read = [Path(f) for f in ('AUXBEAMS', 'PHASESHIFTS')]
        assert expect_files_read == mocks['files_read']
        mocks['writeAUXBEAMS'].assert_not_called()

    @pytest.mark.usefixtures('mocks')
    def test_files_read_end_with_newline(self, rpars, mocker):
        """Check that static input files always end with a newline."""
        no_newline = '''
Missing
newline
at the end'''
        mocker.patch('pathlib.Path.read_text', return_value=no_newline)
        result = collect_static_input_files('mock_slab', rpars)
        expect_args = (no_newline + '\n',) * 2
        assert result[-2:] == expect_args

    def test_success(self, rpars, mocks):
        """Check the successful collection of all input files."""
        result = collect_static_input_files('mock_slab', rpars)
        mocks['delta_basic'].assert_called_once_with('mock_slab', rpars)
        expect = (
            mocks['delta_basic'].return_value,
            'fake text\n',
            'fake text\n',
            )
        assert result == expect

    def test_write_auxbeams_fails(self, rpars, mocks, mocker):
        """Check that failures in writeAUXBEAMS are propagated."""
        mocks['writeAUXBEAMS'].side_effect = Exception('writeAUXBEAMS failed')
        with pytest.raises(Exception, match='writeAUXBEAMS failed'):
            self.test_auxbeams_from_supp(rpars, mocks, mocker)
