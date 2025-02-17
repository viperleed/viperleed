"""Tests for module cleanup of viperleed.calc.section."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-03'
__license__ = 'GPLv3+'

from pathlib import Path

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.constants import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES
from viperleed.calc.sections.cleanup import OPTIONAL_INPUT_FILES
from viperleed.calc.sections.cleanup import preserve_original_inputs

from .conftest import _MODULE


class TestPreserveOriginalInputs:
    """Tests for the preserve_original_inputs function."""

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Return a fake shutil.copy2 and a logger."""
        def _mock(**kwargs):
            mock_copy = mocker.patch('shutil.copy2', **kwargs)
            mock_logger = mocker.patch(f'{_MODULE}.logger')
            return mock_copy, mock_logger
        return _mock

    @fixture(name='make_inputs')
    def fixture_make_inputs(self, tmp_path):
        """Create the given input files in a temporary directory."""
        def _make(files):
            for file in files:
                (tmp_path / file).touch()
            return tmp_path
        return _make

    @fixture(name='run_preserve')
    def fixture_run_preserve(self, rpars, make_inputs, mock_implementation):
        """Execute preserve_original_inputs."""
        def _run(input_files, **kwargs):
            tmp_path = make_inputs(input_files)
            mocked = mock_implementation(**kwargs)
            with execute_in_dir(tmp_path):
                preserve_original_inputs(rpars)
            assert (tmp_path/ORIGINAL_INPUTS_DIR_NAME).is_dir()
            return (*mocked, rpars)
        return _run

    _skip_expbeams = {  # skip: keep
        'EXPBEAMS': 'EXPBEAMS.csv',
        'EXPBEAMS.csv': 'EXPBEAMS',
        }

    @parametrize('skip,expect', _skip_expbeams.items())
    def test_expbeams_selection(self, skip, expect, run_preserve, tmp_path):
        """Test that the first existing EXPBEAMS file is selected."""
        input_files = ALL_INPUT_FILES - {skip}
        nr_files_copied = len(input_files)
        copy, logger, *_ = run_preserve(input_files)

        assert copy.call_count == nr_files_copied
        logger.warning.assert_not_called()
        assert (tmp_path/expect).is_file()

        # Check that we copied the expected one to original_inputs
        # Notice that we can't just use assert_any_call, as the
        # arguments may be OS-dependent with Paths (e.g., on WSL)
        calls = {
            (Path(src).name, Path(dst).name)
            for (src, dst), *_ in copy.call_args_list
            }
        assert (expect, ORIGINAL_INPUTS_DIR_NAME) in calls

    def test_copy_fails(self, run_preserve):
        """Check expected complaints when failing to copy files."""
        _, logger, rpars = run_preserve(ALL_INPUT_FILES, side_effect=OSError)
        assert logger.warning.call_count == len(ALL_INPUT_FILES) - 1
        assert rpars.halt == 1

    def test_mkdir_fails(self, rpars, tmp_path, mocker):
        """Check complaints when failing to create original_inputs."""
        mocker.patch('pathlib.Path.mkdir', side_effect=OSError)
        with pytest.raises(OSError):
            preserve_original_inputs(rpars)
        assert not (tmp_path / ORIGINAL_INPUTS_DIR_NAME).exists()

    def test_required_missing(self, run_preserve):
        """Check that a missing mandatory input file causes warnings."""
        missing_required = ('PARAMETERS',)
        copy, logger, rpars = run_preserve(missing_required)

        logger.warning.assert_called()
        assert copy.call_count == 1
        assert rpars.halt == 1

    cwd_inputs = {
        'simple': ALL_INPUT_FILES,
        'missing optional': ALL_INPUT_FILES - set(OPTIONAL_INPUT_FILES),
        }

    @parametrize(input_files=cwd_inputs.values(), ids=cwd_inputs)
    def test_success(self, input_files, run_preserve):
        """Check execution of a successful storage."""
        # There's one less because we have two acceptable EXPBEAMS
        nr_files_copied = len(input_files) - 1
        copy, logger, *_ = run_preserve(input_files)

        assert copy.call_count == nr_files_copied
        logger.warning.assert_not_called()
