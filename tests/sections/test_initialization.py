"""Tests for section initialization.

Created on 2023-07-19

@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)
"""

import pytest


class TestSetup:
    """Basic tests for pre-running preparation of the work environment."""

    def test_work_path_exists(self, init_files):
        """Check that work_path was created properly."""
        assert init_files.work_path.is_dir()

    def test_files_copied_correctly(self, init_files):
        """Check if all files were copied correctly."""
        input_files = set(f.name for f in init_files.test_path.glob('*'))
        source_copied = all(f in input_files
                            for f in init_files.required_files)
        input_files = set(f.name for f in init_files.work_path.glob('*'))
        work_copied = all(f in input_files for f in init_files.required_files)
        assert source_copied and work_copied


class TestInitialization:                                                       # TODO: find a way to inject the exit-code test (e.g., a class decorator?)
    """Collection of tests for a successful INITIALIZATION run."""

    def test_successful_run(self, init_files):
        """Check that initialization exits without errors."""
        assert not init_files.failed

    _expected_files = 'IVBEAMS', 'BEAMLIST', 'VIBROCC', 'PARAMETERS'

    @pytest.mark.parametrize('expected_file', _expected_files)
    def test_init_files_present(self, init_files, expected_file):
        """Ensure the expected files are present after initialization."""
        assert init_files.expected_file_exists(expected_file)

    def test_parameters_was_updated(self, init_files):
        """Check that PARAMETERS file was updated."""
        parameters = init_files.work_path / 'PARAMETERS'
        with parameters.open('r', encoding='utf-8') as param_file:
            param_content = param_file.read()
        assert 'line commented out automatically' in param_content
