"""Tests for module context of viperleed.calc.lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-04'
__license__ = 'GPLv3+'

from pathlib import Path

import pytest
from pytest_cases import fixture

from viperleed.calc.lib.context import execute_in_dir


class TestExecuteInDir:
    """Tests for the execute_in_dir context manager."""

    @fixture(name='test_dir')
    def fixture_test_dir(self, tmp_path):
        """Create a temporary test directory."""
        test_dir = tmp_path / 'test_dir'
        test_dir.mkdir()
        return test_dir

    def test_always_reverts_on_error(self, test_dir):
        """Check exceptions in code block still revert to the original path."""
        original_dir = Path.cwd()
        with pytest.raises(BaseException, match='Forced error'):
            with execute_in_dir(test_dir):
                assert Path.cwd() == test_dir
                # About the disable: we really want to make sure that
                # ANY exception is caught, not only a specific one.
                # pylint: disable-next=broad-exception-raised
                raise BaseException('Forced error')
        # Ensure it reverts back even after an error
        assert Path.cwd() == original_dir

    def test_existing_directory(self, test_dir):
        """Test execution in an existing directory."""
        original_dir = Path.cwd()
        with execute_in_dir(test_dir):
            assert Path.cwd() == test_dir
        assert Path.cwd() == original_dir  # Ensure it reverts back

    def test_non_existing_directory(self, tmp_path):
        """Check complaints when called with a non-existing directory."""
        non_existing_dir = tmp_path / 'does_not_exist'
        assert not non_existing_dir.exists()
        raises = pytest.raises(ValueError,
                               match='Directory .* does not exist.*mkdir=True')
        with raises:
            with execute_in_dir(non_existing_dir, mkdir=False):
                pytest.fail('Should not reach here')

    def test_non_existing_directory_with_mkdir(self, tmp_path):
        """Test creation of a new directory when mkdir=True."""
        new_dir = tmp_path / 'new_dir'
        assert not new_dir.exists()
        with execute_in_dir(new_dir, mkdir=True):
            assert Path.cwd() == new_dir
            assert new_dir.exists()   # Directory should be created
        assert Path.cwd() != new_dir  # Ensure it reverts back

    def test_not_a_directory(self, tmp_path):
        """Test that an error is raised if the path is not a directory."""
        file_path = tmp_path / 'not_a_directory.txt'
        file_path.touch()  # Create an empty file
        with pytest.raises(ValueError, match='is not a directory'):
            with execute_in_dir(file_path):
                pytest.fail('Should not reach here')

    def test_path_changed_only_in_context(self, test_dir):
        """Ensure that creating the context maintains the path unchanged."""
        here = Path.cwd()
        context = execute_in_dir(test_dir)
        assert Path.cwd() == here
        with context:
            assert Path.cwd() == test_dir
