"""Tests for helper functions of a multi-domain initialization."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-07'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture

from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.initialization import _make_domain_workdir


class TestMakeDomainWork:
    """Tests for the _make_domain_workdir function."""

    @fixture(name='make_work')
    def fixture_make_workdir(self, tmp_path):
        """Call _make_domain_workdir in tmp_path."""
        def _make(*args):
            with execute_in_dir(tmp_path):
                return _make_domain_workdir(*args)
        return _make

    def test_already_exists(self, make_work, tmp_path, caplog):
        """Check successful execution when the work path exists already."""
        name = 'some_domain'
        expect_work = tmp_path / f'Domain_{name}'
        expect_work.mkdir()
        work = make_work(name)
        assert 'already exists' in caplog.text

    def test_mkdir_fails(self, make_work, mocker):
        """Check exceptions raised by mkdir are not caught."""
        mock_mkdir = mocker.patch('pathlib.Path.mkdir', side_effect=OSError)
        with pytest.raises(OSError):
            make_work('fails')
        mock_mkdir.assert_called_once()

    def test_success(self, make_work):
        """Check successful creation of a domain work directory."""
        name = 'some_domain'
        work = make_work(name)
        assert work.is_dir()
        assert work.name == f'Domain_{name}'
