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
        """Call _make_domain_workdir in a subfolder of tmp_path."""
        main_work = tmp_path/'main_work'
        main_work.mkdir()
        def _make(**kwargs):
            kwargs.setdefault('src', None)
            kwargs.setdefault('calc_started_at', None)
            with execute_in_dir(main_work):
                return _make_domain_workdir(**kwargs)
        return _make

    def test_already_exists(self, make_work, tmp_path, caplog):
        """Check successful execution when the work path exists already."""
        name = 'some_domain'
        expect_work = tmp_path / f'main_work/Domain_{name}'
        expect_work.mkdir()
        work = make_work(name=name)
        assert 'already exists' in caplog.text

    def test_mkdir_fails(self, make_work, mocker):
        """Check exceptions raised by mkdir are not caught."""
        mock_mkdir = mocker.patch('pathlib.Path.mkdir', side_effect=OSError)
        with pytest.raises(OSError):
            make_work(name='fails')
        mock_mkdir.assert_called_once()

    def test_read_from_zip(self, make_work, tmp_path):
        """Check outcome when the domain was read from a tensor ZIP."""
        domain_src = tmp_path/'domain_source.zip'
        domain_src.touch()
        name = 'some_domain'
        work = make_work(name=name, src=domain_src, calc_started_at=tmp_path)
        assert work.is_dir()
        assert work.name == f'Domain_{name}'

    def test_read_from_calc_subfolder(self, make_work, tmp_path):
        """Check outcome when the domain source is not a subfolder of calc."""
        domain_src = tmp_path/'domain_source'
        domain_src.mkdir()
        name = 'some_domain'
        work = make_work(name=name, src=domain_src, calc_started_at=tmp_path)
        assert work.is_dir()
        assert work.name == domain_src.name

    def test_read_from_somewhere_else(self, make_work, tmp_path):
        """Check outcome when the domain source is not a subfolder of calc."""
        domain_src = tmp_path/'some_other_folder'/'domain_source'
        domain_src.mkdir(parents=True)
        name = 'some_domain'
        work = make_work(name=name, src=domain_src, calc_started_at=tmp_path)
        assert work.is_dir()
        assert work.name == f'Domain_{name}'

    def test_success(self, make_work):
        """Check successful creation of a domain work directory."""
        name = 'some_domain'
        work = make_work(name=name)
        assert work.is_dir()
        assert work.name == f'Domain_{name}'
