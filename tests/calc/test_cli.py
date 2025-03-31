"""Tests for the viperleed.calc command-line interface."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

import os

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.cli import ViPErLEEDCalcCLI
from viperleed.calc.cli import _copy_files_from_manifest
from viperleed.calc.lib.context import execute_in_dir

from ..helpers import filesystem_from_dict
from ..helpers import filesystem_to_dict


@fixture(name='calc_parser')
def fixture_calc_parser():
    """Return a CLI argument parser for viperleed.calc."""
    return ViPErLEEDCalcCLI().parser


class TestCalcParser:
    """Tests for parsing CLI arguments of viperleed.calc."""

    def test_parse_version(self, calc_parser):
        """Check that requesting the version exits afterwards."""
        with pytest.raises(SystemExit):
            calc_parser.parse_args(['--version'])

    def test_parse_work(self, calc_parser, tmp_path):
        """Check interpretation of -w flag."""
        parsed = calc_parser.parse_args(['-w', str(tmp_path)])
        assert parsed.work == str(tmp_path)

    @parametrize(v_flag=('-v', '--verbose'))
    def test_parse_verbose(self, calc_parser, v_flag):
        """Check interpretation of -v flag."""
        assert calc_parser.parse_args([v_flag,]).verbose

    @parametrize(v_flag=('-vv', '--very-verbose'))
    def test_parse_very_verbose(self, calc_parser, v_flag):
        """Check interpretation of -vv flag."""
        assert calc_parser.parse_args([v_flag,]).very_verbose


class TestCopyFilesFromManifest:
    """Tests for the _copy_files_from_manifest helper function."""

    @fixture(name='manifest')
    def fixture_manifest(self, tmp_path):
        """Create a manifest file and its contents at `tmp_path`."""
        manifest = 'file1.txt \nfile2  \n  \n\n  folder\n'
        manifest += '''
[folder at one]
one/subfile.txt

[Domain 3 at two]
two/domain_file
'''
        copied = {
            'file1.txt': 'Test file 1',
            'file2': 'Test file 2',
            'folder': {},
            'one': {'subfile.txt': 'Subfile contents'},
            'two': {'domain_file': 'domain file contents'},
            }
        stay = {
            'manifest': manifest,
            'file_not_in_manifest': None,
            'folder_not_in_manifest': {}
            }
        tree = {**copied, **stay}
        filesystem_from_dict(tree, tmp_path)
        return copied, stay

    def test_no_manifest_file(self, tmp_path):
        """Check that no resources are copied if manifest does not exist."""
        dest = tmp_path/'dest'
        dest.mkdir()
        with execute_in_dir(tmp_path):
            _copy_files_from_manifest(dest)
        copied = filesystem_to_dict(dest)
        assert not copied

    def test_copy_successful(self, manifest, tmp_path):
        """Check the successful copy of files/folders."""
        dest = tmp_path/'dest'
        dest.mkdir()
        with execute_in_dir(tmp_path):
            _copy_files_from_manifest(dest)
        expect_copy, expect_stay = manifest
        copied = filesystem_to_dict(dest)
        assert copied == expect_copy
        assert not any(s in copied for s in expect_stay)

    def test_copy_fails(self, manifest, tmp_path, capsys):
        """Check complaints are printed if copying resources fails."""
        manifest_file = tmp_path/'manifest'
        with manifest_file.open('a', encoding='utf-8') as file:
            file.write('two/this_does_not_exist\n')
        self.test_copy_successful(manifest, tmp_path)
        expect_print = f'Error copying two{os.sep}this_does_not_exist'
        stdout = capsys.readouterr().out
        assert expect_print in stdout
