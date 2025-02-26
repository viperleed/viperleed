"""Tests for cleanup.cleanup (and helpers) of viperleed.calc.section."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-03'
__license__ = 'GPLv3+'

import re

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.domain_params import DomainParameters
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.files.manifest import ManifestFile
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.cleanup import cleanup
from viperleed.calc.sections.cleanup import _organize_all_work_directories
from viperleed.calc.sections.cleanup import _write_final_log_messages
from viperleed.calc.sections.cleanup import _write_manifest_file

from ....helpers import CustomTestException
from ....helpers import raises_test_exception
from .conftest import _MODULE


class TestCleanup:
    """Tests for the cleanup function."""

    mocked = (  # All the stuff that is called inside cleanup
        '_LOGGER.info',
        '_organize_all_work_directories',
        '_write_manifest_file',
        '_write_final_log_messages',
        'close_all_handlers',
        'logging.shutdown',
        )

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details of cleanup with mocks."""
        return {f: mocker.patch(f'{_MODULE}.{f}') for f in self.mocked}

    @fixture(name='manifest')
    def mock_manifest(self):
        """Return the contents of a sample manifest file."""
        return ManifestFile('file1.txt', 'file2.txt', 'dir1')

    def test_no_rpars(self, manifest, mock_implementation, mocker):
        """Check calls when no Rparams is passed."""
        # Create a "singleton" that we can use to check that
        # cleanup creates an empty Rparams in this case
        fake_rpars = mocker.MagicMock(spec=Rparams)
        mocker.patch(f'{_MODULE}.Rparams', return_value=fake_rpars)

        cleanup(manifest)
        calls = {
            '_LOGGER.info': '\nStarting cleanup...',
            '_organize_all_work_directories': fake_rpars,
            '_write_manifest_file': fake_rpars,
            '_write_final_log_messages': fake_rpars,
            'close_all_handlers': None,
            'logging.shutdown': None,
            }
        for func, arg in calls.items():
            mock = mock_implementation[func]
            if arg is None:
                mock.assert_called_once()
            else:
                mock.assert_called_once_with(arg)

    @parametrize(mock_name=mocked)
    def test_raises(self, mock_name, rpars, mock_implementation):
        """Check complaints when the implementation raises exceptions."""
        mock = mock_implementation[mock_name]
        mock.side_effect = CustomTestException
        with pytest.raises(CustomTestException):
            cleanup(rpars)

    def test_success(self, rpars, mock_implementation):
        """Check a successful execution of cleanup."""
        cleanup(rpars)
        calls = {
            '_LOGGER.info': '\nStarting cleanup...',
            '_organize_all_work_directories': rpars,
            '_write_manifest_file': rpars,
            '_write_final_log_messages': rpars,
            'close_all_handlers': None,
            'logging.shutdown': None,
            }
        for func, arg in calls.items():
            mock = mock_implementation[func]
            if arg is None:
                mock.assert_called_once()
            else:
                mock.assert_called_once_with(arg)

    @parametrize(arg=(None, 'str', set(), {}, tuple()))
    def test_typeerror(self, arg):
        """Check complaints for an invalid argument type."""
        with pytest.raises(TypeError):
            cleanup(arg)


class TestOrganizeAllWorkDirectories:
    """Tests for the _organize_all_work_directories function."""

    @fixture(name='organize')
    def mock_organize_workdir(self, mocker):
        """Replace the organize_workdir function with a mock."""
        return mocker.patch(f'{_MODULE}.organize_workdir')

    def test_domains(self, rpars, organize, mocker):
        """Check the expected calls for a multi-domain calculation."""
        rpars.closePdfReportFigs = mocker.MagicMock()
        main_call = mocker.call(rpars=rpars,
                                path='',
                                delete_unzipped=True,
                                tensors=False,
                                deltas=False)
        # Add some fake domains
        domain_one = mocker.MagicMock()
        domain_one.workdir = mocker.MagicMock()
        domain_one.rpars = Rparams()  # No Tensors/Deltas in manifest
        call_one = mocker.call(rpars=domain_one.rpars,
                               path=domain_one.workdir,
                               delete_unzipped=True,
                               tensors=False,
                               deltas=False)
        domain_two = mocker.MagicMock()
        domain_two.workdir = mocker.MagicMock()
        domain_two.rpars = Rparams()
        domain_two.rpars.manifest = {DEFAULT_TENSORS, DEFAULT_DELTAS}
        call_two = mocker.call(rpars=domain_two.rpars,
                               path=domain_two.workdir,
                               delete_unzipped=True,
                               tensors=True,
                               deltas=True)
        rpars.domainParams = domain_one, domain_2
        calls = main_call, call_one, call_two

        _organize_all_work_directories(rpars)
        rpars.closePdfReportFigs.assert_called_once()
        assert organize.call_count == len(calls)
        organize.assert_has_calls(calls)

    def test_no_domains(self, rpars, organize, mocker):
        """Check the expected calls for a single-domain calculation."""
        rpars.closePdfReportFigs = mocker.MagicMock()
        rpars.manifest.add(DEFAULT_TENSORS)  # But not deltas

        kwargs = {
            'rpars': rpars,
            'path': '',
            'delete_unzipped': True,
            'tensors': True,  # Added to manifest
            'deltas': False,
            }
        _organize_all_work_directories(rpars)
        organize.assert_called_once_with(**kwargs)

    def test_rpars_empty(self, rpars, organize):
        """Check the expected calls with an empty Rparams."""
        _organize_all_work_directories(rpars)
        kwargs = {
            'rpars': rpars,
            'path': '',
            'delete_unzipped': True,
            'tensors': False,
            'deltas': False,
            }
        organize.assert_called_once_with(**kwargs)


class TestWriteFinalLogMessages:
    """Tests for the _write_final_log_messages function."""

    @staticmethod
    def check_has_records(caplog, records):
        """Ensure caplog has only certain records."""
        logged = tuple(r.getMessage() for r in caplog.records)
        print(logged)
        print(records)
        assert len(logged) == len(records)
        for log, expect in zip(logged, records):
            if isinstance(expect, str):
                assert log == expect
            else:
                assert expect.match(log)

    @fixture(name='rpars_filled')
    def fixture_rpars_filled(self, rpars, mocker):
        """Return an Rparams with relevant attributes set."""
        rpars.timer = mocker.MagicMock()
        rpars.timer.how_long.return_value = '1h 30m'
        rpars.runHistory = [1, 2, 3]
        rpars.stored_R = {'refcalc': [0.5, 0.4, 0.3],
                          'superpos': None}
        rpars.checklist = ['Check convergence', 'Verify inputs']
        return rpars

    def test_crashed_early(self, rpars, caplog):
        """Check logging messages when cleanup is called early."""
        caplog.set_level(0)  # All messages
        rpars.timer = None   # Like cleanup if called with None
        _write_final_log_messages(rpars)
        expect = (
            re.compile(r'\nFinishing execution at .*'
                       r'\nTotal elapsed time: unknown\n'),
            '',
            )
        self.check_has_records(caplog, expect)

    def test_domains(self, rpars_filled, caplog):
        """Check logging messages for a multi-domain calculation."""
        caplog.set_level(0)  # All messages
        rpars_filled.domainParams.append(1.234)
        _write_final_log_messages(rpars_filled)
        expect = (
            re.compile(r'\nFinishing execution at .*'
                       r'\nTotal elapsed time: 1h 30m\n'),
            'Executed segments: 1 2 3',
            'Final R (refcalc): 0.5000 (0.4000 / 0.3000)',
            re.compile(r'Domain.*bookkeeper.*--archive.'),
            '',
            '# The following issues should be checked before starting again:',
            '- Check convergence',
            '- Verify inputs',
            '',
            )
        self.check_has_records(caplog, expect)

    def test_no_checklist(self, rpars_filled, caplog):
        """Check that no checklist-related messages are emitted."""
        caplog.set_level(0)  # All messages
        rpars_filled.checklist = []
        rpars_filled.stored_R['refcalc'] = [0.5, 0, 0.5]
        _write_final_log_messages(rpars_filled)
        expect = (
            re.compile(r'\nFinishing execution at .*'
                       r'\nTotal elapsed time: 1h 30m\n'),
            'Executed segments: 1 2 3',
            'Final R (refcalc): 0.5000',
            '',
            )
        self.check_has_records(caplog, expect)

    def test_no_domains(self, rpars_filled, caplog):
        """Check logging messages for a single-domain calculation."""
        caplog.set_level(0)  # All messages
        _write_final_log_messages(rpars_filled)
        expect = (
            re.compile(r'\nFinishing execution at .*'
                       r'\nTotal elapsed time: 1h 30m\n'),
            'Executed segments: 1 2 3',
            'Final R (refcalc): 0.5000 (0.4000 / 0.3000)',
            '',
            '# The following issues should be checked before starting again:',
            '- Check convergence',
            '- Verify inputs',
            '',
            )
        self.check_has_records(caplog, expect)


class TestWriteManifest:
    """Tests for the _write_manifest_file function."""

    def test_fails(self, rpars, tmp_path, mocker, caplog):
        """Test logging when opening manifest fails."""
        with execute_in_dir(tmp_path):
            rpars.manifest = manifest = ManifestFile()
            mocker.patch.object(manifest, 'write', side_effect=OSError)
            _write_manifest_file(rpars)
        expect_log = 'Failed to write manifest file.'
        assert not (tmp_path/'manifest').is_file()
        assert expect_log in caplog.text

    def test_collects_domains(self, rpars, tmp_path):
        """Check that contents of domains are collected."""
        rpars.manifest = manifest = ManifestFile(path=tmp_path)
        for i in range(5):
            d_path = tmp_path/f'dd_{i}'
            d_rpars = Rparams()
            d_rpars.manifest = ManifestFile(path=d_path)
            domain = DomainParameters(d_path, f'name {i}')
            domain.rpars = d_rpars
            rpars.domainParams.append(domain)
        with execute_in_dir(tmp_path):
            _write_manifest_file(rpars)
        # Check that now all sub-manifests are present
        subpaths = set(manifest.paths)
        expect_paths = (
            tmp_path,
            *(d.rpars.manifest.path for d in rpars.domainParams),
            )
        assert all(p in subpaths for p in expect_paths)

        # And verify that contents have been written correctly
        contents = (manifest.path/'manifest').read_text().splitlines()
        headers = [line for line in contents if line.startswith('[')]
        expect_headers = [f'[domain name {i} at dd_{i}]' for i in range(5)]
        assert headers == expect_headers

    def test_raises(self, rpars):
        """Check that only OSError is caught cleanly."""
        rpars.manifest = manifest = ManifestFile()
        with raises_test_exception(manifest, 'write'):
            _write_manifest_file(rpars)

    _success = {
        'unique': ('file1.txt', 'file2.txt', 'file3.txt'),
        'duplicates': ('file1.txt', 'file2.txt', 'file1.txt'),
        }

    @parametrize(contents=_success.values(), ids=_success)
    def test_success(self, rpars, contents, tmp_path, caplog):
        """Check successful writing to the manifest file."""
        caplog.set_level(0)  # All messages
        with execute_in_dir(tmp_path):
            rpars.manifest = manifest = ManifestFile(*contents)
            _write_manifest_file(rpars)
        manifest = tmp_path/'manifest'
        assert manifest.is_file()

        # Check contents
        written = sorted(line
                         for line in manifest.read_text().splitlines()
                         if line)
        assert written == sorted(set(contents))

        # Check logging
        expect_log = 'Wrote manifest file successfully.'
        assert expect_log in caplog.text
