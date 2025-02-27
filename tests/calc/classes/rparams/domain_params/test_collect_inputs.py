"""Tests for the collection of domain inputs in section initialization."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-19'
__license__ = 'GPLv3+'

import logging
from pathlib import Path
from zipfile import BadZipFile

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.domain_params import _DOMAIN_INPUT_FILES

from .....helpers import filesystem_from_dict
from .....helpers import filesystem_to_dict

_MODULE = 'viperleed.calc.classes.rparams.domain_params'


@fixture(name='domain')
def fixture_domain(make_domain, tmp_path):
    """Return an initialized DomainParameters."""
    return make_domain(tmp_path/'fake_workdir', 'fake_domain')


@fixture(name='src_dir')
def fixture_src_dir(tmp_path):
    """Return a temporary source directory for domain inputs."""
    src = tmp_path / 'source_dir'
    src.mkdir()
    return src


@fixture(name='src_tensor')
def fixture_src_tensor(tmp_path):
    """Return the path to a temporary Tensor file as domain source."""
    src = tmp_path / 'source_tensor.zip'
    src.touch()
    return src


class TestCollectInputsForDomain:
    """Tests for the collect_input_files method."""

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, domain, mocker):
        """Replace implementation details of collect_input_files."""
        def _mock(is_tensor=False, is_dir=False,
                  dir_raises=None, tensor_raises=None):
            mocker.patch('pathlib.Path.is_file', return_value=is_tensor)
            mocker.patch('pathlib.Path.is_dir', return_value=is_dir)
            return {
            'dir': mocker.patch.object(domain,
                                       '_collect_inputs_from_directory',
                                       side_effect=dir_raises),
            'tensor' : mocker.patch.object(domain,
                                           '_collect_inputs_from_tensor_file',
                                           side_effect=tensor_raises),
            }
        return _mock

    @fixture(name='collect')
    def fixture_collect(self, domain, mock_implementation, caplog):
        """Call collect_input_files at a given path."""
        def _call(src, **kwargs):
            caplog.set_level(logging.INFO)
            mocks = mock_implementation(**kwargs)
            domain.collect_input_files(src)
            return mocks
        return _call

    @staticmethod
    def check_log_written(caplog):
        """Ensure the expected log messages are emitted."""
        expect_log = 'Fetching input files for domain fake_domain'
        assert expect_log in caplog.text

    def test_fetch_folder(self, collect, src_dir, caplog):
        """Check dispatch to _collect_inputs_from_directory."""
        mocks = collect(src_dir, is_dir=True)
        mocks['dir'].assert_called_once_with(src_dir)
        mocks['tensor'].assert_not_called()
        self.check_log_written(caplog)

    def test_fetch_folder_fails(self, collect, src_dir):
        """Check complaints when collection fails."""
        _raises = pytest.raises(RuntimeError,
                                match='Error getting domain input files')
        with _raises:
            collect(src_dir, is_dir=True, dir_raises=OSError)

    def test_fetch_tensor(self, collect, src_tensor, caplog):
        """Check dispatch to _collect_inputs_from_tensor_file."""
        mocks = collect(src_tensor, is_tensor=True)
        mocks['tensor'].assert_called_once_with(src_tensor)
        mocks['dir'].assert_not_called()
        self.check_log_written(caplog)

    def test_fetch_nowhere(self, collect, caplog):
        """Check no dispatch if src is neither a folder nor a directory."""
        mocks = collect(Path.cwd())
        for mock in mocks.values():
            mock.assert_not_called()
        self.check_log_written(caplog)

    def test_fetch_tensor_fails(self, collect, src_tensor):
        """Check complaints when collection fails."""
        _raises = pytest.raises(RuntimeError,
                                match='Error getting domain input files')
        with _raises:
            collect(src_tensor, is_tensor=True, tensor_raises=OSError)


class TestCollectFromDirectory:
    """Tests for the _collect_inputs_from_directory method."""

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details with mocks."""
        def _mock(max_index=1, mock_copy=True, copy_raises=None):
            mocks = {
                'index': mocker.patch(f'{_MODULE}.leedbase.getMaxTensorIndex',
                                      return_value=max_index),
                'get_tensor': mocker.patch(
                    f'{_MODULE}.iotensors.fetch_unpacked_tensor'
                    ),
                }
            if mock_copy or copy_raises:
                mocks['copy'] = mocker.patch('shutil.copy2',
                                             side_effect=copy_raises)
            return mocks
        return _mock

    @fixture(name='collect')
    def fixture_collect(self, domain, src_dir, mock_implementation):
        """Call _collect_inputs_from_directory."""
        def _call(mock=True, **kwargs):
            mocks = mock_implementation(**kwargs) if mock or kwargs else {}
            # pylint: disable-next=protected-access       # OK in tests
            domain._collect_inputs_from_directory(src_dir)
            return mocks
        return _call

    @parametrize(exc=(OSError, FileNotFoundError, BadZipFile))
    def test_input_file_fails(self, exc, collect):
        """Check complaints when failing to fetch an input file."""
        with pytest.raises(exc):
            collect(max_index=0, copy_raises=exc)

    @parametrize(exc=(OSError, FileNotFoundError))
    def test_phaseshifts_fails(self, exc, domain, collect):
        """Check no complaints if PHASESHIFTS can't be fetched from src_dir."""
        def _copy_raises(src, _):
            # pylint: disable-next=magic-value-comparison
            if src.name == 'PHASESHIFTS':
                raise exc
        mocks = collect(max_index=0, copy_raises=_copy_raises)
        assert domain.refcalc_required  # No exception
        assert mocks['copy'].call_count == len(_DOMAIN_INPUT_FILES)

    def test_refcalc_required(self, domain, collect):
        """Check successful fetching of inputs when no Tensor is found."""
        mocks = collect(max_index=0)
        mocks['get_tensor'].assert_not_called()
        assert domain.refcalc_required
        assert mocks['copy'].call_count == len(_DOMAIN_INPUT_FILES)

    def test_tensor_found(self, domain, collect):
        """Check successful collection when a Tensor file is found."""
        # Mock the result of iotensors.fetch_unpacked_tensor
        tensor_input_files = _DOMAIN_INPUT_FILES + ('IVBEAMS',)
        unzipped_tensor = {'Tensors': {'Tensors_001': {
            f: None for f in tensor_input_files
            }}}
        filesystem_from_dict(unzipped_tensor, domain.workdir)
        mocks = collect()
        assert not domain.refcalc_required
        for mock in mocks.values():
            mock.assert_called()
        assert mocks['copy'].call_count == len(tensor_input_files)

    @parametrize(exc=(OSError, BadZipFile))
    def test_tensor_found_but_fetch_fails(self,
                                          exc,
                                          domain,
                                          mock_implementation,
                                          collect):
        """Check fallback to source when failing to fetch a Tensor file."""
        mocks = mock_implementation()
        mocks['get_tensor'].side_effect = exc
        collect(mock=False)
        assert domain.refcalc_required
        for mock in mocks.values():
            mock.assert_called()
        assert mocks['copy'].call_count == len(_DOMAIN_INPUT_FILES)

    def test_tensor_found_but_file_missing(self, domain, src_dir, collect):
        """Check fallback to source if input files in Tensor are missing."""
        # No inputs in the unzipped tensor...
        unzipped_tensor = {'Tensors': {'Tensors_001': {}}}
        filesystem_from_dict(unzipped_tensor, domain.workdir)
        # ... but they're there in the source directory
        src_inputs = dict.fromkeys(_DOMAIN_INPUT_FILES, '')
        filesystem_from_dict(src_inputs, src_dir)
        mocks = collect(mock_copy=False)
        assert domain.refcalc_required
        for mock in mocks.values():
            mock.assert_called()
        work_tree = filesystem_to_dict(domain.workdir)
        assert work_tree == {**unzipped_tensor, **src_inputs}


class TestCollectFromZip:
    """Tests for the _collect_inputs_from_tensor_file method."""

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details with mocks."""
        def _mock(mock_copy=True):
            mocks = {
                'index': mocker.patch(f'{_MODULE}.leedbase.getMaxTensorIndex',
                                      return_value=1),
                'zip': mocker.patch(f'{_MODULE}.iotensors.unpack_tensor_file'),
                }
            if mock_copy:
                mocks['copy'] = mocker.patch('shutil.copy2')
            return mocks
        return _mock

    @fixture(name='collect')
    def fixture_collect(self, mock_implementation, domain, src_tensor):
        """Call _collect_inputs_from_tensor_file."""
        def _call(mock=True, **kwargs):
            mocks = mock_implementation(**kwargs) if mock else {}
            # pylint: disable-next=protected-access       # OK in tests
            domain._collect_inputs_from_tensor_file(src_tensor)
            return mocks
        return _call

    def test_success(self, collect, domain):
        """Check the successful collection of inputs from a Tensor file."""
        # Fake the result of the patched unpack_tensor_file
        tensor_input_files = _DOMAIN_INPUT_FILES + ('IVBEAMS',)
        unzipped_tensor = {'Tensors': {'Tensors_002': {
            f: None for f in tensor_input_files
            }}}
        filesystem_from_dict(unzipped_tensor, domain.workdir)
        mocks = collect()
        for mock in mocks.values():
            mock.assert_called()
        assert mocks['copy'].call_count == len(tensor_input_files)

    def test_inputs_missing(self, collect):
        """Check complaints when no Tensors exist."""
        with pytest.raises(FileNotFoundError):
            collect(mock_copy=False)

    @parametrize(exc=(OSError, BadZipFile))
    def test_unzip_fails(self, exc, mock_implementation, collect):
        """Check complaints if unzipping Tensors fails."""
        mocks = mock_implementation()
        mocks['zip'].side_effect = exc
        with pytest.raises(exc):
            collect(mock=False)
