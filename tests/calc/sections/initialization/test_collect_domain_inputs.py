"""Tests for the collection of domain inputs in section initialization."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-19'
__license__ = 'GPLv3+'

import logging
from pathlib import Path

import pytest
from pytest_cases import fixture

from viperleed.calc.classes.rparams.domain_params import DomainParameters
from viperleed.calc.sections.initialization import (
    _DOMAIN_INPUT_FILES,
    _collect_inputs_for_domain,
    _collect_inputs_for_domain_from_directory,
    _collect_inputs_for_domain_from_tensor_file,
    )

from ....helpers import filesystem_from_dict

_MODULE = 'viperleed.calc.sections.initialization'


@fixture(name='domain')
def fixture_domain(tmp_path):
    """Return an initialized DomainParameters."""
    return DomainParameters(workdir=tmp_path/'fake_workdir',
                            name='fake_domain')


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
    """Tests for the _collect_inputs_for_domain function."""

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details of _collect_inputs_for_domain."""
        return {
        'dir': mocker.patch(
            f'{_MODULE}._collect_inputs_for_domain_from_directory',
            ),
        'tensor' : mocker.patch(
            f'{_MODULE}._collect_inputs_for_domain_from_tensor_file',
            ),
        }

    @fixture(name='collect')
    def fixture_collect(self, domain, mock_implementation, caplog):
        """Call _collect_inputs_for_domain at a given path."""
        def _call(src):
            caplog.set_level(logging.INFO)
            _collect_inputs_for_domain(domain, src)
            return mock_implementation
        return _call

    @staticmethod
    def check_log_written(caplog):
        """Ensure the expected log messages are emitted."""
        expect_log = 'Fetching input files for domain fake_domain'
        assert expect_log in caplog.text

    # pylint: disable-next=too-many-arguments  # All fixtures
    def test_fetch_folder(self, collect, domain, src_dir, mocker, caplog):
        """Check dispatch to _collect_inputs_for_domain_from_directory."""
        mocker.patch('pathlib.Path.is_dir', return_value=True)
        mocker.patch('pathlib.Path.is_file', return_value=False)
        mocks = collect(src_dir)
        mocks['dir'].assert_called_once_with(domain, src_dir)
        mocks['tensor'].assert_not_called()
        self.check_log_written(caplog)

    # pylint: disable-next=too-many-arguments  # All fixtures
    def test_fetch_tensor(self, collect, domain, src_tensor, mocker, caplog):
        """Check dispatch to _collect_inputs_for_domain_from_tensor_file."""
        mocker.patch('pathlib.Path.is_dir', return_value=False)
        mocker.patch('pathlib.Path.is_file', return_value=True)
        mocks = collect(src_tensor)
        mocks['tensor'].assert_called_once_with(domain, src_tensor)
        mocks['dir'].assert_not_called()
        self.check_log_written(caplog)

    def test_fetch_nowhere(self, collect, mocker, caplog):
        """Check no dispatch if src is neither a folder nor a directory."""
        mocker.patch('pathlib.Path.is_dir', return_value=False)
        mocker.patch('pathlib.Path.is_file', return_value=False)
        mocks = collect(Path.cwd())
        for mock in mocks.values():
            mock.assert_not_called()
        self.check_log_written(caplog)


class TestCollectFromDirectory:
    """Tests for the _collect_inputs_for_domain_from_directory function."""

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details with mocks."""
        def _mock(max_index=1):
            return {
                'index': mocker.patch(f'{_MODULE}.leedbase.getMaxTensorIndex',
                                      return_value=max_index),
                'get_tensor': mocker.patch(f'{_MODULE}.iotensors.getTensors'),
                'copy': mocker.patch('shutil.copy2'),
                }
        return _mock

    @fixture(name='collect')
    def fixture_collect(self, domain, src_dir, mock_implementation):
        """Call _collect_inputs_for_domain_from_directory."""
        def _call(mock=True, **kwargs):
            mocks = mock_implementation(**kwargs) if mock or kwargs else {}
            _collect_inputs_for_domain_from_directory(domain, src_dir)
            return mocks
        return _call

    def test_copy_file_fails(self, mock_implementation, collect, mocker):
        """Check complaints when failing to fetch an input file."""
        mocker.patch('pathlib.Path.is_file', return_value=True)
        mocks = mock_implementation(max_index=0)
        mocks['copy'].side_effect = OSError
        with pytest.raises(RuntimeError):
            collect(mock=False)

    def test_copy_phaseshifts_fails(self, domain, mock_implementation,
                                    collect, mocker):
        """Check complaints when failing to fetch an input file."""
        def _ps_copy_fails(src, _):
            # pylint: disable-next=magic-value-comparison
            if src.name == 'PHASESHIFTS':
                raise OSError
        mocks = mock_implementation(max_index=0)
        mocks['copy'] = mocker.patch('shutil.copy2',
                                     side_effect=_ps_copy_fails)
        mocker.patch('pathlib.Path.is_file', return_value=True)
        collect(mock=False)
        assert domain.refcalcRequired  # No exception
        assert mocks['copy'].call_count == len(_DOMAIN_INPUT_FILES)

    def test_input_file_missing(self, collect):
        """Check complaints when failing to fetch an input file."""
        with pytest.raises(RuntimeError):
            collect(max_index=0)

    def test_phaseshifts_missing(self, domain, collect, mocker):
        """Check complaints when failing to fetch an input file."""
        mocker.patch('pathlib.Path.is_file',
                     # pylint: disable-next=magic-value-comparison
                     lambda p: p.name != 'PHASESHIFTS')
        collect(max_index=0)
        assert domain.refcalcRequired  # No exception

    def test_refcalc_required(self, domain, collect, mocker):
        """Check successful fetching of inputs when no Tensor is found."""
        mocker.patch('pathlib.Path.is_file', return_value=True)
        mocks = collect(max_index=0)
        mocks['get_tensor'].assert_not_called()
        assert not domain.tensorDir
        assert domain.refcalcRequired
        assert mocks['copy'].call_count == len(_DOMAIN_INPUT_FILES)

    def test_tensor_found(self, domain, collect):
        """Check successful collection when a Tensor file is found."""
        # getTensors unzips the tensor at domain.workdir
        tensor_input_files = _DOMAIN_INPUT_FILES + ('IVBEAMS',)
        unzipped_tensor = {'Tensors': {'Tensors_001': {
            f: None for f in tensor_input_files
            }}}
        filesystem_from_dict(unzipped_tensor, domain.workdir)
        mocks = collect()
        assert domain.tensorDir == domain.workdir/'Tensors'/'Tensors_001'
        assert not domain.refcalcRequired
        for mock in mocks.values():
            mock.assert_called()
        assert mocks['copy'].call_count == len(tensor_input_files)

    def test_tensor_found_but_fetch_fails(self,
                                          domain,
                                          mock_implementation,
                                          collect,
                                          mocker):
        """Check fallback to source when failing to fetch a Tensor file."""
        mocks = mock_implementation()
        mocks['get_tensor'].side_effect = Exception
        mocker.patch('pathlib.Path.is_file', return_value=True)
        collect(mock=False)
        assert not domain.tensorDir
        assert domain.refcalcRequired
        for mock in mocks.values():
            mock.assert_called()
        assert mocks['copy'].call_count == len(_DOMAIN_INPUT_FILES)

    def test_tensor_found_but_file_missing(self, domain, src_dir, collect):
        """Check fallback to source if input files in Tensor are missing."""
        # No inputs in the unzipped tensor...
        unzipped_tensor = {'Tensors': {'Tensors_001': {}}}
        filesystem_from_dict(unzipped_tensor, domain.workdir)
        # ... but they're there in the source directory
        src_inputs = dict.fromkeys(_DOMAIN_INPUT_FILES)
        filesystem_from_dict(src_inputs, src_dir)
        mocks = collect()
        assert not domain.tensorDir
        assert domain.refcalcRequired
        for mock in mocks.values():
            mock.assert_called()
        assert mocks['copy'].call_count == len(_DOMAIN_INPUT_FILES)


class TestCollectFromZip:
    """Tests for the _collect_inputs_for_domain_from_tensor_file function."""

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details with mocks."""
        def _mock():
            return {
                'index': mocker.patch(f'{_MODULE}.leedbase.getMaxTensorIndex',
                                      return_value=1),
                'zip': mocker.patch(f'{_MODULE}.ZipFile'),
                'copy': mocker.patch('shutil.copy2'),
                }
        return _mock

    @fixture(name='collect')
    def fixture_collect(self, mock_implementation, domain, src_tensor):
        """Call _collect_inputs_for_domain_from_tensor_file."""
        def _call(mock=True):
            mocks = mock_implementation() if mock else {}
            _collect_inputs_for_domain_from_tensor_file(domain, src_tensor)
            return mocks
        return _call

    def test_success(self, collect, domain):
        """Check the successful collection of inputs from a Tensor file."""
        # Fake the result of the patched ZipFile.extractall
        tensor_input_files = _DOMAIN_INPUT_FILES + ('IVBEAMS',)
        unzipped_tensor = {'Tensors': {'Tensors_002': {
            f: None for f in tensor_input_files
            }}}
        filesystem_from_dict(unzipped_tensor, domain.workdir)
        mocks = collect()
        assert domain.tensorDir == domain.workdir/'Tensors/Tensors_002'
        for mock in mocks.values():
            mock.assert_called()
        assert mocks['copy'].call_count == len(tensor_input_files)

    def test_inputs_missing(self, collect, domain):
        """Check the successful collection of inputs from a Tensor file."""
        with pytest.raises(RuntimeError):
            collect()
        assert not domain.tensorDir

    def test_max_tensor_fails(self,
                              mock_implementation,
                              collect,
                              domain,
                              mocker):
        """Check the successful collection of inputs from a Tensor file."""
        mocks = mock_implementation()
        mocks['index'].side_effect = Exception
        mocker.patch('pathlib.Path.is_file', return_value=True)
        collect(mock=False)
        assert domain.tensorDir == domain.workdir/'Tensors/Tensors_001'
        for mock in mocks.values():
            mock.assert_called()

    def test_unzip_fails(self, mock_implementation, collect, domain):
        """Check the successful collection of inputs from a Tensor file."""
        mocks = mock_implementation()
        mocks['zip'].side_effect = OSError
        with pytest.raises(RuntimeError):
            collect(mock=False)
        assert not domain.tensorDir
