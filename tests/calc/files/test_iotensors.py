"""Tests for the iotensors module of viperleed.calc.files."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-20'
__license__ = 'GPLv3+'

from pathlib import Path
import shutil
from zipfile import BadZipFile
from zipfile import ZipFile

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.files.iotensors import getTensors
from viperleed.calc.files.iotensors import unpack_tensor_file
from viperleed.calc.lib.context import execute_in_dir

_MODULE = 'viperleed.calc.files.iotensors'


@fixture(name='patched_zip')
def fixture_patch_zip(mocker):
    """Replace ZipFile with a fake one."""
    fake_instance = mocker.MagicMock(spec=ZipFile)
    mock_zip = mocker.patch(f'{_MODULE}.ZipFile', return_value=fake_instance)
    fake_instance.__enter__.return_value = fake_instance
    return mock_zip, fake_instance


class TestGetTensors:
    """Tests for the getTensors function."""

    @fixture(name='base_dir')
    def fixture_base_dir(self, tmp_path):
        """Return a temporary directory from which Tensors are fetched."""
        base_dir = tmp_path / 'base'
        base_dir.mkdir()
        return base_dir

    @fixture(name='target_dir')
    def mock_target_dir(self, tmp_path):
        """Return a temporary directory where Tensors should go to."""
        target_dir = tmp_path / 'target'
        target_dir.mkdir()
        return target_dir

    def test_no_folder_or_zip(self, base_dir):
        """Check complaints when no Tensor is ound at base_dir."""
        _raises = pytest.raises(FileNotFoundError,
                                match='No Tensors folder/zip file')
        with _raises:
            getTensors(1, base_dir=base_dir)

    def test_unzip_success(self, base_dir, target_dir, mocker):
        """Check successful unzipping of a Tensor file."""
        mock_unpack = mocker.patch(f'{_MODULE}.unpack_tensor_file')
        tensor_zip = base_dir / DEFAULT_TENSORS / 'Tensors_001.zip'
        unpack_path = target_dir / DEFAULT_TENSORS / 'Tensors_001'
        tensor_zip.parent.mkdir(parents=True)
        tensor_zip.touch()
        getTensors(1, base_dir=base_dir, target_dir=target_dir)
        mock_unpack.assert_called_once_with(tensor_zip, unpack_path)

    @parametrize(exc=(OSError, BadZipFile))
    def test_unzip_fails(self, exc, base_dir, target_dir, mocker):
        """Check complaints if unpacking a Tensor fails."""
        mock_unpack = mocker.patch(f'{_MODULE}.unpack_tensor_file',
                                   side_effect=exc)
        tensor_zip = base_dir / DEFAULT_TENSORS / 'Tensors_001.zip'
        tensor_zip.parent.mkdir(parents=True)
        tensor_zip.touch()
        with pytest.raises(exc):
            getTensors(1, base_dir=base_dir, target_dir=target_dir)

    def test_copy_folder_successful(self, base_dir, target_dir, mocker):
        """Test successful copying of a tensor folder."""
        mock_copy = mocker.patch(f'{_MODULE}.copytree_exists_ok')
        tensor_folder = base_dir / DEFAULT_TENSORS / 'Tensors_001'
        unpack_path = target_dir / DEFAULT_TENSORS / 'Tensors_001'
        tensor_folder.mkdir(parents=True)
        getTensors(1, base_dir=base_dir, target_dir=target_dir)
        mock_copy.assert_called_once_with(tensor_folder, unpack_path)

    def test_no_need_to_copy_folder(self, tmp_path, mocker):
        """Test that no copying is done if target == source."""
        mock_copy = mocker.patch(f'{_MODULE}.copytree_exists_ok')
        tensor_folder = tmp_path / DEFAULT_TENSORS / 'Tensors_001'
        tensor_folder.mkdir(parents=True)
        with execute_in_dir(tmp_path):
            getTensors(1)
        mock_copy.assert_not_called()

    def test_copy_folder_fails(self, base_dir, target_dir, mocker):
        """Check complaints if copying a tensor folder fails."""
        mock_copy = mocker.patch(f'{_MODULE}.copytree_exists_ok',
                                 side_effect=OSError('Copy failed'))
        tensor_folder = base_dir / DEFAULT_TENSORS / 'Tensors_001'
        tensor_folder.mkdir(parents=True)
        with pytest.raises(OSError, match='Copy failed'):
            getTensors(1, base_dir=base_dir, target_dir=target_dir)

    def test_basedir_is_tensors(self, base_dir, target_dir, mocker):
        """Check correct behavior when base_dir is the Tensors folder."""
        mock_copy = mocker.patch(f'{_MODULE}.copytree_exists_ok')
        base_tensors = base_dir / DEFAULT_TENSORS
        tensor_folder = base_tensors / 'Tensors_001'
        tensor_folder.mkdir(parents=True)
        unpack_path = target_dir / DEFAULT_TENSORS / 'Tensors_001'
        getTensors(1, base_dir=base_tensors, target_dir=target_dir)
        mock_copy.assert_called_once_with(tensor_folder, unpack_path)


class TestUnpackTensorFile:
    """Tests for the unpack_tensor_file function."""

    @fixture(name='call')
    def fixture_call(self, tmp_path):
        """Call unpack_tensor_file on a tensor."""
        def _call(tensor_name='Tensor.zip'):
            tensor_zip = tmp_path/tensor_name
            unpack_path = tmp_path/'unpack'
            tensor_zip.touch()
            unpack_tensor_file(tensor_zip, unpack_path)
            assert unpack_path.is_dir()
            return tensor_zip, unpack_path
        return _call

    def test_fails_to_open_archive(self, call):
        """Check complaints if opening a Tensor zip fails."""
        with pytest.raises(BadZipFile):
            call('not_a_zip.zip')

    def test_success(self, call, patched_zip):
        """Check successful unzipping of a Tensor file."""
        tensor_zip, unpack_path = call()
        mock_zip, mock_archive = patched_zip
        mock_zip.assert_called_once_with(tensor_zip, 'r')
        mock_archive.extractall.assert_called_once_with(unpack_path)

    def test_unzip_fails(self, call, patched_zip):
        """Check complaints if unpacking a Tensor fails."""
        _, mock_archive = patched_zip
        mock_archive.extractall.side_effect = OSError('Extraction failed')
        with pytest.raises(OSError, match='Extraction failed'):
            call()
