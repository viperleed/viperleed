"""Tests for cleanup.organize_workdir (+ helpers) of viperleed.calc.section."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-03'
__license__ = 'GPLv3+'

import logging
import shutil
from zipfile import ZipFile

from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import LOG_PREFIX
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.cleanup import _OUT_FILES
from viperleed.calc.sections.cleanup import _SUPP_DIRS
from viperleed.calc.sections.cleanup import _SUPP_FILES
from viperleed.calc.sections.cleanup import organize_workdir
from viperleed.calc.sections.cleanup import _collect_delta_files
from viperleed.calc.sections.cleanup import _collect_out_contents
from viperleed.calc.sections.cleanup import _collect_supp_contents
from viperleed.calc.sections.cleanup import _copy_files_and_directories
from viperleed.calc.sections.cleanup import _zip_folder
from viperleed.calc.sections.cleanup import _zip_subfolders

from ....helpers import filesystem_from_dict
from ....helpers import filesystem_to_dict
from .conftest import _MODULE


class TestCollectDeltaFiles:
    """Tests for the _collect_delta_files function."""

    @fixture(name='run')
    def fixture_run(self, workdir):
        """Execute _collect_delta_files in workdir."""
        def _run(tensor_index, tree=None):
            if tree:
                filesystem_from_dict(tree, workdir)
            with execute_in_dir(workdir):
                _collect_delta_files(tensor_index)
            return filesystem_to_dict(workdir)
        return _run

    def test_delta_folder_exists(self, run, workdir, caplog):
        """Check no complaints if the target folder is there already."""
        (workdir/'Deltas'/'Deltas_123').mkdir(parents=True)
        self.test_success(run, caplog)

    def test_file_copy_fails(self, run, caplog, mocker):
        """Check complaints when moving files fails."""
        mock_move = mocker.patch('pathlib.Path.replace', side_effect=OSError)
        tree = {'DEL_1': '', 'DEL_2': ''}
        collected = run(75, tree)
        expect = {'Deltas': {'Deltas_075': {}},  # Fails to move all
                  **tree}
        expect_log = 'Error moving Delta'
        assert collected == expect
        assert mock_move.call_count == len(tree)
        assert expect_log in caplog.text

    def test_mkdir_fails(self, run, caplog, mocker):
        """Check warnings when failing to create the Deltas directory."""
        mocker.patch('pathlib.Path.mkdir', side_effect=OSError)
        tree = {'DEL_1': 'contents'}
        unchanged = run(5, tree)
        expect_log = 'Failed to create'
        assert unchanged == tree  # No changes
        assert expect_log in caplog.text

    def test_nothing_to_collect(self, run, workdir, caplog):
        """Check that no folder is created if no files should be collected."""
        run(9)
        assert not (workdir/'Deltas'/'Deltas_009').exists()
        assert not caplog.text

    def test_success(self, run, caplog):
        """Check results of a successful execution."""
        tree = {f: f'{f} contents'
                for f in ('DEL_1', 'DEL_two', 'DEL_tafile')}
        collected = run(123, tree)
        expect = {'Deltas': {'Deltas_123': tree}}
        assert collected == expect
        assert not caplog.text


class TestCollectSuppOut:
    """Tests for _collect_out_contents and _collect_supp_contents."""

    @property
    def generated(self):
        """Files that were generated/edited during a run, they go to OUT."""
        return {f: f'{f} contents' for f in (
            'PARAMETERS',
            'POSCAR',
            'POSCAR_mincell',  # Name overlap with SUPP
            'VIBROCC',
            )}

    @property
    def to_out(self):
        """Files that should be stored in OUT."""
        to_out = {f: f'{f} contents' for f in _OUT_FILES}
        to_out.update(self.generated)
        to_out.update({f'R_{s}_OUT_R=0.123': f'R factor {s} contents'
                       for s in ('refcalc', 'superpos', 'some-other-section')})
        return to_out

    @property
    def to_supp(self):
        """Files that should be stored in SUPP."""
        generated = self.generated
        to_supp = {f: f'{f} contents'
                   for f in _SUPP_FILES
                   if f not in generated}
        to_supp.update({d: {'test': f'test file in {d}'} for d in _SUPP_DIRS})
        to_supp[f'not-a-{LOG_PREFIX}-log.log'] = 'some log contents'
        return to_supp

    @property
    def tree(self):
        """Return a dict of files/folders in the main work tree."""
        # Stuff that is untouched
        stays = {
            f'{LOG_PREFIX}-main-log.log': 'fake calc log contents',
            'a-log-for-a-compile-step.log': 'fake compilation log',
            'another_folder': {'with-contents.txt': ''},
            }
        return {**stays, **self.to_out, **self.to_supp}

    @fixture(name='work_tree')
    def fixture_work_tree(self, workdir):
        """Create files to be moved to SUPP/OUT at workdir."""
        filesystem_from_dict(self.tree, workdir)
        return workdir

    @fixture(name='run')
    def fixture_run(self, rpars, work_tree):
        """Execute _organize_supp_out."""
        rpars.files_to_out.update(self.generated)
        def _run(*which):
            expect = {**self.tree,
                      DEFAULT_SUPP: self.to_supp,
                      DEFAULT_OUT: self.to_out}
            with execute_in_dir(work_tree):
                for func in which:
                    func(rpars)
            organized = filesystem_to_dict(work_tree)
            return expect, organized
        return _run

    _ordered = {
        'out,supp': (_collect_out_contents, _collect_supp_contents),
        'supp,out': (_collect_supp_contents, _collect_out_contents),
        }

    @parametrize(which=_ordered.values(), ids=_ordered)
    def test_success(self, which, run, caplog):
        """Check a successful copy of contents to SUPP/OUT."""
        expect, organized = run(*which)
        assert organized == expect
        assert not caplog.text


class TestCopyFilesAndDirectories:
    """Tests for the _copy_files_and_directories helper."""

    tree = {
        'test.txt': 'test file contents',
        'test_directory': {'test_subfolder': {'test_file': 'another file'}},
        }

    @fixture(name='run')
    def fixture_run(self, workdir):
        """Run _copy_files_and_directories in workdir."""
        def _run(*args, make_tree=True):
            if make_tree:
                filesystem_from_dict(self.tree, workdir)
            with execute_in_dir(workdir):
                _copy_files_and_directories(*args)
            return filesystem_to_dict(workdir)
        return _run

    def test_copy_directories(self, run, workdir):
        """Check a successful copy of directories."""
        target = workdir / 'test_target'
        target.mkdir()
        dirs = [workdir/d for d, c in self.tree.items() if isinstance(c, dict)]
        copied = run([], dirs, target)
        expect = self.tree.copy()
        expect[target.name] = {d.name: self.tree[d.name] for d in dirs}
        assert copied == expect

    def test_copy_fails(self, run, workdir, caplog, mocker):
        """Check logging messages when copying fails."""
        copy_file = mocker.patch('shutil.copy2', side_effect=OSError)
        copy_dir = mocker.patch(f'{_MODULE}.copytree_exists_ok',
                                side_effect=OSError)
        expect = self.tree.copy()
        target = workdir / 'test_target'
        expect[target.name] = {}
        files = [workdir/f for f, c in self.tree.items() if isinstance(c, str)]
        dirs = [workdir/d for d, c in self.tree.items() if isinstance(c, dict)]
        copied = run(files, dirs, target)
        assert copied == expect
        assert copy_file.call_count == len(files)
        assert copy_dir.call_count == len(dirs)
        log_msg_file = 'Error moving test_target file {}'
        log_msg_dir = 'Error moving test_target directory {}'
        log = caplog.text
        assert all(log_msg_file.format(f.name) in log for f in files)
        assert all(log_msg_dir.format(d.name) in log for d in dirs)

    def test_copy_files(self, run, workdir):
        """Check a successful copy of files."""
        target = workdir / 'test_target'
        files = [workdir/f for f, c in self.tree.items() if isinstance(c, str)]
        copied = run(files, [], target)
        expect = self.tree.copy()
        expect[target.name] = {f.name: self.tree[f.name] for f in files}
        assert copied == expect

    def test_copy_non_existing(self, run, workdir):
        """Check no complaints when copying non-existing files/folders."""
        target = workdir / 'test_target'
        expect = {target.name: {}}
        files = (workdir/'non_existing_file.txt',)
        dirs = (workdir/'non_existing_dir',)
        copied = run(files, dirs, target, make_tree=False)
        assert copied == expect

    def test_mkdir_fails(self, run, workdir, caplog, mocker):
        """Check warnings when making the target fails."""
        files = [workdir/f for f, c in self.tree.items() if isinstance(c, str)]
        target = workdir / 'test_target'
        filesystem_from_dict(self.tree, workdir)  # BEFORE patching
        mocker.patch('pathlib.Path.mkdir', side_effect=OSError)
        expect_log = 'Error creating test_target folder'
        unchanged = run(files, (), target, make_tree=False)
        assert unchanged == self.tree
        assert expect_log in caplog.text


class TestOrganizeWorkdir:
    """Tests for the organize_workdir function."""

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details of organize_workdir with mocks."""
        helpers = (
            '_collect_delta_files',
            '_zip_subfolders',
            '_collect_supp_contents',
            '_collect_out_contents',
            )
        return {helper: mocker.patch(f'{_MODULE}.{helper}')
                for helper in helpers}

    @fixture(name='run')
    def fixture_run(self, rpars, workdir):
        """Execute organize_workdir with rpars and path=workdir."""
        def _run(**kwargs):
            organize_workdir(rpars, workdir, **kwargs)
        return _run

    def test_archive_and_delete(self, rpars, run, mock_implementation, mocker):
        """Check that the implementation is executed as expected."""
        rpars.TENSOR_INDEX = mocker.MagicMock()
        rpars.ZIP_COMPRESSION_LEVEL = level = mocker.MagicMock()
        run(delete_unzipped=False, tensors=True, deltas=True)
        calls = {
            '_collect_delta_files': (mocker.call(rpars.TENSOR_INDEX),),
            '_zip_subfolders': (
                mocker.call(DEFAULT_TENSORS, True, False, level),
                mocker.call(DEFAULT_DELTAS, True, False, level),
                ),
            '_collect_supp_contents': (mocker.call(rpars),),
            '_collect_out_contents': (mocker.call(rpars),),
            }
        for helper, mock in mock_implementation.items():
            mock.assert_has_calls(calls[helper])

    def test_archive_deltas(self, run, rpars, mock_implementation, mocker):
        """Check that the implementation is executed as expected."""
        run(delete_unzipped=False, tensors=False, deltas=True)
        calls = {
            '_collect_delta_files': (
                mocker.call(0),  # TENSOR_INDEX is None by default
                ),
            '_zip_subfolders': (
                mocker.call(DEFAULT_DELTAS, True, False, 2),
                ),
            '_collect_supp_contents': (mocker.call(rpars),),
            '_collect_out_contents': (mocker.call(rpars),),
            }
        for helper, mock in mock_implementation.items():
            mock.assert_has_calls(calls[helper])

    def test_archive_tensors(self, run, rpars, mock_implementation, mocker):
        """Check that the implementation is executed as expected."""
        run(delete_unzipped=False, tensors=True, deltas=False)
        calls = {
            '_collect_delta_files': (
                mocker.call(0),  # TENSOR_INDEX is None by default
                ),
            '_zip_subfolders': (
                mocker.call(DEFAULT_TENSORS, True, False, 2),
                ),
            '_collect_supp_contents': (mocker.call(rpars),),
            '_collect_out_contents': (mocker.call(rpars),),
            }
        for helper, mock in mock_implementation.items():
            mock.assert_has_calls(calls[helper])

    def test_delete_only(self, run, rpars, mock_implementation, mocker):
        """Check that the implementation is executed as expected."""
        run(delete_unzipped=True, tensors=False, deltas=False)
        calls = {
            '_collect_delta_files': (
                mocker.call(0),  # TENSOR_INDEX is None by default
                ),
            '_zip_subfolders': (
                mocker.call(DEFAULT_TENSORS, False, True, 2),
                mocker.call(DEFAULT_DELTAS, False, True, 2),
                ),
            '_collect_supp_contents': (mocker.call(rpars),),
            '_collect_out_contents': (mocker.call(rpars),),
            }
        for helper, mock in mock_implementation.items():
            mock.assert_has_calls(calls[helper])


class TestZipFolder:
    """Tests for the _zip_folder function."""

    def test_does_not_pack_itself(self, tmp_path):
        """Check that an existing archive is not included."""
        folder = 'tst_folder'
        contents = {'file': 'file contents'}
        filesystem_from_dict({folder: contents}, tmp_path)
        # Create also an existing zip file
        archive = (tmp_path/folder).with_suffix('.zip')
        with ZipFile(archive, 'w') as zip_:
            zip_.writestr('pre-existing', 'contents were already there')

        _zip_folder(tmp_path/folder, 2)

        # Check that the archive contains what it should
        extract_into = tmp_path/'extracted'
        shutil.unpack_archive(archive, extract_into)
        extracted = filesystem_to_dict(extract_into)
        expect = {**contents,
                  'pre-existing': 'contents were already there'}
        assert extracted == expect


class TestZipSubfolders:
    """Tests for the _zip_subfolders function."""

    folder = 'test_zipping'
    contents = {
        'test_zipping_001': {'file_1': 'file_1 contents'},
        'test_zipping_9991': {'file_9991': 'file_9991 contents'},
        'not-a-directory': 'some stray file in test_zipping',
        }
    tree = {folder: contents}

    @property
    def packed_all(self):
        """Return the expected tree after packing self.tree."""
        files = {f: c for f, c in self.contents.items() if isinstance(c, str)}
        files.update({f'{f}.zip': f'{f}.zip is a binary file'
                      for f in self.contents
                      if f not in files})
        return {self.folder: files}

    @fixture(name='work_tree')
    def fixture_work_tree(self, workdir):
        """Create a temporary tree at workdir."""
        filesystem_from_dict(self.tree, workdir)
        return workdir

    @fixture(name='run')
    def fixture_run(self, work_tree):
        """Execute _zip_subfolders in workdir."""
        def _run(**kwargs):
            with execute_in_dir(work_tree):
                _zip_subfolders(at_path=self.folder, **kwargs)
            return filesystem_to_dict(work_tree)
        return _run

    def test_delete_fails(self, run, caplog, mocker):
        """Check complaints when deleting an unzipped folder fails."""
        mocker.patch('shutil.rmtree', side_effect=OSError)
        clean = run(delete_unzipped=True,
                    archive=True,
                    compression_level=2)
        expect = self.packed_all
        expect[self.folder].update(self.contents)
        expect_log = 'Error deleting unzipped'
        assert clean == expect
        assert expect_log in caplog.text

    def test_delete_only(self, run, caplog):
        """Check no deletion of unzipped directories."""
        clean = run(delete_unzipped=True,
                    archive=False,
                    compression_level=2)
        expect = {self.folder: {
            f: c for f, c in self.contents.items() if isinstance(c, str)
            }}
        assert clean == expect
        assert not caplog.text

    def test_do_nothing(self, run, caplog):
        """Check correct behavior for not packing nor deleting."""
        unchanged = run(delete_unzipped=False,
                        archive=False,
                        compression_level=2)
        assert unchanged == self.tree
        assert not caplog.text

    def test_folder_does_not_exist(self, work_tree, caplog):
        """Check that nothing happens if called with a non-existing folder."""
        _zip_subfolders('does-not-exist', True, True, 2)
        unchanged = filesystem_to_dict(work_tree)
        assert unchanged == self.tree
        assert not caplog.text

    def test_invalid_name_format(self, run, workdir, caplog):
        """Check that a folder with invalid format is not processed."""
        # Increase log level: There are INFO for the valid folders
        caplog.set_level(logging.WARNING)
        tree = {self.folder: {'test_zipping_001-invalid-fmt': {}}}
        filesystem_from_dict(tree, workdir)
        clean = run(delete_unzipped=True,
                    archive=True,
                    compression_level=2)
        expect = self.packed_all
        # The non-matching folder should sty where it is
        expect[self.folder].update(tree[self.folder])
        assert clean == expect
        assert not caplog.text

    def test_not_a_valid_folder(self, run, workdir, caplog):
        """Check that a stray subfolder is not packed."""
        # Increase log level: There are INFO for the valid folders
        caplog.set_level(logging.WARNING)
        tree = {self.folder: {'not-a-test_zipping_001': {}}}
        filesystem_from_dict(tree, workdir)
        clean = run(delete_unzipped=True,
                    archive=True,
                    compression_level=2)
        expect = self.packed_all
        expect[self.folder].update(tree[self.folder])
        assert clean == expect
        assert not caplog.text

    def test_pack_and_delete(self, run, workdir, caplog):
        """Check packing and deleting of both Tensors and Deltas."""
        # Increase log level: There are INFO messages about packing
        caplog.set_level(logging.WARNING)
        clean = run(delete_unzipped=True,
                    archive=True,
                    compression_level=2)
        assert clean == self.packed_all
        assert not caplog.text

        # Check the contents of an archive
        archived = next(f for f, c in self.contents.items()
                        if isinstance(c, dict))
        target = workdir/'extracted'
        target.mkdir()
        archive_p = (workdir/self.folder/archived).with_suffix('.zip')
        shutil.unpack_archive(archive_p, target)
        extracted = filesystem_to_dict(target)
        assert extracted == self.contents[archived]

    def test_pack_fails(self, run, caplog, mocker):
        """Check nothing is deleted if packing is not successful."""
        mocker.patch(f'{_MODULE}.ZipFile', side_effect=OSError)
        unchanged = run(delete_unzipped=True,  # not honored
                        archive=True,
                        compression_level=2)
        expect_log = 'Error packing'
        assert unchanged == self.tree
        assert expect_log in caplog.text

    def test_pack_only(self, run, caplog):
        """Check no deletion of unzipped directories."""
        # Increase log level: There are INFO messages about packing
        caplog.set_level(logging.WARNING)
        clean = run(delete_unzipped=False,
                    archive=True,
                    compression_level=2)
        expect = self.packed_all
        expect[self.folder].update(self.contents)
        assert clean == expect
        assert not caplog.text
