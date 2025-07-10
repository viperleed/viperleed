"""Tests for module deltas of viperleed.calc.sections.

This module collects tests for various helper functions
in the deltas module.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-07-08'
__license__ = 'GPLv3+'

from pathlib import Path
import logging

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.deltas import DeltasError
from viperleed.calc.sections.deltas import _ensure_tensors_loaded
from viperleed.calc.sections.deltas import _prepare_log_file
from viperleed.calc.sections.deltas import _remove_old_param_file
from viperleed.calc.sections.deltas import _sort_current_deltas_by_element


class TestEnsureTensorsLoaded:
    """Tests for the _ensure_tensors_loaded helper function."""

    @fixture(name='mocks')
    def fixture_mocks(self, mocker):
        """Replace implementation details with mocks."""
        return {
            'is_dir': mocker.patch('pathlib.Path.is_dir', return_value=True),
            'fetch tensor': mocker.patch(
                'viperleed.calc.files.iotensors.fetch_unpacked_tensor'
                ),
            'load refcalc': mocker.patch(
                'viperleed.calc.files.iotensors.getTensorOriStates',
                ),
            }

    def test_loads_refacalc_state(self, rpars, mocks, mocker, caplog):
        """Check loading of the most recent refcalc state into the slab."""
        caplog.set_level(logging.DEBUG)
        mock_slab = mocker.MagicMock()
        rpars.runHistory = [0, 2, 3, 12, 99, 57]  # No refcalc (== 1)
        _ensure_tensors_loaded(mock_slab, rpars)

        tensor_path = (
            Path(DEFAULT_TENSORS)
            / f'{DEFAULT_TENSORS}_{rpars.TENSOR_INDEX}'
            )
        mocks['fetch tensor'].assert_called_once_with(rpars.TENSOR_INDEX)
        mocks['load refcalc'].assert_called_once_with(mock_slab, tensor_path)
        mock_slab.restoreOriState.assert_called_once_with(keepDisp=True)
        expect_log = 'Running without reference calculation'
        assert expect_log in caplog.text

    def test_raises_without_tensors_folder(self, mocks, caplog):
        """Check complaints when deltas is called without a Tensors folder."""
        mocks['is_dir'].return_value = False
        with pytest.raises(FileNotFoundError):
            _ensure_tensors_loaded('mock_slab', 'mock_rpars')
        expect_log = 'No Tensors directory found.'
        assert expect_log in caplog.text

    def test_refcalc_already_done(self, rpars, mocks, caplog):
        """Check that tensors are not reloaded when a refcalc was done."""
        caplog.set_level(0)   # All messages
        rpars.runHistory.append(1)
        _ensure_tensors_loaded('mock_slab', rpars)
        assert not caplog.text
        mocks['fetch tensor'].assert_called_once_with(rpars.TENSOR_INDEX)
        mocks['load refcalc'].assert_not_called()



class TestPrepareLogFile:
    """Tests for the _prepare_log_file helper."""

    def check_file_contents(self, file_path):
        """Ensure the delta log file has the expected contents."""
        log_header = ('Logs from multiple delta calculations are collected '
                      'here. Their order may not be preserved.\n')
        contents = file_path.read_text(encoding='utf-8')
        assert contents == log_header

    @fixture(name='mock_rpars')
    def fixture_rpars(self, mocker):
        """Return a fake Rparams."""
        return mocker.MagicMock(timestamp='20240625')

    @fixture(name='call')
    def fixture_call(self, mock_rpars, tmp_path):
        """Call _prepare_log_file in a temporary directory."""
        def _call(**kwargs):
            with execute_in_dir(tmp_path):
                result = _prepare_log_file(mock_rpars, **kwargs)
            assert result == f'delta-{mock_rpars.timestamp}.log'
            return tmp_path/result
        return _call

    @parametrize(subdomain=(True, False))
    def test_sucess(self, subdomain, call, caplog):
        """Check successful file creation and logging messages."""
        caplog.set_level(logging.INFO)
        log_file = call(subdomain=subdomain)
        self.check_file_contents(log_file)
        logging_msgs = (
            'Generating delta files...',
            'collected in delta-20240625.log',
            )
        if subdomain:
            assert not any(msg in caplog.text for msg in logging_msgs)
        else:
            assert all(msg in caplog.text for msg in logging_msgs)

    def test_oserror_warns(self, call, mocker, caplog):
        """Check that failure to create the log file does not raise."""
        caplog.set_level(logging.WARNING)
        mocker.patch('pathlib.Path.write_text', side_effect=OSError)
        call(subdomain=False)
        expect_log = (
            'Error creating delta log file',
            'will not affect execution',
            )
        assert all(msg in caplog.text for msg in expect_log)


class TestRemoveOldParamFile:
    """Tests for the _remove_old_param_file helper function."""

    @fixture(name='call')
    def fixture_call(self, tmp_path):
        """Call _remove_old_param_file in a temporary directory."""
        def _call():
            old_param = tmp_path/'PARAM'
            renamed = old_param.with_name('PARAM-old')
            contents = 'test'
            old_param.write_text(contents)
            with execute_in_dir(tmp_path):
                _remove_old_param_file()
            return old_param, renamed, contents
        return _call

    def test_deletes_param_if_rename_fails(self, call, mocker):
        """Check deletion of an existing PARAM when renaming fails."""
        mocker.patch('pathlib.Path.replace', side_effect=OSError)
        param, renamed, _ = call()
        assert not param.exists()
        assert not renamed.exists()

    def test_does_nothing_if_no_param_file(self, tmp_path):
        """Check that no action is take if PARAM does not exist."""
        with execute_in_dir(tmp_path):
            _remove_old_param_file()
        assert not (tmp_path / 'PARAM').exists()
        assert not (tmp_path / 'PARAM-old').exists()

    def test_renames_existing_param_file(self, call):
        """Check successful renaming of a PARAM file."""
        param, renamed, contents = call()
        assert not param.exists()
        assert renamed.exists()
        assert renamed.read_text() == contents

    def test_warns_if_remove_fails(self, call, mocker, caplog):
        """Check warnings are emitted if removal of PARAM fails."""
        mocker.patch('pathlib.Path.replace', side_effect=OSError)
        mocker.patch('pathlib.Path.unlink', side_effect=OSError)
        param, *_ = call()
        expect_log = 'Cannot rename/remove old PARAM file'
        assert expect_log in caplog.text
        assert param.exists()


class TestSortCurrentDeltasByElement:
    """Tests for the _sort_current_deltas_by_element helper."""

    @fixture(name='mock_atoms')
    def fixture_mock_atoms(self, mocker):
        """Return a few fake Atoms."""
        fake_atom_1 = mocker.MagicMock(
            disp_occ=dict.fromkeys(('A', 'B', '00c')),
            current_deltas=[
                'DEL_1_00C_02',
                'DEL_1_B_01',
                'DEL_1_a_01',
                ],
            )
        fake_atom_2 = mocker.MagicMock(
            disp_occ=dict.fromkeys(('O', 'fe')),
            current_deltas=[
                'DEL_2_vac_1',
                'DEL_2_O_1',
                'DEL_2_Fe_2',
                ],
            )
        return [fake_atom_1, fake_atom_2], [fake_atom_2]

    def test_sucess(self, mock_atoms):
        """Check the expected result of sorting atoms."""
        (atom_1, atom_2), _ = mock_atoms
        _sort_current_deltas_by_element(*mock_atoms)
        assert atom_1.current_deltas == [
            'DEL_1_a_01',
            'DEL_1_B_01',
            'DEL_1_00C_02',
            ]
        assert atom_2.current_deltas == [
            'DEL_2_O_1',
            'DEL_2_Fe_2',
            'DEL_2_vac_1',
            ]

    def test_raises_without_delta(self, mock_atoms):
        """Ensure complaints when an element has no known delta."""
        (failing_atom, _), _ = mock_atoms
        failing_atom.disp_occ['NoDeltasForThis'] = None
        with pytest.raises(DeltasError, match='Found no delta files'):
            _sort_current_deltas_by_element(*mock_atoms)

    def test_raises_without_vacancy_delta(self, mock_atoms):
        """Ensure complaints when a vacancy has no known delta."""
        (_, failing_atom), _ = mock_atoms
        failing_atom.current_deltas.remove('DEL_2_vac_1')
        with pytest.raises(DeltasError, match='Found no delta files'):
            _sort_current_deltas_by_element(*mock_atoms)

    def test_raises_multiple_deltas(self, mock_atoms):
        """Ensure complaints when there are multiple deltas for one element."""
        (failing_atom, _), _ = mock_atoms
        failing_atom.current_deltas.append('DEL_1_a_another')
        with pytest.raises(DeltasError, match='Found multiple delta files'):
            _sort_current_deltas_by_element(*mock_atoms)

    def test_raises_multiple_vacancy_deltas(self, mock_atoms):
        """Ensure complaints when there are multiple deltas for a vacancy."""
        (_, failing_atom), _ = mock_atoms
        failing_atom.current_deltas.append('DEL_2_vac_another')
        with pytest.raises(DeltasError, match='Found multiple delta files'):
            _sort_current_deltas_by_element(*mock_atoms)
