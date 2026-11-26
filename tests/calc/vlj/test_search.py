"""Tests for module viperleed.calc.vlj.search."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-26'
__license__ = 'GPLv3+'

import logging
from types import SimpleNamespace
import subprocess

import numpy as np
import pytest

from viperleed.calc.files import vibrocc
from viperleed.calc.symmetry import enforceSymmetry
from viperleed.calc.symmetry import findSymmetry
from viperleed.calc.vlj import search

viperleed_jax = pytest.importorskip('viperleed_jax')

# pylint: disable=wrong-import-position,wrong-import-order
# The plugin is optional: import it only after the skip above
from viperleed.calc.files.new_displacements.file import DisplacementsFile
from viperleed_jax.from_objects import setup_tl_parameter_space
from viperleed_jax.transformation_tree import tree
# pylint: enable=wrong-import-position,wrong-import-order

_DISPLACEMENTS = '''\
== SEARCH z

= GEO_DELTA
Ag L(1-4)   z = -0.2 0.2 0.01
'''


@pytest.fixture(name='ag100_parameter_space')
def fixture_ag100_parameter_space(ag100, data_path, tmp_path):
    """Return a plugin ParameterSpace for Ag(100) with a search segment."""
    slab, rpars, *_ = ag100
    vibrocc.readVIBROCC(rpars, slab,
                        str(data_path / 'Ag(100)/mergeDisp/VIBROCC'))
    slab.full_update(rpars)
    findSymmetry(slab, rpars)
    enforceSymmetry(slab, rpars)
    parameter_space = setup_tl_parameter_space(slab, rpars)
    displacements = tmp_path / 'DISPLACEMENTS'
    displacements.write_text(_DISPLACEMENTS)
    disp_file = DisplacementsFile()
    disp_file.read(displacements)
    parameter_space.apply_search_segment(disp_file.next(0.9))
    return parameter_space


class TestExportParameterSpace:
    """Tests for the graphical export of the plugin's ParameterSpace."""

    def test_export_called(self, ag100_parameter_space, monkeypatch):
        """Check that the plugin's export writes the expected file."""
        written = []
        monkeypatch.setattr(tree.UniqueDotExporter, 'to_picture',
                            lambda _, filename: written.append(filename))
        search._export_parameter_space(ag100_parameter_space, 'z')
        assert written == ['parameter_space_z.pdf']

    _dot_errors = {
        'dot missing': FileNotFoundError('dot'),
        'dot failed': subprocess.CalledProcessError(1, 'dot'),
        }

    @pytest.mark.parametrize('error', _dot_errors.values(), ids=_dot_errors)
    def test_dot_failure_warns(self, error, ag100_parameter_space,
                               monkeypatch, caplog):
        """Check that a failure of the `dot` tool only logs a warning."""
        def _fail(*_):
            raise error
        monkeypatch.setattr(tree.UniqueDotExporter, 'to_picture', _fail)
        with caplog.at_level(logging.WARNING):
            search._export_parameter_space(ag100_parameter_space, 'z')
        assert 'Could not export parameter space tree view' in caplog.text

    def test_interface_mismatch_raises(self):
        """Check that a missing export method is not swallowed."""
        with pytest.raises(AttributeError):
            search._export_parameter_space(object(), 'z')


class TestOptimizerChainBestResult:
    """Check that vlj_search keeps the best result of the optimizer chain."""

    @staticmethod
    def _result(best_r, best_x):
        """Return a fake optimizer result."""
        return SimpleNamespace(best_R=best_r,
                               best_x=np.array(best_x),
                               write_to_file=lambda _: None)

    @pytest.fixture(name='run_chain')
    def fixture_run_chain(self, tmp_path, monkeypatch, mocker):
        """Return a function that runs vlj_search with a fake chain."""
        def _run(chain):
            monkeypatch.chdir(tmp_path)
            (tmp_path / 'Tensors/Tensors_001').mkdir(parents=True)
            rpars = mocker.MagicMock(TENSOR_INDEX=1, last_R=None,
                                     stored_R={}, files_to_out=set())
            slab = mocker.MagicMock(name='slab')
            calculator = mocker.MagicMock(name='calculator')
            calculator.R.side_effect = (
                lambda x, groups=None, **_: 0.5 if groups is None
                else (0.9, 0.8)
                )
            calculator.expand_params.return_value = (np.array([0.0]),)
            iterator = mocker.MagicMock()
            iterator.__iter__.return_value = iter(
                (SimpleNamespace(name=name), result)
                for name, result in chain
                )
            patches = {
                'VLJ_AVAILABLE': True,
                'setup_tl_parameter_space': mocker.MagicMock(),
                'setup_tl_calculator': mocker.MagicMock(
                    return_value=calculator
                    ),
                'OptimizerIterator': mocker.MagicMock(return_value=iterator),
                'benchmark_calculator': mocker.MagicMock(),
                'format_benchmark_results': mocker.MagicMock(return_value=''),
                'determine_integer_or_fractional': mocker.MagicMock(),
                'writeVIBROCC': mocker.MagicMock(),
                '_export_parameter_space': mocker.MagicMock(),
                }
            for name, value in patches.items():
                monkeypatch.setattr(search, name, value, raising=False)
            monkeypatch.setattr(search.poscar, 'write', mocker.MagicMock())
            monkeypatch.setattr(search.copy, 'deepcopy',
                                lambda _: mocker.MagicMock(name='tmp_slab'))
            search.vlj_search(slab, rpars)
            final_x = [call.args[2] for call in
                       calculator.apply_to_slab.call_args_list
                       if call.args[0] is slab]
            return rpars, final_x
        return _run

    def test_last_worse(self, run_chain):
        """Check that an earlier, better result is kept."""
        chain = (('CMAES', self._result(0.11, [0.1, 0.2])),
                 ('SLSQP', self._result(0.17, [0.3, 0.4])))
        rpars, final_x = run_chain(chain)
        assert rpars.stored_R['superpos'] == (0.11, 0.9, 0.8)
        assert len(final_x) == 1
        np.testing.assert_array_equal(final_x[0], [0.1, 0.2])

    def test_last_better(self, run_chain):
        """Check that a later, better result replaces an earlier one."""
        chain = (('CMAES', self._result(0.17, [0.1, 0.2])),
                 ('SLSQP', self._result(0.11, [0.3, 0.4])))
        rpars, final_x = run_chain(chain)
        assert rpars.stored_R['superpos'] == (0.11, 0.9, 0.8)
        assert len(final_x) == 1
        np.testing.assert_array_equal(final_x[0], [0.3, 0.4])
