"""Tests for the optimizer chain of viperleed.calc.vlj.search.

Everything around the chain is mocked: these tests do not need the
viperleed_jax plugin.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-26'
__license__ = 'GPLv3+'

import sys
from types import ModuleType
from types import SimpleNamespace

import numpy as np
import pytest

from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.vlj import search


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
        def _run(chain, r_factor_type=RFactorType.PENDRY):
            monkeypatch.chdir(tmp_path)
            (tmp_path / 'Tensors/Tensors_001').mkdir(parents=True)
            rpars = mocker.MagicMock(TENSOR_INDEX=1, last_R=None,
                                     R_FACTOR_TYPE=r_factor_type,
                                     stored_R={}, files_to_out=set())
            slab = mocker.MagicMock(name='slab')
            calculator = mocker.MagicMock(name='calculator')
            _run.calculator = calculator
            calculator.R.side_effect = (
                # Integer and fractional R are the x they were called at
                lambda x, groups=None, **_: 0.5 if groups is None
                else tuple(x)
                )
            calculator.expand_params.return_value = (np.array([0.0]),)
            iterator = mocker.MagicMock()
            iterator.__iter__.return_value = iter(
                (SimpleNamespace(name=name), result)
                for name, result in chain
                )
            fake_plugin_modules = {
                'viperleed_jax.from_objects': {
                    'setup_tl_parameter_space': mocker.MagicMock(),
                    'setup_tl_calculator': mocker.MagicMock(
                        return_value=calculator
                        ),
                    },
                'viperleed_jax.optimization.iterator': {
                    'OptimizerIterator': mocker.MagicMock(
                        return_value=iterator
                        ),
                    },
                'viperleed_jax.utils': {
                    'benchmark_calculator': mocker.MagicMock(),
                    'format_benchmark_results': mocker.MagicMock(
                        return_value=''
                        ),
                    },
                }
            for mod_name, attrs in fake_plugin_modules.items():
                module = ModuleType(mod_name)
                module.__dict__.update(attrs)
                monkeypatch.setitem(sys.modules, mod_name, module)
            patches = {
                'VLJ_AVAILABLE': True,
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
        assert rpars.stored_R['superpos'] == (0.11, 0.1, 0.2)
        assert len(final_x) == 1
        np.testing.assert_array_equal(final_x[0], [0.1, 0.2])

    def test_last_better(self, run_chain):
        """Check that a later, better result replaces an earlier one."""
        chain = (('CMAES', self._result(0.17, [0.1, 0.2])),
                 ('SLSQP', self._result(0.11, [0.3, 0.4])))
        rpars, final_x = run_chain(chain)
        assert rpars.stored_R['superpos'] == (0.11, 0.3, 0.4)
        assert len(final_x) == 1
        np.testing.assert_array_equal(final_x[0], [0.3, 0.4])

    def test_three_first_best(self, run_chain):
        """Check that each result is compared with the best so far."""
        chain = (('CMAES', self._result(0.11, [0.1, 0.2])),
                 ('BFGS', self._result(0.17, [0.3, 0.4])),
                 ('SLSQP', self._result(0.13, [0.5, 0.6])))
        rpars, final_x = run_chain(chain)
        assert rpars.stored_R['superpos'] == (0.11, 0.1, 0.2)
        np.testing.assert_array_equal(final_x[0], [0.1, 0.2])

    def test_three_middle_best(self, run_chain):
        """Check that a best result in the middle of the chain is kept."""
        chain = (('CMAES', self._result(0.17, [0.1, 0.2])),
                 ('BFGS', self._result(0.11, [0.3, 0.4])),
                 ('SLSQP', self._result(0.13, [0.5, 0.6])))
        rpars, final_x = run_chain(chain)
        assert rpars.stored_R['superpos'] == (0.11, 0.3, 0.4)
        np.testing.assert_array_equal(final_x[0], [0.3, 0.4])

    def test_equal_keeps_earlier(self, run_chain):
        """Check that a later result with equal R does not replace it."""
        chain = (('CMAES', self._result(0.11, [0.1, 0.2])),
                 ('SLSQP', self._result(0.11, [0.3, 0.4])))
        rpars, final_x = run_chain(chain)
        assert rpars.stored_R['superpos'] == (0.11, 0.1, 0.2)
        np.testing.assert_array_equal(final_x[0], [0.1, 0.2])

    def test_first_nan(self, run_chain):
        """Check that a NaN R of an earlier result counts as the worst."""
        chain = (('BFGS', self._result(np.nan, [0.1, 0.2])),
                 ('SLSQP', self._result(0.11, [0.3, 0.4])))
        rpars, final_x = run_chain(chain)
        assert rpars.stored_R['superpos'] == (0.11, 0.3, 0.4)
        np.testing.assert_array_equal(final_x[0], [0.3, 0.4])

    def test_last_nan(self, run_chain):
        """Check that a later NaN R does not replace a finite one."""
        chain = (('CMAES', self._result(0.11, [0.1, 0.2])),
                 ('SLSQP', self._result(np.nan, [0.3, 0.4])))
        rpars, final_x = run_chain(chain)
        assert rpars.stored_R['superpos'] == (0.11, 0.1, 0.2)
        np.testing.assert_array_equal(final_x[0], [0.1, 0.2])

    @pytest.mark.parametrize('r_type', RFactorType)
    def test_rfactor_name_passed(self, run_chain, r_type):
        """Check that the calculator is given the name of the R factor."""
        chain = (('CMAES', self._result(0.11, [0.1, 0.2])),)
        run_chain(chain, r_factor_type=r_type)
        run_chain.calculator.set_rfactor.assert_called_once_with(r_type.name)
