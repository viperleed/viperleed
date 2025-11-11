"""Tests for module run_sections of viperleed.calc.sections."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-28'
__license__ = 'GPLv3+'

from collections import defaultdict

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.domain_params import DomainParameters
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.rparams.special.max_tl_displacement import (
    MaxTLAction,
    )
from viperleed.calc.sections.run_sections import run_section
from viperleed.calc.sections.run_sections import section_loop


_MODULE = 'viperleed.calc.sections.run_sections'


class TestRunSection:
    """Tests for the run_section function."""

    def test_search_marks_domains_as_require_refcalc(self, mocker):
        """Check that all domains require a refcalc after a search."""
        rpars = Rparams()
        for file in rpars.fileLoaded:
            # Fake that all files have been read already to skip
            # all the code that loads files before running.
            rpars.fileLoaded[file] = True
        rpars.domainParams = [mocker.MagicMock(refcalc_required=False)]
        mock_search = mocker.patch('viperleed.calc.sections.search.search')
        run_section(3, mocker.MagicMock(name='slab'), rpars)
        mock_search.assert_called_once()
        assert all(d.refcalc_required for d in rpars.domainParams)


class TestSectionLoop:
    """Tests for the section_loop function."""

    @fixture(name='patch_common', autouse=True)
    def patch_common(self, mocker):
        """Replace implementation details with mocks."""
        for func in ('viperleed.calc.files.parameters.update',
                     f'{_MODULE}.cleanup',
                     f'{_MODULE}.move_oldruns'):
            mocker.patch(func, return_value=None)

    @fixture(name='dummy_slab')
    def _make_dummy_slab(self, mocker):
        slab = mocker.MagicMock()
        slab.preprocessed = False
        slab.__iter__.return_value = []
        return slab

    @fixture(name='dummy_rp')
    def _make_dummy_rp(self):
        rpars = Rparams()
        rpars.HALTING = 3    # don't halt unless specified
        rpars.fileLoaded['EXPBEAMS'] = True
        return rpars

    def make_fake_run_section(self, actions_per_index=None):
        """Factory that returns a stateful fake_run_section.

        Parameters
        ----------
        actions_per_index: dict, optional
            Keys are section indices, values are sequences of callables
            representing actions to perform when the corresponding
            section is executed. Each action is a callable taking
            a single `rpars` argument.
        """
        call_counts = defaultdict(int)

        def fake_run_section(index, _slab, rpars):
            call_counts[index] += 1
            rpars.runHistory.append(index)

            if actions_per_index and index in actions_per_index:
                action_list = actions_per_index[index]
                call_index = call_counts[index] - 1  # 0-based
                if call_index < len(action_list):
                    action = action_list[call_index]
                    action(rpars)

        return fake_run_section

    def make_attr_setter_list(attr_name, values):
        """
        Returns a list of functions, each of which sets `rp.attr_name` to a
        specific value.

        Useful for use in `actions = {index: [func1, func2, ...]}` where each
        func is called per invocation.
        """
        return [lambda rp, val=val: setattr(rp, attr_name, val)
                for val in values]

    def make_is_too_far_by_section(which_calls):
        """Return a function that returns True only if it's called for the Nth
        time and N is in the which_calls argument."""
        call_count = {'count': -1}

        def is_too_far(_atom):
            # simulate one call per section
            call_count['count'] += 1
            # pylint: disable-next=unsupported-membership-test  # dict
            return call_count['count'] in which_calls
        return is_too_far

    def set_multi_search(rpars):
        rpars.disp_blocks = [('dummy_lines', 'search 1'),
                             ('dummy_lines', 'search 2')]

    def set_looped_search(rpars):
        rpars.disp_blocks = [('dummy_lines', 'search 1')]
        rpars.disp_loops = [(0, 0)]

    def set_nested_looped_search(rpars):
        """Inner loop at beginning"""
        rpars.disp_blocks = [('dummy_lines', 'search 1'),
                             ('dummy_lines', 'search 2'),]
        rpars.disp_loops = [(0, 0), (0, 1)]

    def set_nested_looped_search_2(rpars):
        """Inner loop at end"""
        rpars.disp_blocks = [('dummy_lines', 'search 1'),
                             ('dummy_lines', 'search 2'),]
        rpars.disp_loops = [(0, 1), (1, 1)]

    def raise_fake_error(_rpars):
        raise RuntimeError('Testing Error')

    def raise_keyboard_interrupt(_rpars):
        raise KeyboardInterrupt

    valid_runs = {  # name: (RUN, actions, expect)
        'init': ([0], {}, [0]),
        'halt after init': ([0, 1], {
            0: [lambda rp: setattr(rp, 'HALTING', 2)]
            }, [0]),
        'refcalc, with EXPBEAMS': ([0, 1], {}, [0, 1, 11]),
        'one search block': ([0, 1, 2, 3], {}, [0, 1, 11, 2, 3, 31, 12]),
        'explicit rfac': ([0, 1, 11, 2, 3, 31, 12], {},
                          [0, 1, 11, 2, 3, 31, 12]),
        '2 blocks': ([0, 1, 2, 3], {0: [set_multi_search]},
                     [0, 1, 11, 2, 3, 31, 12, 2, 3, 31, 12]),
        '2 blocks explicit': ([0, 1, 2, 3, 2, 3], {0: [set_multi_search]},
                              [0, 1, 11, 2, 3, 31, 12, 2, 3, 31, 12]),
        'loop': ([0, 1, 2, 3], {
                0: [set_looped_search],
                3: make_attr_setter_list('last_R', [0.2, 0.2])
                },
            [0, 1, 11] + [2, 3, 31, 12]*2),
        'loop with improvement': ([0, 1, 2, 3], {
            0: [set_looped_search],
            3: make_attr_setter_list('last_R', [0.3, 0.2, 0.2])
            },
            [0, 1, 11] + [2, 3, 31, 12]*3),
        'nested loop trivial': ([0, 1, 2, 3], {
            0: [set_nested_looped_search],
            3: make_attr_setter_list('last_R', [0.2]*5)
            },
            [0, 1, 11] + [2, 3, 31, 12]*5),
        'nested loop improvement': ([0, 1, 2, 3], {
            0: [set_nested_looped_search],
            3: make_attr_setter_list(
                'last_R', [0.3, 0.25, 0.25] + [0.24]*4
                ),
            },
            [0, 1, 11] + [2, 3, 31, 12]*7),
        'nested loop trivial 2': ([0, 1, 2, 3], {
            0: [set_nested_looped_search_2],
            3: make_attr_setter_list('last_R', [0.2]*5)
            },
            [0, 1, 11] + [2, 3, 31, 12]*5),
        'nested loop improvement 2': ([0, 1, 2, 3], {
            0: [set_nested_looped_search_2],
            3: make_attr_setter_list(
                'last_R', [0.3, 0.25] + [0.24]*4
                ),
            },
            [0, 1, 11] + [2, 3, 31, 12]*6),
        }

    domain_runs = {  # name: (RUN, actions, expect)
        'init': ([0], {}, [0]),
        'refcalc': ([0, 1], {}, [0, 1]),
        'one search block': ([0, 1, 2, 3], {}, [0, 1, 2, 3, 31, 12]),
        'explicit rfac': ([0, 1, 2, 3, 31, 12], {},
                          [0, 1, 2, 3, 31, 12]),
        '2 blocks': ([0, 1, 2, 3], {0: [set_multi_search]},
                     [0, 1, 2, 3, 31, 12, 2, 3, 31, 12]),
        'loop': ([0, 1, 2, 3], {
            0: [set_looped_search],
            3: make_attr_setter_list('last_R', [0.2, 0.2])
            },
            [0, 1] + [2, 3, 31, 12]*2),
        'loop with improvement': ([0, 1, 2, 3], {
            0: [set_looped_search],
            3: make_attr_setter_list('last_R', [0.3, 0.2, 0.2])
            },
            [0, 1] + [2, 3, 31, 12]*3),
        }

    max_disp_runs = {   # name: (RUN, actions, which, expect_code, expect)
        'default': ([0, 1, 2, 3], {0: [set_multi_search]}, MaxTLAction.REFCALC,
                    [0, 1, 11, 2, 3, 31, 12, 1, 11, 2, 3, 31, 12]),
        'ignore': ([0, 1, 2, 3], {0: [set_multi_search]}, MaxTLAction.IGNORE,
                   [0, 1, 11, 2, 3, 31, 12, 2, 3, 31, 12]),
        'stop': ([0, 1, 2, 3], {0: [set_multi_search]}, MaxTLAction.STOP,
                 [0, 1, 11, 2, 3, 31, 12]),
        'timeout': ([0, 1, 2, 3], {
            0: [set_multi_search],
            1: [lambda rp: setattr(rp, 'last_refcalc_time', 1e5)]
                }, MaxTLAction.REFCALC, [0, 1, 11, 2, 3, 31, 12]),
        'known time': ([0, 1, 2, 3], {
            0: [set_multi_search],
            1: [lambda rp: setattr(rp, 'last_refcalc_time', 1)]
                }, MaxTLAction.REFCALC, [0, 1, 11, 2, 3, 31, 12,
                                         1, 11, 2, 3, 31, 12]),
        'loop with ignore': ([0, 1, 2, 3], {
            0: [set_looped_search],
            3: make_attr_setter_list('last_R', [0.2, 0.2])
            },
            MaxTLAction.IGNORE,
            [0, 1, 11] + [2, 3, 31, 12]*2),
        }

    stopped_runs = {  # name: (RUN, actions, expect_code, expect_history)
        'user STOP': ([0, 1, 2, 3], {
            1: [lambda rp: setattr(rp, 'STOP', True)]}, 0, [0, 1, 11]),
        'halting': ([0, 1, 2, 3], {
            1: [lambda rp: setattr(rp, 'halt', 3)]}, 0, [0, 1]),
        'KeyboardInterrupt': ([0, 1, 2, 3], {
            1: [raise_keyboard_interrupt]}, 1, [0, 1]),
        'Exception': ([0, 1, 2, 3], {
            1: [raise_fake_error]}, 3, [0, 1]),
        }

    @parametrize('run,actions,expect', valid_runs.values(), ids=valid_runs)
    def test_section_loop_valid(self, run, actions, expect,
                                dummy_rp, dummy_slab, monkeypatch):
        dummy_rp.RUN = run[:]

        monkeypatch.setattr(f'{_MODULE}.run_section',
                            self.make_fake_run_section(actions))
        exit_code, _ = section_loop(dummy_rp, dummy_slab)

        assert exit_code == 0
        assert dummy_rp.runHistory == expect

    @parametrize('run,actions,expect', domain_runs.values(), ids=domain_runs)
    def test_section_loop_valid_with_domains(self, run, actions, expect,
                                             dummy_rp, dummy_slab,
                                             monkeypatch):
        dummy_rp.RUN = run[:]
        domain = DomainParameters('.', 'dummy_name')
        domain.rpars = Rparams()
        domain.rpars.HALTING = 3
        domain.slab = dummy_slab
        dummy_rp.domainParams = [domain]

        monkeypatch.setattr(f'{_MODULE}.run_section',
                            self.make_fake_run_section(actions))
        exit_code, _ = section_loop(dummy_rp, dummy_slab)

        assert exit_code == 0
        assert dummy_rp.runHistory == expect

    @pytest.mark.timeout(5)
    @parametrize('run,actions,which,expect', max_disp_runs.values(),
                 ids=max_disp_runs)
    def test_section_loop_valid_max_tl_disp(self, run, actions, which, expect,
                                            dummy_rp, dummy_slab, monkeypatch):
        dummy_rp.RUN = run[:]
        dummy_rp.MAX_TL_DISPLACEMENT.action = which

        # always exceed MAX_TL_DISPLACEMENT limit
        monkeypatch.setattr(
            f'{_MODULE}._check_exceeds_tl_limit',
            lambda *_args, **_kwargs: True)
        monkeypatch.setattr(f'{_MODULE}.run_section',
                            self.make_fake_run_section(actions))
        exit_code, _ = section_loop(dummy_rp, dummy_slab)

        assert exit_code == 0
        assert dummy_rp.runHistory == expect

    @parametrize('run,actions,expect_code,expect_history',
                 stopped_runs.values(), ids=stopped_runs)
    def test_section_loop_stops(self, run, actions,
                                expect_code, expect_history, dummy_rp,
                                dummy_slab, monkeypatch):
        dummy_rp.RUN = run[:]

        monkeypatch.setattr(f'{_MODULE}.run_section',
                            self.make_fake_run_section(actions))
        exit_code, _ = section_loop(dummy_rp, dummy_slab)

        assert exit_code == expect_code
        assert dummy_rp.runHistory == expect_history
