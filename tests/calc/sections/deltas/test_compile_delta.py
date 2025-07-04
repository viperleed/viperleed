"""Tests for module deltas of viperleed.calc.sections.

This module collects tests for the compile_delta function.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

from viperleed.calc.sections.deltas import compile_delta

from ..test_refcalc import TestCompileRefcalc


# Notice that the compile_delta and compile_refcalc functions are
# virtually identical, except for a few error messages and a few
# variable names. There is no need to add more tests. When doing
# #43, the tests can be given to the base-class method!
class TestCompileDelta(TestCompileRefcalc):
    """Tests for the compile_delta function."""

    compile_func = compile_delta
    compiler_cls_name = 'DeltaCompileTask'
    default_sources = 'delta.f', 'lib.tleed.f', 'lib.delta.f', 'GLOBAL'
    section_name = 'delta-amplitudes'
