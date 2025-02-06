"""Tests for viperleed.calc section deltas."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

from viperleed.calc.sections.deltas import compile_delta

from .test_refcalc import TestCompileRefcalc


class TestDeltasAg100:
    """Test the successful outcome of a delta-amplitude calculation."""

    def test_successful_run(self, delta_files_ag100):
        """Check that delta-amplitude calculation exits without errors."""
        assert not delta_files_ag100.failed
        assert delta_files_ag100.records is not None
        assert delta_files_ag100.records.get_last_state_for_section('delta')

    def test_delta_input_written(self, delta_files_ag100):
        """Check that an input file was correctly written."""
        assert delta_files_ag100.expected_file_exists('delta-input')

    def test_deltas_zip_created(self, delta_files_ag100):
        """Check that an archive with delta-amplitude files was created."""
        assert delta_files_ag100.expected_file_exists('Deltas/Deltas_001.zip')


# Notice that the compile_delta and compile_refcalc functions are
# virtually identical, except for a few error messages and a few
# variable names. There is no need to add more tests. When doing
# #43, the tests can be given to the base-class method!
# pylint: disable-next=too-few-public-methods
class TestCompileDelta(TestCompileRefcalc):
    """Tests for the compile_delta function."""

    compile_func = compile_delta
    compiler_cls_name = 'DeltaCompileTask'
    default_sources = 'delta.f', 'lib.tleed.f', 'lib.delta.f', 'GLOBAL'
    section_name = 'delta-amplitudes'
