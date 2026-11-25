"""Tests for the final beam correspondence at initialization."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-25'
__license__ = 'GPLv3+'

import logging

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.beam_correspondence import BeamCorrespondence
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.files.beam_equivalence import read_beam_equivalence
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.calc_section import BEAM_EQUIVALENCE_NAME
from viperleed.calc.sections.calc_section import BEAM_EQUIVALENCE_USED_NAME
from viperleed.calc.sections.cleanup import _collect_supp_contents
from viperleed.calc.sections.initialization import BeamEquivalenceError
from viperleed.calc.sections.initialization import (
    _finalize_beam_correspondence,
    )

# pylint: disable-next=unused-import
from .test_beam_equivalence_checks import factory_make_rpars  # noqa: F401

_LEEDBASE = 'viperleed.calc.lib.leedbase'
_IVBEAMS = [(1, 0), (0, 1), (1, 1)]
_ASSIGNMENTS = {
    'weighted': {(1, 0): {(1, 0): 1, (0, 1): 2}},
    'shared': {(1, 0): {(1, 0): 1}, (1, 1): {(1, 0): 1, (0, 1): 1}},
    'equal weights': {(1, 0): {(1, 0): 2, (0, 1): 2}, (1, 1): {(1, 1): 3}},
    }
_RUNS = {
    'refcalc': [1],
    'search': [3],
    'superpos': [31],
    'error calc': [4],
    'no R factor': [2, 5],
    }


@fixture(name='rpars_for')
def factory_rpars_for(make_rpars):
    """Return a factory of Rparams with EXPBEAMS and a file."""
    def _make(assignments, run):
        rpars = make_rpars(assignments,
                           expbeams=list(assignments),
                           ivbeams=_IVBEAMS)
        rpars.RUN = list(run)
        return rpars
    return _make


@fixture(name='in_workdir')
def fixture_in_workdir(tmp_path):
    """Run the test in a temporary work directory."""
    with execute_in_dir(tmp_path):
        yield tmp_path


class TestWeightCheck:
    """Tests for the check of weights and shared calculated beams."""

    _run = staticmethod(_finalize_beam_correspondence)

    @parametrize('kind', ('weighted', 'shared'))
    @parametrize('run', _RUNS)
    @parametrize('legacy', (True, False))
    def test_error(self, kind, run, legacy, rpars_for, in_workdir, caplog):
        """Check an error when an R factor runs, even if not legacy."""
        rpars = rpars_for(_ASSIGNMENTS[kind], _RUNS[run])
        rpars.R_FACTOR_LEGACY = legacy
        expect_error = run != 'no R factor'
        if expect_error:
            with pytest.raises(BeamEquivalenceError):
                self._run(None, rpars)
        else:
            self._run(None, rpars)
        errors = [r.message for r in caplog.records
                  if r.levelno == logging.ERROR]
        assert bool(errors) == expect_error
        if expect_error:
            (message,) = errors
            assert message == self._expected_message(kind)

    @staticmethod
    def _expected_message(kind):
        cause = {
            'weighted': 'assignment weights',
            'shared': 'calculated beams in more than one beam assignment',
            }[kind]
        lines = {
            'weighted': '  (1 0) = 1/3*(1 0) + 2/3*(0 1)',
            'shared': '  (1 0) = (1 0)\n  (1 1) = (1 0) + (0 1)',
            }[kind]
        return (f'BEAM_EQUIVALENCE contains {cause}, which the R factor '
                f'to be run does not support. Offending lines:\n{lines}')

    @parametrize('run', _RUNS)
    def test_equal_weights_ok(self, run, rpars_for, in_workdir, caplog):
        """Check that equal weights within each line are no problem."""
        rpars = rpars_for(_ASSIGNMENTS['equal weights'], _RUNS[run])
        self._run(None, rpars)
        assert not [r for r in caplog.records if r.levelno >= logging.ERROR]

    def test_both_causes(self, rpars_for, in_workdir):
        """Check that the message names both causes and only bad lines."""
        assignments = {
            (1, 0): {(1, 0): 1, (0, 1): 2},
            (1, 1): {(0, 1): 1},
            (2, 0): {(1, 1): 1},
            }
        rpars = rpars_for(assignments, [1])
        with pytest.raises(BeamEquivalenceError) as exc_info:
            self._run(None, rpars)
        message = str(exc_info.value)
        assert ('assignment weights and calculated beams in more than one '
                'beam assignment') in message
        assert '(1 0) = 1/3*(1 0) + 2/3*(0 1)' in message
        assert '(1 1) = (0 1)' in message
        assert '(2 0)' not in message


class TestUsedFile:
    """Tests for the correspondence written for SUPP."""

    def test_file_from_beam_equivalence(self, rpars_for, in_workdir):
        """Check the written file for a file-based correspondence."""
        rpars = rpars_for(_ASSIGNMENTS['equal weights'], [2])
        _finalize_beam_correspondence(None, rpars)
        used = read_beam_equivalence(in_workdir / BEAM_EQUIVALENCE_USED_NAME)
        assert used == rpars.beam_correspondence

    @parametrize('sl', (None, object()), ids=('domains', 'single'))
    def test_automatic_to_supp(self, sl, rpars_for, in_workdir, mocker):
        """Check SUPP contents after cleanup for an automatic grouping."""
        rpars = rpars_for(_ASSIGNMENTS['equal weights'], [2])
        rpars.fileLoaded['BEAM_EQUIVALENCE'] = False
        rpars.beam_correspondence = None
        mocker.patch(f'{_LEEDBASE}.getSymEqBeams', return_value=[])
        _finalize_beam_correspondence(sl, rpars)
        with execute_in_dir(in_workdir):
            _collect_supp_contents(rpars)
        supp = read_beam_equivalence(
            in_workdir / DEFAULT_SUPP / BEAM_EQUIVALENCE_NAME
            )
        assert supp == BeamCorrespondence({(1, 0): {(1, 0): 1},
                                           (1, 1): {(1, 1): 1}})
        assert DEFAULT_SUPP in rpars.manifest
        assert BEAM_EQUIVALENCE_USED_NAME not in rpars.manifest

    @parametrize('sl', (None, object()), ids=('domains', 'single'))
    def test_file_to_supp(self, sl, rpars_for, in_workdir):
        """Check SUPP contents, and that the input file is untouched."""
        input_file = in_workdir / BEAM_EQUIVALENCE_NAME
        input_file.write_text('# input\n(1 0) = (1 0)\n')
        rpars = rpars_for(_ASSIGNMENTS['equal weights'], [2])
        _finalize_beam_correspondence(sl, rpars)
        _collect_supp_contents(rpars)
        supp = read_beam_equivalence(
            in_workdir / DEFAULT_SUPP / BEAM_EQUIVALENCE_NAME
            )
        assert supp == rpars.beam_correspondence
        assert input_file.read_text() == '# input\n(1 0) = (1 0)\n'

    def test_empty_not_written(self, rpars_for, in_workdir, mocker, caplog):
        """Check that an empty correspondence is not written."""
        rpars = rpars_for(_ASSIGNMENTS['equal weights'], [2])
        rpars.fileLoaded['BEAM_EQUIVALENCE'] = False
        rpars.beam_correspondence = None
        rpars.ivbeams = []
        mocker.patch(f'{_LEEDBASE}.getSymEqBeams', return_value=[])
        _finalize_beam_correspondence(None, rpars)
        assert not rpars.expbeams
        assert not (in_workdir / BEAM_EQUIVALENCE_USED_NAME).exists()
        assert 'No experimental beam has a calculated beam' in caplog.text

    def test_no_expbeams(self, make_rpars, in_workdir):
        """Check that nothing is done without EXPBEAMS."""
        rpars = make_rpars({(1, 0): {(1, 0): 1}})
        _finalize_beam_correspondence(None, rpars)
        assert not (in_workdir / BEAM_EQUIVALENCE_USED_NAME).exists()
