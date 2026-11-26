"""Equivalence gate for BEAM_EQUIVALENCE files at initialization.

A BEAM_EQUIVALENCE file that reproduces the automatic beam
correspondence must give the same R-factor inputs as no file at all.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-25'
__license__ = 'GPLv3+'

import math
import shutil

import numpy as np
import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.beam import Beam
from viperleed.calc.classes.beam_correspondence import BeamCorrespondence
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_WORK
from viperleed.calc.files.beam_equivalence import read_beam_equivalence
from viperleed.calc.files.beam_equivalence import write_beam_equivalence
from viperleed.calc.files.iorfactor import writeWEXPEL
from viperleed.calc.files.iosearch import writeRfInfo
from viperleed.calc.files.tenserleed import CURRENT_TL_VERSION
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.leedbase import getBeamCorrespondence
from viperleed.calc.run import run_calc
from viperleed.calc.sections.calc_section import BEAM_EQUIVALENCE_NAME
from viperleed.calc.sections.run_sections import section_loop

pytestmark = pytest.mark.tensorleed

_IR_ENERGIES = range(50, 301, 3)
_IR_BEAMS = ('1| 0', '0| 1', '1/2| 0', '3/2| 0', '1/2| 1', '1| 1')


def _ag100_inputs(folder, data_path):
    """Copy the Ag(100) input files to folder."""
    for source in (data_path / 'Ag(100)' / 'initialization').iterdir():
        if source.is_file():
            shutil.copy(source, folder)


def _ir100_2x1_o_inputs(folder, data_path):
    """Write the inputs of Ir(100)-(2x1)-O (fractional beams)."""
    shutil.copy(data_path / 'POSCARs' / 'POSCAR_Ir(100)-(2x1)-O',
                folder / 'POSCAR')
    shutil.copy(data_path / 'parameters' / 'PARAMETERS_Ir(100)-(2x1)-O',
                folder / 'PARAMETERS')
    header = 'E,' + ','.join(f'({label})' for label in _IR_BEAMS)
    rows = [
        f'{energy:.2f},' + ','.join(
            f'{100 + 50*math.sin(energy/(10 + n)):.5E}'
            for n, _ in enumerate(_IR_BEAMS)
            )
        for energy in _IR_ENERGIES
        ]
    (folder / 'EXPBEAMS.csv').write_text('\n'.join([header, *rows]) + '\n')


_CASES = {
    'Ag(100)': _ag100_inputs,
    'Ir(100)-(2x1)-O': _ir100_2x1_o_inputs,
    }


@fixture(name='run_initialization')
def factory_run_initialization(mocker, data_path, make_section_tempdir,
                               tensorleed_path):
    """Return a function running initialization in a new folder.

    The function copies the inputs of a case, plus the files it is
    given, runs initialization only, and returns the Slab, the Rparams,
    and the work folder (with the SUPP folder).
    """
    loop = mocker.patch('viperleed.calc.run.section_loop',
                        wraps=section_loop)

    def _run(case, extra_files=None):
        home = make_section_tempdir(case, 'beam_equivalence_gate')
        work = home / DEFAULT_WORK
        work.mkdir()
        _CASES[case](home, data_path)
        for name, text in (extra_files or {}).items():
            (home / name).write_text(text)
        for source in home.iterdir():
            if source.is_file():
                shutil.copy(source, work)
        with execute_in_dir(work):
            failed, _ = run_calc(
                system_name=case,
                source=tensorleed_path,
                preset_params={'RUN': [0], 'TL_VERSION': CURRENT_TL_VERSION},
                home=home,
                )
        assert not failed
        rpars, slab = loop.call_args.args
        return slab, rpars, work
    return _run


def _r_factor_inputs(slab, rpars, folder):
    """Return the vector, rf.info, and WEXPEL for the current beams."""
    start, stop, step = rpars.THEO_ENERGIES.as_floats()
    theobeams = [Beam(beam.hk) for beam in rpars.ivbeams]
    for beam in theobeams:
        beam.intens = dict.fromkeys(np.arange(start, stop + step/2, step), 1.0)
    vector = getBeamCorrespondence(slab, rpars)
    with execute_in_dir(folder):  # writeAUXEXPBEAMS writes to cwd
        rf_info = writeRfInfo(slab, rpars, folder / 'rf.info')
        writeWEXPEL(slab, rpars, theobeams, folder / 'WEXPEL')
    return vector, rf_info, (folder / 'WEXPEL').read_text()


def _swapped(correspondence):
    """Return a correspondence with two assignments exchanged."""
    (key1, key2, *_), values = (
        list(correspondence.assignments),
        list(correspondence.assignments.values()),
        )
    assignments = dict(correspondence.assignments)
    assignments[key1], assignments[key2] = values[1], values[0]
    return BeamCorrespondence(assignments)


@parametrize(case=_CASES)
class TestEquivalenceGate:
    """Tests for the R-factor inputs with and without a file."""

    def test_reproducing_file(self, case, run_initialization):
        """Check identical R-factor inputs for a reproducing file."""
        slab, rpars, work = run_initialization(case)
        assert not rpars.fileLoaded['BEAM_EQUIVALENCE']
        automatic = _r_factor_inputs(slab, rpars, work)
        supp_file = work / DEFAULT_SUPP / BEAM_EQUIVALENCE_NAME

        slab, rpars, work = run_initialization(
            case,
            {BEAM_EQUIVALENCE_NAME: supp_file.read_text()}
            )
        assert rpars.fileLoaded['BEAM_EQUIVALENCE']
        from_file = _r_factor_inputs(slab, rpars, work)
        assert from_file[0] == automatic[0]
        assert from_file[1] == automatic[1]
        assert from_file[2] == automatic[2]

    def test_file_is_used(self, case, run_initialization, tmp_path):
        """Check that a file other than the automatic one is used."""
        slab, rpars, work = run_initialization(case)
        automatic = _r_factor_inputs(slab, rpars, work)
        supp_file = work / DEFAULT_SUPP / BEAM_EQUIVALENCE_NAME
        changed = _swapped(read_beam_equivalence(supp_file))
        write_beam_equivalence(changed, tmp_path / BEAM_EQUIVALENCE_NAME)

        slab, rpars, work = run_initialization(
            case,
            {BEAM_EQUIVALENCE_NAME:
             (tmp_path / BEAM_EQUIVALENCE_NAME).read_text()}
            )
        from_file = _r_factor_inputs(slab, rpars, work)
        assert rpars.beam_correspondence == changed
        assert from_file[0] != automatic[0]
        assert from_file[1] != automatic[1]
