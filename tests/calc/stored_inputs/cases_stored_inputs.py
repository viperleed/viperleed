"""Cases for the tests of tests.calc.stored_inputs.

Each case function returns the StoredInputs of a case, i.e., where its
input files come from, and is tagged with the calculation steps whose
inputs are checked for it: REFCALC (reference calculation), DELTAS
(delta amplitudes) and SEARCH. The tests of each step select their
cases by tag. The id of a case is the name of the folder of its
references, with the name of the domain appended for domain cases.

The inputs come from either
- a folder of stored inputs (StoredInputs.source) in the test data,
  all of which are copied to the folders of the calculation, or
- a slab and Rparams of the POSCAR presets (StoredInputs.preset, see
  tests/calc/poscar_slabs.py), with synthetic PHASESHIFTS and the
  PRESET_FILES taken from the folder StoredInputs.source.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-13'
__license__ = 'GPLv3+'

import dataclasses
import sys

from pytest_cases import case
from pytest_cases import get_all_cases
from pytest_cases import get_case_id

from ...helpers import TEST_DATA

REFCALC = 'refcalc'
DELTAS = 'deltas'
SEARCH = 'search'
DISPLACEMENTS = 'DISPLACEMENTS'
PRESET_FILES = 'IVBEAMS', 'VIBROCC'
AG100 = 'Ag(100)/initialization'
AG100_DISPLACEMENTS = TEST_DATA / 'Ag(100)/displacements'
PTRH_INPUTS = 'stored_inputs/inputs/ptrh_100_1x3_o'


@dataclasses.dataclass(frozen=True)
class StoredInputs:
    """Where the input files of a case come from.

    Attributes
    ----------
    source : str
        A folder in the test data. If `preset` is not given, all the
        files in it are copied to the folders of the calculation.
        Otherwise, only the PRESET_FILES are.
    preset : str, optional
        The name of a case method of CasePOSCARSlabs that returns the
        slab and Rparams of the case. Default is no preset.
    extra : tuple, optional
        (destination, source) pairs of files of the test data that
        are added on top of those in `source`. Default is no files.
    domains : tuple, optional
        For domain cases, the names of the domains, i.e., of the
        sub-folders of `source` to copy. Default is a single-surface
        case.
    displacements : Path or None, optional
        The DISPLACEMENTS file to use for the deltas and the search.
        Default is None, i.e., use the one in `source`.
    n_param : int, optional
        The number of delta PARAM files, i.e., of different array
        dimensions, of each domain. Default is 1.
    """

    source: str
    preset: str = ''
    extra: tuple = ()
    domains: tuple = ()
    displacements: object = None
    n_param: int = 1

    @property
    def domain_names(self):
        """Return the domains, or (None,) for a single-surface case."""
        return self.domains or (None,)

    @property
    def initialization(self):
        """Return the inputs that the initialization depends on."""
        return dataclasses.replace(self, displacements=None, n_param=1)


def tagged_inputs(tag):
    """Return the ids and StoredInputs of the cases tagged with `tag`.

    For tests whose parametrization depends on the inputs of the case
    (e.g., one test per domain), which parametrize_with_cases cannot
    express.

    Parameters
    ----------
    tag : str
        The tag of the cases to return.

    Returns
    -------
    list
        (case id, StoredInputs) pairs, in the order of the cases.
    """
    cases = get_all_cases(cases=sys.modules[__name__], has_tag=tag)
    return [(get_case_id(case_fn), case_fn()) for case_fn in cases]


@case(tags=REFCALC)
def case_ag100():
    """Return Ag(100), as stored."""
    return StoredInputs(AG100)


@case(tags=REFCALC)
def case_ag100_el_rename():
    """Return Ag(100) with a renamed chemical element."""
    return StoredInputs('Ag(100)_el_rename/initialization')


@case(tags=(DELTAS, SEARCH))
def case_ag100_z():
    """Return Ag(100) with displacements along z."""
    return StoredInputs(
        AG100,
        displacements=AG100_DISPLACEMENTS / 'DISPLACEMENTS_z',
        )


@case(tags=(DELTAS, SEARCH))
def case_ag100_vib():
    """Return Ag(100) with vibration displacements."""
    return StoredInputs(
        AG100,
        displacements=AG100_DISPLACEMENTS / 'DISPLACEMENTS_vib',
        )


@case(id='ag100_z+vib', tags=(DELTAS, SEARCH))
def case_ag100_z_vib():
    """Return Ag(100) with displacements along z and vibrations."""
    return StoredInputs(
        AG100,
        displacements=AG100_DISPLACEMENTS / 'DISPLACEMENTS_z+vib',
        )


@case(tags=(DELTAS, SEARCH))
def case_ag100_mixed():
    """Return Ag(100) with VIBROCC offsets and mixed displacements.

    The VIBROCC offsets and DISPLACEMENTS only matter after a search:
    the inputs of the reference calculation are those of case_ag100.
    """
    return StoredInputs(
        AG100,
        extra=(('VIBROCC', 'Ag(100)/mergeDisp/VIBROCC'),
               (DISPLACEMENTS, 'Ag(100)/mergeDisp/DISPLACEMENTS_mixed')),
        n_param=3,
        )


@case(tags=(REFCALC, DELTAS, SEARCH))
def case_telluride_stacking():
    """Return a telluride with two stacking domains."""
    return StoredInputs(
        'domains/telluride_stacking/deltas',
        domains=('STA', 'FAU'),
        n_param=2,
        )


@case(tags=(REFCALC, DELTAS, SEARCH))
def case_ptrh_100_1x3_o():
    """Return the PtRh(100)-(1x3)-O preset."""
    return StoredInputs(
        PTRH_INPUTS,
        preset='case_poscar_ptrh_100_1x3_o',
        displacements=TEST_DATA / PTRH_INPUTS / DISPLACEMENTS,
        n_param=3,
        )
