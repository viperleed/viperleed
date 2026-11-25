"""Test configuration for tests.calc.sections.initialization.

Fixtures
--------
make_rpars (factory)
    Rparams with a BEAM_EQUIVALENCE loaded.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-25'
__license__ = 'GPLv3+'

from pytest_cases import fixture

from viperleed.calc.classes.beam_correspondence import BeamCorrespondence
from viperleed.calc.classes.rparams.rparams import Rparams

_BEAMLIST = ['  0  0\n', '  1  0\n', '  0  1\n', '  1  1\n']


@fixture(name='make_rpars')
def factory_make_rpars(make_beams):
    """Return a factory of Rparams with a BEAM_EQUIVALENCE loaded."""
    def _make(assignments, expbeams=(), ivbeams=()):
        rpars = Rparams()
        rpars.beamlist = list(_BEAMLIST)
        rpars.fileLoaded['BEAM_EQUIVALENCE'] = True
        rpars.beam_correspondence = BeamCorrespondence(assignments)
        rpars.expbeams = make_beams(*expbeams)
        rpars.fileLoaded['EXPBEAMS'] = bool(expbeams)
        rpars.ivbeams = make_beams(*ivbeams)
        return rpars
    return _make
