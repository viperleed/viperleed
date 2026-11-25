"""Tests for beam files with a BEAM_EQUIVALENCE file loaded."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-25'
__license__ = 'GPLv3+'

from pytest_cases import fixture, parametrize

from viperleed.calc.classes.beam import Beam
from viperleed.calc.classes.beam_correspondence import BeamCorrespondence
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.files.beams import checkEXPBEAMS
from viperleed.calc.files.beams import sortIVBEAMS

_MODULE = 'viperleed.calc.files.beams'

# Beams of a p4mm surface: (1 0) and (0 1) are symmetry-equivalent
_P4MM_SYMEQ = (
    ((0.0, 0.0), 0),
    ((1.0, 0.0), 1), ((0.0, 1.0), 1), ((-1.0, 0.0), 1), ((0.0, -1.0), 1),
    )
_BEAMLIST = ['  0  0\n', '  1  0\n', '  0  1\n']


@fixture(name='rpars')
def fixture_rpars(mocker):
    """Return an Rparams with a p4mm symmetry list and a BEAMLIST."""
    mocker.patch(f'{_MODULE}.leedbase.getSymEqBeams',
                 return_value=list(_P4MM_SYMEQ))
    rpars = Rparams()
    rpars.HALTING = 3
    rpars.beamlist = list(_BEAMLIST)
    return rpars


def _use_file(rpars, use_file):
    """Load a BEAM_EQUIVALENCE into rpars if `use_file`."""
    if not use_file:
        return
    rpars.fileLoaded['BEAM_EQUIVALENCE'] = True
    rpars.beam_correspondence = BeamCorrespondence(
        {(1, 0): {(1, 0): 1}, (0, 1): {(0, 1): 1}}
        )


class TestSymmetryMergingSkipped:
    """Tests for the symmetry-equivalent beams of BEAM_EQUIVALENCE."""

    @parametrize(use_file=(True, False))
    def test_check_expbeams(self, use_file, rpars):
        """Check that equivalent experimental beams survive with a file."""
        _use_file(rpars, use_file)
        rpars.expbeams = [Beam((1, 0)), Beam((0, 1))]
        checkEXPBEAMS(None, rpars)
        assert len(rpars.expbeams) == (2 if use_file else 1)

    @parametrize(use_file=(True, False))
    def test_sort_ivbeams_both_in_beamlist(self, use_file, rpars):
        """Check that both p4mm-equivalent beams of BEAMLIST survive."""
        _use_file(rpars, use_file)
        rpars.ivbeams = [Beam((1, 0)), Beam((0, 1))]
        sorted_beams = sortIVBEAMS(None, rpars)
        assert [b.hkfrac for b in sorted_beams] == [
            (1, 0), (0, 1),
            ]

    @parametrize(use_file=(True, False))
    def test_sort_ivbeams_not_in_beamlist(self, use_file, rpars, caplog):
        """Check that an equivalent beam is renamed only without a file."""
        _use_file(rpars, use_file)
        rpars.beamlist = ['  0  0\n', '  1  0\n']
        rpars.ivbeams = [Beam((0, 1))]
        with caplog.at_level('DEBUG'):
            sorted_beams = sortIVBEAMS(None, rpars)
        if use_file:
            assert not sorted_beams
            assert 'not in BEAMLIST and will be dropped' in caplog.text
            assert 'renaming' not in caplog.text
        else:
            assert [b.hkfrac for b in sorted_beams] == [(1, 0)]
            assert 'renaming' in caplog.text

