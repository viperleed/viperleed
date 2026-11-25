"""Tests for beam files with a BEAM_EQUIVALENCE file loaded."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-25'
__license__ = 'GPLv3+'

from pytest_cases import fixture
from pytest_cases import parametrize
from quicktions import Fraction

from viperleed.calc.classes.beam import Beam
from viperleed.calc.classes.beam_correspondence import BeamCorrespondence
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.files.beams import checkEXPBEAMS
from viperleed.calc.files.beams import readIVBEAMS
from viperleed.calc.files.beams import sortIVBEAMS
from viperleed.calc.files.beams import writeIVBEAMS

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
        """Check equivalent experimental beams survive with a file."""
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
        """Check an equivalent beam is renamed only without a file."""
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


class TestWriteIvbeamsFromFile:
    """Tests for the IVBEAMS generated from BEAM_EQUIVALENCE."""

    @parametrize(domains=(True, False))
    def test_union_of_right_hand_sides(self, domains, rpars, mocker, tmp_path):
        """Check IVBEAMS holds all calculated beams of the file."""
        rpars.fileLoaded['BEAM_EQUIVALENCE'] = True
        rpars.beam_correspondence = BeamCorrespondence({
            (1, 0): {(1, 0): 1, (0, 1): 1},
            (1, 1): {(0, 1): 1, (Fraction(1, 2), 0): 1},
            })
        rpars.expbeams = [Beam((1, 0)), Beam((1, 1))]
        leed_dict = mocker.patch(f'{_MODULE}.leedbase.getLEEDdict')
        filename = tmp_path/'IVBEAMS'
        ivbeams = writeIVBEAMS(None, rpars, filename=filename,
                               domains=domains)
        expect = [(1, 0), (0, 1), (Fraction(1, 2), 0)]
        assert [b.hkfrac for b in ivbeams] == expect
        assert [b.hkfrac for b in readIVBEAMS(filename)] == expect
        assert 'BEAM_EQUIVALENCE' in filename.read_text().splitlines()[0]
        leed_dict.assert_not_called()

    def test_sorted_in_beamlist_order(self, rpars, tmp_path):
        """Check that the generated beams are sorted by sortIVBEAMS."""
        rpars.fileLoaded['BEAM_EQUIVALENCE'] = True
        rpars.beam_correspondence = BeamCorrespondence(
            {(0, 1): {(0, 1): 1}, (1, 0): {(1, 0): 1}}
            )
        rpars.ivbeams = writeIVBEAMS(None, rpars, filename=tmp_path/'IVBEAMS')
        assert [b.hkfrac for b in rpars.ivbeams] == [(0, 1), (1, 0)]
        sorted_beams = sortIVBEAMS(None, rpars)
        assert [b.hkfrac for b in sorted_beams] == [(1, 0), (0, 1)]
