"""Tests for viperleed.calc.files.displacements."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from pytest import approx
from pytest_cases import fixture
from pytest_cases import parametrize_with_cases

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.files.displacements import NoDisplacementsError
from viperleed.calc.files.displacements import readDISPLACEMENTS
from viperleed.calc.lib.context import execute_in_dir


@fixture(name='write_displacements')
def fixture_write_displacements(tmp_path):
    """Write a temporary DISPLACEMENTS file with given contents."""
    file = tmp_path/'DISPLACEMENTS'
    def _write(contents):
        file.write_text(contents, encoding='utf-8')
        return file
    return _write


@fixture(name='domains')
def fixture_make_domains(mocker):
    """Return an Rparams with domain info."""
    rpars = Rparams()
    for ind in range(3):
        mock_domain = mocker.MagicMock()
        mock_domain.name = f'Domain {ind}'
        rpars.domainParams.append(mock_domain)
    return rpars


class TestReadDISPLACEMENTS:
    """Test successful reading of a DISPLACEMENTS block."""

    def test_geo(self, displaced_atom):
        """Test successful reading of a geometric DISPLACEMENTS block."""
        displ_read = np.array(displaced_atom.disp_geo['all'])
        assert displ_read[:, 2] == approx([0.2, 0.1, 0.0, -0.1, -0.2])

    def test_geo_inverted_ranges(self, ag100_with_displacements_and_offsets):
        """Test reading geometric DISPLACEMENTS with start > stop."""
        slab, *_ = ag100_with_displacements_and_offsets
        atom = slab.atlist[1]
        displ_read = np.array(atom.disp_geo['all'])
        assert displ_read[:, 2] == approx([-0.2, -0.1, 0.0, 0.1, 0.2])

    def test_occ(self, displaced_atom):
        """Test successful reading of an occupation DISPLACEMENTS block."""
        displ_read = displaced_atom.disp_occ[displaced_atom.el]
        assert displ_read == approx([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

    def test_vib(self, displaced_atom):
        """Test successful reading of a vibration DISPLACEMENTS block."""
        displ_read = displaced_atom.disp_vib['all']
        assert displ_read == approx([-0.1, 0.0, 0.1])


class CasesEmptyFile:
    """Collection of cases of DISPLACEMENTS file with no valid content."""

    def case_comments_only(self, write_displacements):
        """Prepare a DISPLACEMENTS file with no contents."""
        contents = '''
# A comment line with a hash character
! == SEARCH commented with an exclamation mark
% = GEO_DELTA   ! A commented tag
#    Fe* L(1-5) z = -0.5 0.5 0.05   ! And a commented valid line
'''
        return write_displacements(contents)

    def case_empty_lines(self, write_displacements):
        """Prepare a DISPLACEMENTS file with no contents."""
        return write_displacements('\n  \n      \n')

    def case_only_loop_tags(self, write_displacements):
        """Prepare a DISPLACEMENTS file with only <loop> tags."""
        contents = '''
<loop>

</loop>
'''
        return write_displacements(contents)

    def case_no_lines(self, write_displacements):
        """Prepare a DISPLACEMENTS file with no contents."""
        return write_displacements('')


class CasesInvalidDomainDisplacements:
    """Invalid DISPLACEMENTS files for multi-domain calculations."""

    def case_all_domain_names_mismatched(self, domains, write_displacements):
        """Prepare a DISPLACEMENTS with invalid DOMAIN names."""
        contents = '''
== DOMAIN not-among-the-names
  = GEO_DELTA
     * L(1-5) xy[1 0] = -0.5 0.5 0.05'''
        return domains, write_displacements(contents)

    def case_only_blocks_outside_domain(self, domains, write_displacements):
        """Prepare a single-domain DISPLACEMENTS for a DOMAINs calculation."""
        contents = '''
= GEO_DELTA   # Not assigned to a DOMAIN, will be skipped
   * L(1-5) z = 0 0.2 0.05'''
        return domains, write_displacements(contents)

    def case_all_domain_blocks_empty(self, domains, write_displacements):
        """Prepare a DISPLACEMENTS with empty domain DISPLACEMENTS."""
        contents = '''
= DOMAIN 0
= DOMAIN 1
    # comments only
= DOMAIN 2
    = GEO_DELTA  ! empty'''
        return domains, write_displacements(contents)


class TestReadDISPLACEMENTSRaises:
    """Tests for exceptions raised by the readDISPLACEMENTS function."""

    @parametrize_with_cases('file_path', cases=CasesEmptyFile)
    def test_empty_file_fails(self, file_path):
        """Test complaints when an empty DISPLACEMENTS file is read."""
        rpars = Rparams()
        with pytest.raises(NoDisplacementsError):
            readDISPLACEMENTS(rpars, file_path)

    @parametrize_with_cases('args', cases=CasesInvalidDomainDisplacements)
    def test_invalid_domain_names(self, args):
        """Test complaints for invalid DISPLACEMENTS for DOMAINs."""
        with pytest.raises(NoDisplacementsError):
            readDISPLACEMENTS(*args)
