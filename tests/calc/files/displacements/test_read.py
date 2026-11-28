"""Tests for viperleed.calc.files.displacements."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-06-08'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize_with_cases

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.files.displacements import NoDisplacementsError
from viperleed.calc.files.displacements import readDISPLACEMENTS

from .cases_read import CasesEmptyFile
from .cases_read import CasesInvalidDomainDisplacements


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


class TestReadDISPLACEMENTSMissingModeHeader:
    """Tests for blocks with assignments but no GEO/VIB/OCC tag."""

    def test_block_without_mode_header(self, write_displacements):
        """Check complaints for a search block without a *_DELTA tag."""
        file_path = write_displacements(
            '== SEARCH L1-4 x\n'
            '\n'
            '* L(1-4) [1 0] = -0.05 0.05 0.01\n'
            )
        with pytest.raises(NoDisplacementsError, match='no GEO_DELTA'):
            readDISPLACEMENTS(Rparams(), file_path)

    def test_block_without_mode_header_domains(
            self, domains, write_displacements):
        """Check complaints for a DOMAIN block without a *_DELTA tag."""
        file_path = write_displacements(
            '== DOMAIN Domain 0\n'
            '* L(1-4) z = -0.05 0.05 0.01\n'
            )
        with pytest.raises(NoDisplacementsError, match='no GEO_DELTA'):
            readDISPLACEMENTS(domains, file_path)

    def test_valid_block(self, write_displacements):
        """Check that a block with a mode header is read correctly."""
        file_path = write_displacements(
            '== SEARCH L1-4 x\n'
            '= GEO_DELTA\n'
            '* L(1-4) xy[1 0] = -0.05 0.05 0.01\n'
            )
        rpars = Rparams()
        readDISPLACEMENTS(rpars, file_path)
        assert len(rpars.disp_blocks) == 1
