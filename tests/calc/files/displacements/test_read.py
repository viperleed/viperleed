"""Tests for viperleed.calc.files.displacements."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-06-08'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize
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


_PLANE_SHORTHAND_WARNING = 'deprecated shorthand for two search blocks'


class TestPlaneShorthandDeprecated:
    """Tests for the deprecated 'xy'/'ab' two-search-block shorthand."""

    @parametrize(
        ('contents', 'expected_blocks', 'n_warnings'),
        (
            ('== SEARCH s\n= GEO_DELTA\nO 1 xy = -0.1 0.1\n',
             [(['= GEO_DELTA', 'O 1 xy[1 0] = -0.1 0.1'], 's [1 0]'),
              (['= GEO_DELTA', 'O 1 xy[0 1] = -0.1 0.1'], 's [0 1]')],
             1),
            (('== SEARCH s\n= GEO_DELTA\nO 1 ab = -0.1 0.1\n'
              'Fe 2 z = -0.1 0.1\nFe 3 ab = -0.2 0.2\n'),
             [(['= GEO_DELTA', 'O 1 ab[1 0] = -0.1 0.1',
                'Fe 2 z = -0.1 0.1', 'Fe 3 ab[1 0] = -0.2 0.2'],
               's [1 0]'),
              (['= GEO_DELTA', 'O 1 ab[0 1] = -0.1 0.1',
                'Fe 2 z = -0.1 0.1', 'Fe 3 ab[0 1] = -0.2 0.2'],
               's [0 1]')],
             2),
            ('== SEARCH s\n= GEO_DELTA\nO 1 xy[1 0] = -0.1 0.1\n',
             [(['= GEO_DELTA', 'O 1 xy[1 0] = -0.1 0.1'], 's')],
             0),
        ),
        ids=('xy', 'ab twice', 'xy[1 0]'),
    )
    def test_warning(self, contents, expected_blocks, n_warnings,
                     write_displacements, caplog):
        """Check one warning per shorthand line and the blocks made."""
        rpars = Rparams()
        readDISPLACEMENTS(rpars, write_displacements(contents))
        warnings = [r for r in caplog.records
                    if _PLANE_SHORTHAND_WARNING in r.getMessage()]
        assert len(warnings) == n_warnings
        assert all(r.levelname == 'WARNING' for r in warnings)
        assert rpars.disp_blocks == expected_blocks
