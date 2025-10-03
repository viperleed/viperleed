"""Tests for viperleed.calc.files.displacements."""

__authors__ = ('Michele Riva (@michele-riva)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-06-08'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import (
    fixture,
    parametrize_with_cases,
)

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.files.displacements import (
    NoDisplacementsError,
    readDISPLACEMENTS,
)

from .cases_read import (
    CasesEmptyFile,
    CasesInvalidDomainDisplacements,
)


@fixture(name='write_displacements')
def fixture_write_displacements(tmp_path):
    """Write a temporary DISPLACEMENTS file with given contents."""
    file = tmp_path / 'DISPLACEMENTS'

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
