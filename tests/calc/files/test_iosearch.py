"""Tests for module iosearch of viperleed.calc.files."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-09'
__license__ = 'GPLv3+'

import fortranformat as ff
import pytest
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.files import iosearch
from viperleed.calc.lib.version import Version


_COMMENT = '1: use RPe -- 2: use R2 -- 4: use RS\n'
_TL_VERSIONS = ('1.6.0', '1.7.3', '2.0.0', '2.1.0')


def _format_line(name, fmt, tl_version):
    """Return the search.steu line for R-factor type `name`."""
    # pylint: disable-next=protected-access
    return iosearch._format_r_factor_type_line(RFactorType(name),
                                               ff.FortranRecordWriter(fmt),
                                               Version(tl_version))


@parametrize('name,fmt,code', (
    ('pendry', 'I5', '    1'),
    ('r2', 'I5', '    2'),
    ('smooth', 'I5', '    4'),
    ('pendry', 'I3', '  1'),
    ('smooth', 'I3', '  4'),
    ))
def test_r_factor_type_line(name, fmt, code):
    """Check the TensErLEED R-factor code written to search.steu."""
    line = _format_line(name, fmt, '2.1.0')
    assert line == code.ljust(16) + _COMMENT


@parametrize(name=('pendry', 'r2'))
@parametrize(tl_version=_TL_VERSIONS)
def test_r_factor_type_line_any_version(name, tl_version):
    """Check that Pendry and R2 are written for any TensErLEED."""
    code = str(RFactorType(name).value).rjust(5)
    line = _format_line(name, 'I5', tl_version)
    assert line == code.ljust(16) + _COMMENT


_UNSUPPORTED = (  # name, tl_version, match
    *(('zj', version, 'ZJ cannot be used') for version in _TL_VERSIONS),
    ('smooth', '1.7.3', r'SMOOTH needs TensErLEED 2\.1\.0.*is 1\.7\.3'),
    ('smooth', '2.0.0', r'SMOOTH needs TensErLEED 2\.1\.0.*is 2\.0\.0'),
    )


@parametrize('name,tl_version,match', _UNSUPPORTED)
def test_r_factor_type_line_unsupported(name, tl_version, match):
    """Check that the search refuses R-factor types it cannot run."""
    with pytest.raises(ValueError, match=match):
        _format_line(name, 'I5', tl_version)
