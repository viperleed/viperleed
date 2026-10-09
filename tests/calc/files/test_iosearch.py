"""Tests for module iosearch of viperleed.calc.files."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-10-09'
__license__ = 'GPLv3+'

import fortranformat as ff
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.files import iosearch


_COMMENT = '1: use RPe -- 2: use R2\n'


@parametrize('name,fmt,code', (
    ('pendry', 'I5', '    1'),
    ('r2', 'I5', '    2'),
    ('zj', 'I5', '    3'),
    ('smooth', 'I5', '    4'),
    ('pendry', 'I3', '  1'),
    ('smooth', 'I3', '  4'),
    ))
def test_r_factor_type_line(name, fmt, code):
    """Check the TensErLEED R-factor code written to search.steu."""
    # pylint: disable-next=protected-access
    line = iosearch._format_r_factor_type_line(RFactorType(name),
                                               ff.FortranRecordWriter(fmt))
    assert line == code.ljust(16) + _COMMENT
