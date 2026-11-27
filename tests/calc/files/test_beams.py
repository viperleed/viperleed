"""Tests for module viperleed.calc.files.beams."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

from pytest_cases import parametrize

from viperleed.calc.files import beams

_EXPECTED_HK = [(1, 0), (0, 1), (0.5, 0.5)]
_COMMENTED = {
    'bang': """\
h k   ! header
1 0   # trailing hash
# 0 2
0 1 % trailing percent

% 2 0
1/2 1/2!no space
""",
    'blank-header': """\

1 0
#1 1

0 1
1/2 1/2
""",
    'hash-header': """\
# h k
1 0
! 1 1
0 1
1/2 1/2
""",
    'whitespace-header': """\
\t
1 0
0 1  ! comment
1/2 1/2
""",
    }


@parametrize('text', _COMMENTED.values(), ids=_COMMENTED)
def test_read_ivbeams_comments(text, tmp_path):
    """Check that comments in IVBEAMS are ignored."""
    ivbeams = tmp_path / 'IVBEAMS'
    ivbeams.write_text(text, encoding='utf-8')
    assert [b.hk for b in beams.readIVBEAMS(ivbeams)] == _EXPECTED_HK
