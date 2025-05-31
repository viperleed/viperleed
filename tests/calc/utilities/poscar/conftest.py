"""Test configuration for tests/calc/utilities/poscar.

Fixtures
--------
poscar_stream
    Mimic users piping a POSCAR file to sys.stdin.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-03-20'
__license__ = 'GPLv3+'

import io

from pytest_cases import fixture


@fixture
def poscar_stream(mocker, poscars_path):
    """Replace sys.stdin with a readable stream read from a POSCAR file."""
    def _patch(poscar_name):
        poscar_file = poscars_path/poscar_name
        return mocker.patch('sys.stdin', io.StringIO(poscar_file.read_text()))
    return _patch
