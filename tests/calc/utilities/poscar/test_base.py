"""Tests for module base of viperleed.utilities.poscar."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-03-20'
__license__ = 'GPLv3+'

import io

from viperleed.utilities.poscar.project_c_to_z import ProjectCToZCLI
from viperleed.utilities.poscar.sort_by_z import SortByZCLI
from viperleed.utilities.poscar.strip_comments import StripCommentsCLI


def test_pipe_utilities(poscar_stream, mocker, capsys):
    """Check correct piping of multiple 'viperleed poscar' utilities."""
    clis = StripCommentsCLI(), SortByZCLI(), ProjectCToZCLI()
    # We simulate piping by replacing stdin and stdout, ensuring
    # that stdin behaves like the terminal (isatty() == True).
    # (1) we start by reading in an example POSCAR to stdin, then
    stdin = poscar_stream('POSCAR_Ag(100)')
    mocker.patch.object(stdin, 'isatty', return_value=True)
    for cli in clis:
        cli([])
        # (2) We use the stdout of one cli call as stdin for the next
        stdout = capsys.readouterr().out
        stdin = mocker.patch('sys.stdin', io.StringIO(stdout))
        mocker.patch.object(stdin, 'isatty', return_value=True)
