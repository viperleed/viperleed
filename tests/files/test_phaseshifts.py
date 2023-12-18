"""Tests for module viperleed.calc.files.phaseshifts.

Created on 2023-07-28

@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)
"""

from pathlib import Path
import sys

from viperleed.calc.files import phaseshifts


def test_write_phaseshifts(run_phaseshift):
    """Ensure a PHASESHIFTS file is successfully written to disk."""
    param, _, firstline, phaseshift = run_phaseshift
    phaseshifts.writePHASESHIFTS(firstline, phaseshift,
                                 file_path=param.workdir/'PHASESHIFTS')
    assert any(param.workdir.glob('PHASESHIFTS'))
