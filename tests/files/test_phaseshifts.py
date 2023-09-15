"""Tests for module viperleed.tleedmlib.files.phaseshifts.

Created on 2023-07-28

@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)
"""

from pathlib import Path
import sys

# pylint: disable=unused-import
# This is a workaround to make sure we do point to the necessary
# plugin. Will be resolved in installable by adding as a dependency
import pytest_subtests
# pylint: enable=unused-import

VPR_PATH = str(Path(__file__).resolve().parents[3])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Will be fixed in installable version
from viperleed.tleedmlib.files import phaseshifts
# pylint: enable=wrong-import-position


def test_write_phaseshifts(run_phaseshift):
    """Ensure a PHASESHIFTS file is successfully written to disk."""
    param, _, firstline, phaseshift = run_phaseshift
    phaseshifts.writePHASESHIFTS(firstline, phaseshift,
                                 file_path=param.workdir/'PHASESHIFTS')
    assert any(param.workdir.glob('PHASESHIFTS'))
