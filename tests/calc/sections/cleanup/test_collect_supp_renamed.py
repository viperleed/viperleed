"""Tests for storing renamed work-folder files in SUPP."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-25'
__license__ = 'GPLv3+'

import logging

from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.calc_section import BEAM_EQUIVALENCE_NAME
from viperleed.calc.sections.calc_section import BEAM_EQUIVALENCE_USED_NAME
from viperleed.calc.sections.cleanup import _collect_supp_contents


def test_used_beam_equivalence_in_supp(rpars, workdir):
    """Check that the used file lands in SUPP under the input name."""
    (workdir/BEAM_EQUIVALENCE_NAME).write_text('input')
    (workdir/BEAM_EQUIVALENCE_USED_NAME).write_text('used')
    with execute_in_dir(workdir):
        _collect_supp_contents(rpars)
    assert (workdir/DEFAULT_SUPP/BEAM_EQUIVALENCE_NAME).read_text() == 'used'
    assert not (workdir/DEFAULT_SUPP/BEAM_EQUIVALENCE_USED_NAME).exists()
    assert (workdir/BEAM_EQUIVALENCE_NAME).read_text() == 'input'
    assert DEFAULT_SUPP in rpars.manifest
    assert BEAM_EQUIVALENCE_USED_NAME not in rpars.manifest
    assert BEAM_EQUIVALENCE_NAME not in rpars.manifest


def test_used_beam_equivalence_absent(rpars, workdir):
    """Check that nothing is stored without a used file."""
    (workdir/BEAM_EQUIVALENCE_NAME).write_text('input')
    with execute_in_dir(workdir):
        _collect_supp_contents(rpars)
    assert not (workdir/DEFAULT_SUPP/BEAM_EQUIVALENCE_NAME).exists()
    assert DEFAULT_SUPP not in rpars.manifest


def test_used_beam_equivalence_copy_fails(rpars, workdir, mocker, caplog):
    """Check that an OSError while copying is logged, not raised."""
    (workdir/BEAM_EQUIVALENCE_USED_NAME).write_text('used')
    mocker.patch('shutil.copy2', side_effect=OSError)
    with execute_in_dir(workdir), caplog.at_level(logging.ERROR):
        _collect_supp_contents(rpars)
    assert f'file {BEAM_EQUIVALENCE_USED_NAME}' in caplog.text
    assert DEFAULT_SUPP not in rpars.manifest
