"""Tests for loading BEAM_EQUIVALENCE into Rparams."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-25'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.sections.calc_section import BEAM_EQUIVALENCE_NAME


@fixture(name='rpars')
def fixture_rpars():
    """Return a bare Rparams."""
    return Rparams()


@fixture(name='in_tmp_dir')
def fixture_in_tmp_dir(tmp_path, monkeypatch):
    """Run the test in an empty temporary directory."""
    monkeypatch.chdir(tmp_path)
    return tmp_path


class TestTryLoadingBeamEquivalenceFile:
    """Tests for Rparams.try_loading_beam_equivalence_file."""

    def test_default_not_loaded(self, rpars):
        """Check that nothing is loaded by default."""
        assert not rpars.fileLoaded['BEAM_EQUIVALENCE']
        assert rpars.beam_correspondence is None

    def test_file_present(self, rpars, in_tmp_dir):
        """Check that a valid file is read and stored."""
        (in_tmp_dir / BEAM_EQUIVALENCE_NAME).write_text(
            '(1 0) = (1 0)\n(0 1) = (0 1) + (0 -1)\n',
            encoding='utf-8',
            )
        rpars.try_loading_beam_equivalence_file()
        assert rpars.fileLoaded['BEAM_EQUIVALENCE']
        assert len(rpars.beam_correspondence.assignments) == 2

    @pytest.mark.usefixtures('in_tmp_dir')
    def test_file_absent(self, rpars, caplog):
        """Check that a missing file is silently ignored."""
        with caplog.at_level('DEBUG'):
            rpars.try_loading_beam_equivalence_file()
        assert not rpars.fileLoaded['BEAM_EQUIVALENCE']
        assert rpars.beam_correspondence is None
        assert not [r for r in caplog.records if r.levelname != 'DEBUG']

    def test_malformed(self, rpars, in_tmp_dir, caplog):
        """Check that a malformed file is logged as error and raised."""
        (in_tmp_dir / BEAM_EQUIVALENCE_NAME).write_text(
            '(1 0) = (1 0)\nthis is not valid\n',
            encoding='utf-8',
            )
        with pytest.raises(ValueError, match='line 2'):
            rpars.try_loading_beam_equivalence_file()
        errors = [r for r in caplog.records if r.levelname == 'ERROR']
        assert len(errors) == 1
        assert 'line 2' in errors[0].getMessage()
        assert not rpars.fileLoaded['BEAM_EQUIVALENCE']

    def test_already_loaded(self, rpars, in_tmp_dir, mocker):
        """Check that a loaded file is not read again."""
        (in_tmp_dir / BEAM_EQUIVALENCE_NAME).write_text(
            '(1 0) = (1 0)\n', encoding='utf-8',
            )
        rpars.try_loading_beam_equivalence_file()
        reader = mocker.patch(
            'viperleed.calc.classes.rparams.rparams.read_beam_equivalence'
            )
        rpars.try_loading_beam_equivalence_file()
        reader.assert_not_called()
