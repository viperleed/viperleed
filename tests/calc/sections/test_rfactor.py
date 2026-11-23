"""Tests for the choice of R factor in section Rfactor."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-20'
__license__ = 'GPLv3+'

from pathlib import Path
from types import SimpleNamespace

import pytest
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.rparams.special.energy_range import TheoEnergies
from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.files import parameters
from viperleed.calc.lib.version import Version
from viperleed.calc.sections import rfactor as rfactor_section

_RUNNERS = {
    'run_legacy_rfactor': 'legacy',
    'run_new_rfactor': 'python',
    }
_LEGACY_CASES = (
    pytest.param('', 'legacy', id='not given'),
    pytest.param('R_FACTOR_LEGACY = True\n', 'legacy', id='true'),
    pytest.param('R_FACTOR_LEGACY = False\n', 'python', id='false'),
    )


def test_r_factor_legacy_default():
    """Check that the TensErLEED R factor is the default."""
    assert Rparams().R_FACTOR_LEGACY is True


@parametrize('line,expected', _LEGACY_CASES)
def test_rfactor_dispatch(line, expected, monkeypatch, tmp_path):
    """Check which R factor is run for PARAMETERS with `line`."""
    called_runners = []
    for runner, label in _RUNNERS.items():
        def _record(*_, label=label, **__):
            called_runners.append(label)
        monkeypatch.setattr(rfactor_section, runner, _record)
    parameters_file = tmp_path / 'PARAMETERS'
    parameters_file.write_text(f'RUN = 0 1\n{line}', encoding='utf-8')
    rpars = parameters.read(parameters_file)
    parameters.interpret(rpars)
    rpars.THEO_ENERGIES = TheoEnergies.from_sorted_grid((50.0, 60.0, 70.0))
    rpars.refcalc_fdout = 'not empty, no need to fetch spectra'
    rfactor_section.rfactor(None, rpars, index=11)
    assert called_runners == [expected]


# R_FACTOR_TYPE, TL_VERSION, and the R factor TensErLEED computes
_LEGACY_COMPUTED = (
    pytest.param('pendry', '2.1.0', 'pendry', id='pendry'),
    pytest.param('smooth', '2.1.0', 'smooth', id='smooth'),
    pytest.param('smooth', '2.0.0', 'pendry', id='smooth-tl2.0'),
    pytest.param('r2', '2.0.0', 'r2', id='r2-tl2.0'),
    pytest.param('zj', '2.0.0', 'zj', id='zj-tl2.0'),
    )


def _fake_tensorleed_run(rpars, tmp_path, monkeypatch):
    """Replace the TensErLEED R factor with one that writes ROUT."""
    section = rfactor_section
    tl_path = tmp_path / 'tensorleed'
    for sub_dir, f_name in (('lib', 'rfacsb.f'), ('src', 'rfactor.f')):
        (tl_path / sub_dir).mkdir(parents=True)
        (tl_path / sub_dir / f_name).touch()
    rpars.get_tenserleed_directory = lambda: SimpleNamespace(path=tl_path)
    rpars.TL_IGNORE_CHECKSUM = True
    rpars.FORTRAN_COMP = ['gfortran', '']
    for name in ('writeWEXPEL', 'writeRfactPARAM'):
        monkeypatch.setattr(section.iorfactor, name, lambda *_, **__: None)
    monkeypatch.setattr(section.leedbase,
                        'fortran_compile_batch',
                        lambda *_, **__: None)
    monkeypatch.setattr(section.fs_utils, 'move', lambda *_: None)
    monkeypatch.setattr(section.subprocess,
                        'run',
                        lambda *_, **__: Path('ROUT').touch())
    monkeypatch.setattr(section.iorfactor,
                        'readROUT',
                        lambda: ((0.3, 0.3, 0.), 0., [0.3]))


@parametrize('r_factor_type,tl_version,expected', _LEGACY_COMPUTED)
def test_legacy_analysis_r_factor_type(r_factor_type, tl_version, expected,
                                       tmp_path, monkeypatch):
    """Check the R factor of the analysis PDF of the TensErLEED R."""
    monkeypatch.chdir(tmp_path)
    rpars = Rparams()
    rpars.R_FACTOR_TYPE = RFactorType(r_factor_type)
    rpars.TL_VERSION = Version(tl_version)
    rpars.expbeams = [SimpleNamespace(getLabel=lambda **_: ('(1|0)', 5))]
    _fake_tensorleed_run(rpars, tmp_path, monkeypatch)
    pdf_calls = []
    monkeypatch.setattr(rfactor_section.iorfactor,
                        'writeRfactorPdf',
                        lambda *_, **kwargs: pdf_calls.append(kwargs))
    rfactor_section.run_legacy_rfactor(None, rpars,
                                       for_error=False,
                                       name='refcalc',
                                       theobeams=[],
                                       index=11,
                                       only_vary=None)
    assert rpars.last_R == pytest.approx(0.3)
    assert [c['which_r'] for c in pdf_calls] == [RFactorType(expected)]
