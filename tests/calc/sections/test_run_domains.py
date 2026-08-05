"""Tests for an entire multi-domain calculation."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-10-22'
__license__ = 'GPLv3+'


from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.files import tenserleed
from viperleed.calc.lib.version import Version

from .conftest import DomainsCalcFilesSetup as Setup

@fixture(name='domains_full', scope='session')
@parametrize(tl_version=(
    tenserleed.CURRENT_TL_VERSION,
    Version('1.7.5'),    # The first one where search works
    # Earlier versions fail during search for various reasons:
    # - 1.7.1 -- 1.7.4
    #    Race condition in writing SD.TL in fortran code prevents
    #    from reading a full configuration. Fixed in v1.7.5.
    # - 1.6.1
    #    With the IVBEAMS currently used, at least one beam has
    #    too few data for interpolation in the given energy range.
    #    This beam, however, is used in lib.search/COMNEI (though
    #    it should be skipped), giving a SEGFAULT (index error).
    #    Using a smaller set of IVBEAMS (and thus re-running a
    #    refcalc) prevents the segfault but has the same race
    #    condition problem as v1.7.1--1.7.4.
    ))
def fixure_domains_full(tl_version, tensorleed_path, make_section_tempdir):
    """A setup for a full multi-domain calculation."""
    domains = 'telluride_stacking'
    setup = Setup(
        surface_dir=domains,
        tmp_test_path=make_section_tempdir(domains, 'full'),
        copy_dirs=['deltas', 'full_run'],
        )
    setup.run_calc_from_setup(
        source=tensorleed_path,
        preset_params={'RUN': [0, 4],  # Initialization & domains
                       'TL_VERSION': tl_version,}
        )
    return setup


def test_success(domains_full):
    """Check that the calculation ended without errors."""
    assert not domains_full.failed
