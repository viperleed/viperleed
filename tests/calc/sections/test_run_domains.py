"""Tests for an entire multi-domain calculation."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-10-22'
__license__ = 'GPLv3+'


from pytest_cases import fixture

from .conftest import DomainsCalcFilesSetup as Setup
from .conftest import with_tl_versions


@fixture(name='domains_full', scope='session')
@with_tl_versions
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
