"""Test configuration for calc/classes/rparams/domain_params."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-20'
__license__ = 'GPLv3+'

from pytest_cases import fixture

from viperleed.calc.classes.rparams.domain_params import DomainParameters


@fixture(name='make_domain')
def factory_domain():
    """Return an initialized DomainParameters."""
    def _make(*args):
        return DomainParameters(*args)
    return _make
