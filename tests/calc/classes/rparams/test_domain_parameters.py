"""Tests for the domain_params module of viperleed.calc.classes.rparams."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-19'
__license__ = 'GPLv3+'

from pathlib import Path

from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.domain_params import DomainParameters


class TestDomainParameters:
    """Tests for the DomainParameters class."""

    @fixture(name='make_domain')
    def factory_domain(self):
        """Return an initialized DomainParameters."""
        def _make(*args):
            return DomainParameters(*args)
        return _make

    def test_init(self, make_domain):
        """Test initialization of DomainParameters."""
        work, name = 'test_workdir', 'test_domain'
        domain = make_domain(work, name)

        assert domain.workdir == Path(work).resolve()
        assert domain.name == name
        assert domain.sl is None
        assert domain.rp is None
        assert domain.refcalcRequired is False
        assert domain.tensorDir is None

    def test_invalid_path(self, make_domain):
        """Test with a non-existent path to ensure it still resolves."""
        work = 'non_existent_path'
        domain = make_domain(work, 'test_domain')
        assert domain.workdir == Path(work)

    @parametrize(name=('test_domain', ''))
    def test_str(self, name, make_domain):
        """Test the string representation of DomainParameters."""
        domain = make_domain('test_workdir', name)
        assert str(domain) == f'domain {name}'

    def test_path_already_resolved(self, make_domain):
        """Test that the workdir path is correctly resolved."""
        work = Path.cwd()
        domain = make_domain(work, 'test_domain')
        assert domain.workdir == work

    def test_path_resolution(self, make_domain):
        """Test that the workdir path is correctly resolved."""
        domain = make_domain('.', 'test_domain')
        assert domain.workdir == Path.cwd()
