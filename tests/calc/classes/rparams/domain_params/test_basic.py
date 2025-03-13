"""Tests module domain_params of viperleed.calc.classes.rparams.

This module contains tests for the basic implementation details of the
DomainParameters class, especially magic methods. 
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-19'
__license__ = 'GPLv3+'

from pathlib import Path

from pytest_cases import parametrize


class TestDomainParameters:
    """Tests for the DomainParameters class."""

    def test_init(self, make_domain):
        """Test initialization of DomainParameters."""
        work, name = 'test_workdir', 'test_domain'
        domain = make_domain(work, name)

        assert domain.workdir == Path(work).resolve()
        assert domain.name == name
        assert domain.slab is None
        assert domain.rpars is None
        assert domain.refcalc_required is False

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
