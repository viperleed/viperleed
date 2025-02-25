"""Tests for _propagate_to_domains of viperleed.calc.sections.cleanup."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-25'
__license__ = 'GPLv3+'

from pathlib import Path

from pytest_cases import fixture

from viperleed.calc.classes.rparams.domain_params import DomainParameters
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.cleanup import _propagate_to_domains


@fixture(name='make_domain')
def factory_domain(tmp_path):
    """Return a DomainParameters instance."""
    def _make(name):
        domain = DomainParameters(tmp_path/name, name)
        domain.rp = Rparams()
        domain.workdir.mkdir()
        return domain
    return _make


@fixture(name='make_nested_domains')
def factory_nested_domains(make_domain):
    """Return an Rparams with multiple nested domains."""
    def _make(root, nested_domains):
        for name, subdomains in nested_domains.items():
            domain = make_domain(name)
            _make(domain.rp, {f'{name}/{d}': s for d, s in subdomains.items()})
            root.domainParams.append(domain)
    return _make


@fixture(name='track_cwd')
def fixture_track_cwd(tmp_path):
    """Return a callable that remembers its executions directories."""
    @_propagate_to_domains
    def _track(rpars, cache):
        cwd = Path.cwd().relative_to(tmp_path)
        cache[rpars] = str(cwd.as_posix())

    def _call(rpars):
        cache = {}
        with execute_in_dir(tmp_path):
            _track(rpars, cache)
        return cache
    return _call


def test_args_kwargs(rpars, make_nested_domains, mocker):
    """Test that decorated function receives args and kwargs correctly."""
    nested = {'one': {'two': {'three': {}}}}
    make_nested_domains(rpars, nested)
    func_mock = mocker.MagicMock()

    @_propagate_to_domains
    def decorated_func(rpars_, arg, key=None):
        func_mock(rpars_, arg, key)

    decorated_func(rpars, 'test_arg', key='test_value')
    one = rpars.domainParams[0]
    two = one.rp.domainParams[0]
    three = two.rp.domainParams[0]

    for rpars_arg in (rpars, one.rp, two.rp, three.rp):
        func_mock.assert_any_call(rpars_arg, 'test_arg', 'test_value')


def test_nested_domains(rpars, track_cwd, make_nested_domains):
    """Test _propagate_to_domains with nested domains."""
    nested = {'one': {'two': {'three': {}}}}
    make_nested_domains(rpars, nested)
    workdirs = track_cwd(rpars)
    one = rpars.domainParams[0]
    two = one.rp.domainParams[0]
    three = two.rp.domainParams[0]
    expect_workdirs = {rpars: '.',
                       one.rp: 'one',
                       two.rp: 'one/two',
                       three.rp: 'one/two/three'}
    assert workdirs == expect_workdirs


def test_no_domain(rpars, track_cwd):
    """Test _propagate_to_domains with no domains."""
    workdirs = track_cwd(rpars)
    assert workdirs == {rpars: '.'}


def test_simple_domains(rpars, track_cwd, make_nested_domains):
    """Test _propagate_to_domains with no nested domains."""
    flat = {'one': {}, 'two': {}}
    make_nested_domains(rpars, flat)
    domains = rpars.domainParams
    workdirs = track_cwd(rpars)
    expect_workdirs = {rpars: '.',
                       domains[0].rp: 'one',
                       domains[1].rp: 'two'}
    assert workdirs == expect_workdirs
