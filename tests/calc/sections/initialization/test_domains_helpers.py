"""Tests for helper functions of a multi-domain initialization."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-07'
__license__ = 'GPLv3+'

from copy import deepcopy

import pytest
from pytest_cases import fixture

from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.version import Version
from viperleed.calc.sections.initialization import _check_domain_consistent
from viperleed.calc.sections.initialization import _inherit_from_main
from viperleed.calc.sections.initialization import _make_domain_workdir

_MODULE = 'viperleed.calc.sections.initialization'


class TestCheckDomainConsistent:
    """Tests for the _check_domain_consistent domain helper function."""

    @fixture(name='rpars')
    def fixture_mock_rpars(self, mocker):
        """Create a mock Rparams-like object."""
        return mocker.MagicMock(
            LMAX=mocker.MagicMock(max=5),
            THEO_ENERGIES=mocker.MagicMock(),
            THETA=1.0,
            PHI=2.0,
            TL_VERSION=Version('1.6.0'),
            ivbeams=[],
            fileLoaded={'IVBEAMS': True}
            )

    @fixture(name='domain')
    def fixture_mock_domain(self, mocker, rpars):
        """Create a mock DomainParameters-like object."""
        return mocker.MagicMock(
            rpars=deepcopy(rpars),
            refcalc_required=False,
            )

    def test_all_consistent(self, domain, rpars, mocker, caplog):
        """Check expected behavior if everything is consistent."""
        domain.rpars.THEO_ENERGIES.contains.return_value = True
        dom_beam, main_beam = mocker.MagicMock(), mocker.MagicMock()
        dom_beam.isEqual.return_value = True
        domain.rpars.ivbeams = [dom_beam]
        rpars.ivbeams = [main_beam]
        with caplog.at_level('INFO'):
            result = _check_domain_consistent(domain, rpars)
        assert not result
        assert not domain.refcalc_required
        assert not caplog.records

    def test_inconsistent_beam_incidence(self, domain, rpars, caplog):
        """Check expected behavior when BEAM_INCIDENCE is mismatched."""
        domain.rpars.THEO_ENERGIES.contains.return_value = True
        rpars.TL_VERSION = Version('2.0.0')
        domain.rpars.THETA = 99
        with caplog.at_level('INFO'):
            result = _check_domain_consistent(domain, rpars)
        assert result == ['BEAM_INCIDENCE']
        assert domain.refcalc_required
        expect_log = 'BEAM_INCIDENCE is mismatched'
        assert expect_log in caplog.text

    def test_inconsistent_lmax(self, domain, rpars, caplog):
        """Check expected behavior with mismatched LMAX for early versions."""
        domain.rpars.THEO_ENERGIES.contains.return_value = True
        rpars.TL_VERSION = Version('1.6.0')
        domain.rpars.LMAX.max = 7
        with caplog.at_level('INFO'):
            result = _check_domain_consistent(domain, rpars)
        assert result == ['LMAX']
        assert domain.refcalc_required
        expect_log = 'LMAX is mismatched'
        assert expect_log in caplog.text

    def test_inconsistent_theo_energies(self, domain, rpars, caplog):
        """Check expected behavior when THEO_ENERGIES is mismatched."""
        domain.rpars.THEO_ENERGIES.contains.return_value = False
        with caplog.at_level('INFO'):
            result = _check_domain_consistent(domain, rpars)
        assert result == ['THEO_ENERGIES']
        assert domain.refcalc_required
        expect_log = 'Energy range is mismatched'
        assert expect_log in caplog.text

    def test_mismatched_ivbeams_content(self, domain, rpars, mocker, caplog):
        """Check expected behavior when IVBEAMS beams are different."""
        domain.rpars.THEO_ENERGIES.contains.return_value = True
        dom_beam, main_beam = mocker.MagicMock(), mocker.MagicMock()
        dom_beam.isEqual.return_value = False
        domain.rpars.ivbeams = [dom_beam]
        rpars.ivbeams = [main_beam]
        with caplog.at_level('INFO'):
            result = _check_domain_consistent(domain, rpars)
        assert not result
        assert domain.refcalc_required
        expect_log = 'IVBEAMS file mismatched'
        assert expect_log in caplog.text

    def test_mismatched_ivbeams_length(self, domain, rpars, caplog):
        """Check expected behavior for mismatched number of IVBEAMS."""
        domain.rpars.THEO_ENERGIES.contains.return_value = True
        domain.rpars.ivbeams = [1]
        rpars.ivbeams = [1, 2]
        with caplog.at_level('INFO'):
            result = _check_domain_consistent(domain, rpars)
        assert not result
        assert domain.refcalc_required
        expect_log = 'IVBEAMS file mismatched'
        assert expect_log in caplog.text

    def test_missing_ivbeams_file(self, domain, rpars, caplog):
        """Check expected behavior when no IVBEAMS was read."""
        domain.rpars.THEO_ENERGIES.contains.return_value = True
        domain.rpars.fileLoaded['IVBEAMS'] = False
        with caplog.at_level('INFO'):
            result = _check_domain_consistent(domain, rpars)
        assert not result
        assert domain.refcalc_required
        expect_log = 'No IVBEAMS file loaded'
        assert expect_log in caplog.text

    def test_refcalc_already_required(self, domain, rpars, caplog):
        """Check shortcut if domain already requires a refcalc."""
        domain.refcalc_required = True
        with caplog.at_level(0):  # All messages
            result = _check_domain_consistent(domain, rpars)
        assert not result
        assert not caplog.records


class TestInheritFromMain:
    """Tests for the _inherit_from_main domain helper function."""

    expect_inherited = (
        'THEO_ENERGIES',
        'THETA',
        'PHI',
        'N_CORES',
        'ivbeams',
        )

    @fixture(name='main_rpars')
    def mock_main_rpars(self, mocker):
        """Return a fake Rparams-like object for the main calculation."""
        mock = mocker.MagicMock()
        mock.TL_VERSION = Version('1.6.0')
        mock.LMAX.max = 7
        return mock

    @fixture(name='make_domain')
    def fixture_make_domain(self, mocker, tmp_path):
        """Return a factory for domain mocks with a given workdir."""
        def _make(name='d1'):
            domain = mocker.MagicMock()
            domain.workdir = tmp_path / name
            domain.workdir.mkdir()
            domain.rpars = mocker.MagicMock()
            domain.rpars.readParams = {'THEO_ENERGIES'}
            domain.rpars.LMAX.max = 0
            return domain
        return _make

    @fixture(name='mocks')
    def fixture_mocks(self, mocker):
        """Replace implementation details with mocks."""
        mock_execute = mocker.patch(f'{_MODULE}.execute_in_dir')
        mock_execute.return_value.__enter__.return_value = None
        mock_execute.return_value.__exit__.return_value = None
        return {
            'copy': mocker.patch('shutil.copy2'),
            'modify': mocker.patch(f'{_MODULE}.parameters.modify'),
            'execute_in_dir': mock_execute,
            }

    def test_inherit_lmax(self, mocks, mocker, make_domain, main_rpars):
        """Check inheritance of LMAX.max for early TL_VERSION."""
        domain = make_domain('domain')

        inconsistent = {domain: ['THEO_ENERGIES']}
        _inherit_from_main([domain], main_rpars, inconsistent)

        # inherit_from called with expected arguments
        domain.rpars.inherit_from.assert_called_once()
        args, _ = domain.rpars.inherit_from.call_args
        assert all(inherited in args for inherited in self.expect_inherited)

        # copy2 called once per domain
        mocks['copy'].assert_called_once_with('IVBEAMS', domain.workdir)

        # LMAX updated for TL_VERSION <= 1.6.0
        assert domain.rpars.LMAX.max == main_rpars.LMAX.max

        # PARAMETERS modification, done in the work directories
        mocks['execute_in_dir'].assert_called_once_with(domain.workdir)
        expect_write = (
            'THEO_ENERGIES',   # Inconsistent
            'BEAM_INCIDENCE',  # Not explicitly given
            'LMAX',            # Early version
            )
        expect_calls = [mocker.call(domain.rpars, written)
                        for written in expect_write]
        assert len(mocks['modify'].mock_calls) == len(expect_write)
        mocks['modify'].assert_has_calls(expect_calls, any_order=True)

    def test_no_inheritance(self, mocks, make_domain, main_rpars):
        """Check that no parameters are inherited if all are consistent."""
        domain = make_domain('domain')
        main_rpars.TL_VERSION = Version('2.0.0')

        inconsistent = {domain: []}
        # No inconsistent params, all readParams filled
        domain.rpars.readParams = {'THEO_ENERGIES', 'BEAM_INCIDENCE', 'LMAX'}

        _inherit_from_main([domain], main_rpars, inconsistent)

        # LMAX should NOT be updated (since TL_VERSION > 1.6.0)
        assert domain.rpars.LMAX.max != main_rpars.LMAX.max

        # IVBEAMS copied
        mocks['copy'].assert_called_once_with('IVBEAMS', domain.workdir)
        # But no parameters modified
        mocks['modify'].assert_not_called()  # All consistent
        mocks['execute_in_dir'].assert_called_once_with(domain.workdir)

    def test_multiple_domains(self, mocks, mocker, make_domain, main_rpars):
        """Check behavior when multiple domains are given."""
        domain_one = make_domain('domain_one')
        domain_two = make_domain('domain_two')

        inconsistent = {
            domain_one: ['THEO_ENERGIES'],
            domain_two: ['BEAM_INCIDENCE'],
        }

        domains = [domain_one, domain_two]
        for domain in domains:
            domain.rpars.readParams = {'THEO_ENERGIES',
                                       'BEAM_INCIDENCE',
                                       'LMAX'}
        _inherit_from_main(domains, main_rpars, inconsistent)

        # Each domain's inherit_from called separately
        for domain in domains:
            domain.rpars.inherit_from.assert_called_once()
        assert mocks['copy'].call_count == len(domains)
        assert mocks['execute_in_dir'].call_count == len(domains)

        expect_calls = [
            # Only the inconsistent once, as all are simulated as read
            mocker.call(domain_one.rpars, 'THEO_ENERGIES'),
            mocker.call(domain_two.rpars, 'BEAM_INCIDENCE'),
            ]
        mocks['modify'].assert_has_calls(expect_calls, any_order=True)


class TestMakeDomainWork:
    """Tests for the _make_domain_workdir function."""

    @fixture(name='make_work')
    def fixture_make_workdir(self, tmp_path):
        """Call _make_domain_workdir in a subfolder of tmp_path."""
        main_work = tmp_path/'main_work'
        main_work.mkdir()
        def _make(**kwargs):
            kwargs.setdefault('src', None)
            kwargs.setdefault('calc_started_at', None)
            kwargs.setdefault('must_use_auto_name', False)
            with execute_in_dir(main_work):
                return _make_domain_workdir(**kwargs)
        return _make

    def test_already_exists(self, make_work, tmp_path, caplog):
        """Check successful execution when the work path exists already."""
        name = 'some_domain'
        expect_work = tmp_path / f'main_work/Domain_{name}'
        expect_work.mkdir()
        make_work(name=name)
        # pylint: disable-next=magic-value-comparison
        assert 'already exists' in caplog.text

    def test_mkdir_fails(self, make_work, mocker):
        """Check exceptions raised by mkdir are not caught."""
        mock_mkdir = mocker.patch('pathlib.Path.mkdir', side_effect=OSError)
        with pytest.raises(OSError):
            make_work(name='fails')
        mock_mkdir.assert_called_once()

    def test_read_from_zip(self, make_work, tmp_path):
        """Check outcome when the domain was read from a tensor ZIP."""
        domain_src = tmp_path/'domain_source.zip'
        domain_src.touch()
        name = 'some_domain'
        work = make_work(name=name, src=domain_src, calc_started_at=tmp_path)
        assert work.is_dir()
        assert work.name == f'Domain_{name}'

    def test_read_from_calc_subfolder(self, make_work, tmp_path):
        """Check outcome when the domain source is a subfolder of calc."""
        domain_src = tmp_path/'domain_source'
        domain_src.mkdir()
        name = 'some_domain'
        work = make_work(name=name, src=domain_src, calc_started_at=tmp_path)
        assert work.is_dir()
        assert work.name == domain_src.name

    def test_read_from_calc_subfolder_but_cant_use(self, make_work, tmp_path):
        """Check usage of auto-generated name if explicitly requested."""
        domain_src = tmp_path/'domain_source'
        domain_src.mkdir()
        name = 'some_domain'
        work = make_work(name=name,
                         src=domain_src,
                         calc_started_at=tmp_path,
                         must_use_auto_name=True)
        assert work.is_dir()
        assert work.name == f'Domain_{name}'

    def test_read_from_somewhere_else(self, make_work, tmp_path):
        """Check outcome when the domain source is not a subfolder of calc."""
        domain_src = tmp_path/'some_other_folder'/'domain_source'
        domain_src.mkdir(parents=True)
        name = 'some_domain'
        work = make_work(name=name, src=domain_src, calc_started_at=tmp_path)
        assert work.is_dir()
        assert work.name == f'Domain_{name}'

    def test_success(self, make_work):
        """Check successful creation of a domain work directory."""
        name = 'some_domain'
        work = make_work(name=name)
        assert work.is_dir()
        assert work.name == f'Domain_{name}'
