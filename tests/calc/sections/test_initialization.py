"""Tests for section initialization."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-07-19'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.calc.classes.slab import Slab
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.files import parameters
from viperleed.calc.files import poscar
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.initialization import initialization

from ...helpers import filesystem_to_dict
from ...helpers import raises_test_exception


class TestSetup:
    """Basic tests for pre-running preparation of the work environment."""

    def test_work_path_exists(self, init_files):
        """Check that work_path was created properly."""
        assert init_files.work_path.is_dir()

    def test_files_copied_correctly(self, init_files):
        """Check if all files were copied correctly."""
        input_files = set(f.name for f in init_files.test_path.glob('*'))
        source_copied = all(f in input_files
                            for f in init_files.required_files)
        input_files = set(f.name for f in init_files.work_path.glob('*'))
        work_copied = all(f in input_files for f in init_files.required_files)
        assert source_copied and work_copied


class TestInitialization:                                                       # TODO: find a way to inject the exit-code test (e.g., a class decorator?)
    """Collection of tests for a successful INITIALIZATION run."""

    def test_successful_run(self, init_files):
        """Check that initialization exits without errors."""
        assert not init_files.failed
        assert init_files.records is not None
        assert init_files.records.get_last_state_for_section('initialization')

    _expected_files = 'IVBEAMS', 'BEAMLIST', 'VIBROCC', 'PARAMETERS'

    @parametrize(expected_file=_expected_files)
    def test_init_files_present(self, init_files, expected_file):
        """Ensure the expected files are present after initialization."""
        assert init_files.expected_file_exists(expected_file)

    def test_parameters_was_updated(self, init_files):
        """Check that PARAMETERS file was updated."""
        parameters_path = init_files.work_path / 'PARAMETERS'
        with parameters_path.open('r', encoding='utf-8') as param_file:
            param_content = param_file.read()
        # pylint: disable-next=magic-value-comparison
        assert 'line commented out automatically' in param_content

    def test_does_not_write_out_suffixed(self, init_files):
        """Check that no _OUT-suffixed file is generated."""
        out_suffixed = init_files.work_path.rglob('*_OUT*')
        assert not any(out_suffixed)

    def test_vibrocc_generated(self, init_files):
        """Check that VIBROCC_generated is written if needed."""
        had_vibrocc_input = any(any(p.glob('VIBROCC'))
                                for p in init_files.input_files_paths)
        work = init_files.work_path
        for_supp = work/'VIBROCC_generated'
        for_out = work/'VIBROCC'
        in_supp = work/DEFAULT_SUPP/'VIBROCC_generated'
        in_out = work/DEFAULT_OUT/'VIBROCC'
        in_original_inputs = (
            work/DEFAULT_SUPP/ORIGINAL_INPUTS_DIR_NAME/'VIBROCC'
            )
        if had_vibrocc_input:
            assert not for_supp.exists()
            assert not in_supp.exists()
            assert not in_out.exists()
            assert in_original_inputs.is_file()
        else:
            assert for_supp.is_file()
            assert for_supp.read_text() == for_out.read_text()
            assert in_supp.is_file()
            assert in_out.is_file()
            assert not in_original_inputs.exists()
            # If it has generated a VIBROCC, it also should have
            # modified PARAMETERS, which should now be in OUT
            assert (work/DEFAULT_OUT/'PARAMETERS').is_file()


class TestInitializationDomains:
    """Collection of tests for running the initialization for DOMAINS."""

    def check_path_contains(self, path_tree, expect_contents):
        """Check that all expect_contents are in path_tree."""
        for name, contents in expect_contents.items():
            assert name in path_tree
            if isinstance(contents, str):
                assert path_tree[name] == contents
            else:
                self.check_path_contains(path_tree[name], contents)

    def collect_domain_info(self, init_domains):
        """Collect information about the domains in he calculation."""
        rpars = parameters.read(init_domains.test_path/'PARAMETERS')
        domains = rpars.readParams['DOMAIN']
        src_path_to_folder_created = {
            domain.values_str: f'Domain_{domain.flags_str}'
            for domain in domains
            }
        return src_path_to_folder_created

    @pytest.mark.xfail(reason='Issue #301')
    def test_domain_directories_copied(self, init_domains):                     # TODO: check also contents
        """Check that all Domain_xxx work directories are copied to root."""
        work_domains = tuple(init_domains.work_path.glob('Domain_*'))
        root_domains = tuple(init_domains.test_path.glob('Domain_*'))
        assert len(root_domains) == len(work_domains)

    def test_domain_directories_created(self, init_domains):
        """Check that work contains the right Domain_xxx directories."""
        domains = self.collect_domain_info(init_domains)
        work_domains = tuple(init_domains.work_path.glob('Domain_*'))
        assert len(work_domains) == len(domains)

    def test_successful_run(self, init_domains):
        """Check that initialization exits without errors."""
        TestInitialization.test_successful_run(self, init_domains)

    def test_original_inputs_copied(self, init_domains):
        """Check that all original_inputs were copied correctly."""
        # Prepare the contents that we expect
        test_src = next(p for p in init_domains.input_files_paths
                        # pylint: disable-next=magic-value-comparison
                        if p.name == 'initialization')
        domains = self.collect_domain_info(init_domains)
        original_inputs = filesystem_to_dict(test_src)
        expect_original_inputs = {
            # The ones in root
            DEFAULT_SUPP: {
                ORIGINAL_INPUTS_DIR_NAME: {f: c
                                           for f, c in original_inputs.items()
                                           if isinstance(c, str)},
                },
            }
        # And the ones in each domain folder
        for src_path, domain_folder in domains.items():
            expect_original_inputs[domain_folder] = {
                DEFAULT_SUPP: {
                    ORIGINAL_INPUTS_DIR_NAME: original_inputs[src_path],
                    },
                }
        work_tree = filesystem_to_dict(init_domains.work_path)
        self.check_path_contains(work_tree, expect_original_inputs)
        # TODO: Uncomment the following after #301 is fixed!
        # root_tree = filesystem_to_dict(init_domains.test_path)
        # self.check_path_contains(root_tree, expect_original_inputs)

    def test_phaseshifts_generated(self, init_domains):
        """Check that PHASESHIFTS files were generated in the right places."""
        def _collect_ps_firstline(tree):
            first_line_only = {}
            for folder, files in tree.items():
                try:
                    ps_contents = files['PHASESHIFTS']
                except (KeyError, TypeError):
                    continue
                first, *_ = ps_contents.splitlines()
                # skip the comment at the end, convert to numbers
                *items, _ = first.split()
                first_line_only[folder] = (int(items[0]),
                                           *(float(i) for i in items[1:]))
            return first_line_only

        # Collect the expected file contents. Keeping the first
        # line should be enough to identify different elements.
        expected_ps = _collect_ps_firstline(
            filesystem_to_dict(
                init_domains.test_path/'_expect/work'
                )
            )
        work_ps = _collect_ps_firstline(
            filesystem_to_dict(init_domains.work_path)
            )
        for folder, expect in expected_ps.items():
            assert work_ps[folder] == pytest.approx(expect)


class TestInitializationRaises:
    """Tests for checking exceptions raised during initialization."""

    @pytest.fixture(name='ag100_init')
    def fixture_ag100_init(self, ag100, make_section_tempdir, tensorleed_path):
        """Yield slab and rpars ready to execute in a temporary directory."""
        slab, rpars, *_ = ag100
        rpars.paths.tensorleed = tensorleed_path
        tmp = make_section_tempdir('Ag(100)', 'init')
        with execute_in_dir(tmp):  # Not to spam with files
            yield slab, rpars

    def test_bulk_appended_raises(self, ag100_init, caplog):
        """Ensure limited exceptions are caught when appending bulk units."""
        with raises_test_exception(Slab, 'with_extra_bulk_units'):
            initialization(*ag100_init)
        # Notice that here we should also check the logger, as we
        # do call with_extra_bulk_units twice during initialization:
        # once when writing the POSCAR_bulk_appended (which used to
        # be swallowed) and once when generating PHASESHIFTS. Not
        # checking the logging messages would make this test succeed
        # because the exception is raised while making PHASESHIFTS.
        # pylint: disable=magic-value-comparison
        messages = (m for m in caplog.messages if 'exception' in m.lower())
        messages = (m for m in messages if 'bulk_appended' in m)
        log_bulk_appended = next(messages, '')
        assert not log_bulk_appended

    def test_poscar_write_raises(self, ag100_init):
        """Ensure limited exceptions are caught when writing a POSCAR."""
        with raises_test_exception(poscar, 'write'):
            initialization(*ag100_init)
