"""Tests for a multi-domain execution of section initialization."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-07'
__license__ = 'GPLv3+'

from pathlib import Path

import pytest

from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import LOG_PREFIX
from viperleed.calc.constants import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.files import parameters

from ....helpers import filesystem_to_dict
from .test_run_one_domain import TestInitialization as _TestBasic


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
        """Collect information about the domains in the calculation."""
        test_src = next(p for p in init_domains.input_files_paths
                        # pylint: disable-next=magic-value-comparison
                        if p.name == 'initialization')
        rpars = parameters.read(test_src/'PARAMETERS')
        domains = rpars.readParams['DOMAIN']
        src_to_folders_created = {
            domain.values_str: {f'Domain_{domain.flags_str}',
                                f'{Path(domain.values_str).name}'}
            for domain in domains
            }
        work_domains = {
            src: init_domains.work_path/n
            for src, paths in src_to_folders_created.items() for n in paths
            if (init_domains.work_path/n).is_dir()
            }
        return src_to_folders_created, work_domains

    @pytest.mark.xfail(reason='Issue #301')
    def test_domain_directories_copied(self, init_domains):
        """Check that all Domain_xxx work directories are copied to root."""
        _, work_domains = self.collect_domain_info(init_domains)
        assert all((init_domains.test_path/p.name).is_dir()
                   for p in work_domains.values())
        for work_path in work_domains.values():
            work_tree = filesystem_to_dict(work_path)
            root_tree = filesystem_to_dict(
                init_domains.test_path/work_path.name
                )
            assert root_tree[DEFAULT_SUPP] == work_tree[DEFAULT_SUPP]
            assert root_tree[DEFAULT_OUT] == work_tree[DEFAULT_OUT]

    def test_domain_directories_created(self, init_domains):
        """Check that work contains the right Domain_xxx directories."""
        src_domains, work_domains = self.collect_domain_info(init_domains)
        assert len(work_domains) == len(src_domains)

    @pytest.mark.xfail(reason='Issue #301')
    def test_dommain_output_copied(self, init_domains):
        """Check that the results of the main calculation are copied back."""
        # NB: Not actually copied back, since we only run_calc
        manifest = init_domains.read_manifest()
        _, domains = self.collect_domain_info(init_domains)
        for domain_path in domains.values():
            assert f'{domain_path.name}/{DEFAULT_SUPP}' in manifest
            assert f'{domain_path.name}/{DEFAULT_OUT}' in manifest

    def test_domain_paths_edited(self, init_domains):
        """Check that the PARAMETERS file was updated with the new paths."""
        out_path = init_domains.work_path/DEFAULT_OUT
        edited_assignments = parameters.read(out_path/'PARAMETERS').readParams['DOMAIN']
        edited = {a.values_str for a in edited_assignments}
        _, work_domains = self.collect_domain_info(init_domains)
        work_paths = (f'./{p.name}' for p in work_domains.values())
        assert all(new_path in edited for new_path in work_paths)

    def test_main_output_copied(self, init_domains, re_match):
        """Check that the results of the main calculation are copied back."""
        # NB: Not actually copied back, since we only run_calc
        manifest = init_domains.read_manifest()
        assert any(re_match(rf'{LOG_PREFIX}', line) for line in manifest)
        assert DEFAULT_SUPP in manifest
        assert DEFAULT_OUT in manifest

    def test_successful_run(self, init_domains):
        """Check that initialization exits without errors."""
        _TestBasic.test_successful_run(self, init_domains)

    def test_original_inputs_copied(self, init_domains):
        """Check that all original_inputs were copied correctly."""
        # Prepare the contents that we expect
        test_src = next(p for p in init_domains.input_files_paths
                        # pylint: disable-next=magic-value-comparison
                        if p.name == 'initialization')
        _, domains = self.collect_domain_info(init_domains)
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
            expect_original_inputs[domain_folder.name] = {
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
