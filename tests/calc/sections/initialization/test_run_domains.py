"""Tests for a multi-domain execution of section initialization."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-07'
__license__ = 'GPLv3+'

import pytest

from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import LOG_PREFIX
from viperleed.calc.constants import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.files import parameters
from viperleed.calc.lib.string_utils import strip_comments

from ....helpers import filesystem_to_dict
from .test_run_one_domain import TestInitialization as _TestBasic


class TestInitializationDomains:
    """Collection of tests for running the initialization for DOMAINS."""

    def check_path_contains(self, path_tree, expect_contents):
        """Check that all `expect_contents` are in `path_tree`."""
        for name, contents in expect_contents.items():
            assert name in path_tree
            if isinstance(contents, str):
                assert path_tree[name] == contents
            else:
                self.check_path_contains(path_tree[name], contents)

    def test_domain_directories_created(self, init_domains):
        """Check that work contains the right domain directories."""
        assert len(init_domains.src_folders) == len(init_domains.work_domains)

    def test_domain_output_copied(self, init_domains):
        """Check that the results of the main calculation are copied back."""
        # NB: Not actually copied back, since we only run_calc
        manifest = init_domains.read_manifest()
        for domain_path in init_domains.work_domains.values():
            assert f'{domain_path.name}/{DEFAULT_SUPP}' in manifest
            assert f'{domain_path.name}/{DEFAULT_OUT}' in manifest

    def test_domain_paths_edited(self, init_domains):
        """Check that the PARAMETERS file was updated with the new paths."""
        out_path = init_domains.work_path/DEFAULT_OUT
        param_out = out_path/'PARAMETERS'
        edited_assignments = parameters.read(param_out).readParams['DOMAIN']
        edited = {a.values_str for a in edited_assignments}
        work_paths = (f'./{p.name}' for p in init_domains.work_domains.values())
        assert all(new_path in edited for new_path in work_paths)

    def test_ivbeams_generated(self, init_domains):
        """Check that domain folders contain copies of the main IVBEAMS."""
        main_ivbeams = (init_domains.work_path/'IVBEAMS').read_text()
        for domain in init_domains.work_domains.values():
            if (domain/DEFAULT_TENSORS).exists():
                continue   # IVBEAMS copied from tensor
            domain_ivbeams = (domain/'IVBEAMS').read_text()
            assert domain_ivbeams == main_ivbeams

    def test_main_output_copied(self, init_domains, re_match):
        """Check that the results of the main calculation are copied back."""
        # NB: Not actually copied back, since we only run_calc
        manifest = init_domains.read_manifest()
        assert any(re_match(rf'{LOG_PREFIX}', line) for line in manifest)
        assert DEFAULT_SUPP in manifest
        assert DEFAULT_OUT in manifest

    def test_original_inputs_copied(self, init_domains):
        """Check that all original_inputs were copied correctly."""
        # Prepare the contents that we expect
        test_src = next(p for p in init_domains.input_files_paths
                        # pylint: disable-next=magic-value-comparison
                        if p.name == 'initialization')
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
        for src_path, domain_folder in init_domains.work_domains.items():
            expect_original_inputs[domain_folder.name] = {
                DEFAULT_SUPP: {
                    ORIGINAL_INPUTS_DIR_NAME: original_inputs[src_path],
                    },
                }
        work_tree = filesystem_to_dict(init_domains.work_path)
        self.check_path_contains(work_tree, expect_original_inputs)

    def test_parameters_inherited(self, init_domains):
        """Ensure that domain PARAMETERS are updated from the main ones."""
        updated_params = (
            'THEO_ENERGIES',
            'BEAM_INCIDENCE',
            # There would also be LMAX for TL_VERSION <= 1.6.0
            )
        for domain in init_domains.work_domains.values():
            params_lines = [
                line for line in (domain/'PARAMETERS').read_text().splitlines()
                if strip_comments(line)
                ]
            for updated in updated_params:
                assert any(line for line in params_lines if updated in line)

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

    def test_successful_run(self, init_domains):
        """Check that initialization exits without errors."""
        _TestBasic.test_successful_run(self, init_domains)
