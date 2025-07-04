"""Tests for a single-domain execution of section initialization."""

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
from viperleed.calc.files import poscar
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.initialization import initialization

from ....helpers import raises_test_exception


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
