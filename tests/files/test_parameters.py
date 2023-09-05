"""Test module parameters of viperleed.tests.

Created on 2023-06-09

@author: Alexander M. Imre

"""

import pytest
from pathlib import Path
import os, sys
from copy import deepcopy
import numpy as np

vpr_path = str(Path(__file__).parent.parent.parent)
if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))

import viperleed.tleedmlib.files.parameters as parameters
from viperleed.tleedmlib.files.parameters import (readPARAMETERS,
                                                  ParameterInterpreter,
                                                  Assignment, NumericBounds)
from viperleed.tleedmlib.files import poscar
from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tleedmlib.files.parameter_errors import (
    ParameterError, ParameterValueError, ParameterParseError,
    ParameterIntConversionError, ParameterFloatConversionError,
    ParameterBooleanConversionError, ParameterNotRecognizedError,
    ParameterNumberOfInputsError, ParameterRangeError,
    ParameterUnknownFlagError, ParameterNeedsFlagError
    )

_FIXTURES_PATH = Path('tests/fixtures/')                                        # TODO: use conftest functionality?


@pytest.fixture()
def ag100_parameters_example():
    # read Ag(100) POSCAR and PARAMETERS files
    slab = poscar.read(_FIXTURES_PATH / 'Ag(100)' / 'initialization' / 'POSCAR')
    rpars = readPARAMETERS(_FIXTURES_PATH / 'Ag(100)' / 'initialization' / 'PARAMETERS')
    # interpret PARAMETERS file
    interpreter = ParameterInterpreter(rpars)
    interpreter.interpret(slab)
    return (rpars, slab)


@pytest.fixture(scope='function')
def slab_ag100():
    # read Ag(100) POSCAR
    return poscar.read(_FIXTURES_PATH / 'POSCARs' / 'POSCAR_Ag(100)')


@pytest.fixture()
def slab_ir100_2x1_o():
    # read Ir(100)-(2x1)-O POSCAR
    return poscar.read(_FIXTURES_PATH / 'POSCARs' / 'POSCAR_Ir(100)-(2x1)-O')


@pytest.fixture()
def ir100_2x1_o_parameters_example(slab_ir100_2x1_o):
    slab = slab_ir100_2x1_o
    rpars = readPARAMETERS(_FIXTURES_PATH / 'parameters' / 'PARAMETERS_Ir(100)-(2x1)-O')
    interpreter = ParameterInterpreter(rpars)
    interpreter.interpret(slab)
    return (rpars, slab)


@pytest.fixture(scope='function')
def mock_rparams():
    return Rparams()


class TestAg100Parameters():
    def test_read_parameters_for_ag100(self):
        # just check that readPARAMETERS does not crash; not interpreted yet
        filename = 'tests/fixtures/Ag(100)/initialization/PARAMETERS'
        rpars = readPARAMETERS(filename)
        assert rpars

    def test_interpret_parameters_for_ag100(self, ag100_parameters_example):
        rpars, slab = ag100_parameters_example
        # check that the parameters are interpreted correctly based on a few examples
        assert rpars.V0_IMAG == pytest.approx(5.0)
        assert rpars.THEO_ENERGIES == pytest.approx([50, 350, 3])


class TestIr1002x1OParameters():
    # checks reading of parameters:
    # RUN, THEO_ENERGIES, LMAX, BULK_LIKE_BELOW, SITE_DEF, T_DEBYE, T_EXPERIMENT, VIBR_AMP_SCALE
    def test_ir100_2x1_o_interpretation(self, ir100_2x1_o_parameters_example):
        rpars, slab = ir100_2x1_o_parameters_example
        # check that the parameters are interpreted correctly based on a few examples
        assert rpars.RUN == [0, 1, 2, 3]
        assert rpars.THEO_ENERGIES == pytest.approx([49, 700, 3])
        assert rpars.LMAX == [8, 14]
        assert rpars.BULK_LIKE_BELOW == pytest.approx(0.35)
        # check that float assignment works
        assert rpars.T_DEBYE == pytest.approx(420)
        assert rpars.T_EXPERIMENT == pytest.approx(100)
        # check that SITE_DEF assignment works
        assert rpars.SITE_DEF['Ir']['surf'] == [3,2]
        assert rpars.SITE_DEF['O']['ads'] == [1]
        # check that string assignment works for VIBR_AMP_SCALE
        assert rpars.VIBR_AMP_SCALE[0] == '*surf 1.3'
        # check that Pendry R factor is selected
        assert rpars.R_FACTOR_TYPE == 1


# unit tests for parameter interpretation


class TestIntpolDeg:
    def test_interpret_intpol_deg_valid(self, mock_rparams):
        param = 'INTPOL_DEG'
        for val in mock_rparams.get_limits(param):
            interpreter = ParameterInterpreter(mock_rparams)
            interpreter.interpret_intpol_deg(Assignment(val, param))
            assert mock_rparams.INTPOL_DEG == int(val)

    def test_interpret_intpol_def_invalid(self, mock_rparams):
        incompatible_values = ['1', 'text']
        for val in incompatible_values:
            interpreter = ParameterInterpreter(mock_rparams)
            assignment = Assignment(val, 'INTPOL_DEG')
            with pytest.raises(ParameterError):
                        interpreter.interpret_intpol_deg(assignment)


class TestNumericalParameter:
    def test_interpret_numerical_parameter_float(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.01", 'TEST_PARAM')
        bounds = NumericBounds()
        result = interpreter.interpret_numerical_parameter(assignment,
                                                           bounds=bounds)
        assert result == pytest.approx(0.01)
        assert mock_rparams.TEST_PARAM == pytest.approx(0.01)

    def test_interpret_numerical_parameter_int(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("100", 'TEST_PARAM')
        bounds = NumericBounds(type_=int)
        result = interpreter.interpret_numerical_parameter(assignment,
                                                           bounds=bounds)
        assert result == 100
        assert mock_rparams.TEST_PARAM == 100

    def test_interpret_numerical_parameter_int_at_bound_ok(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("100", 'TEST_PARAM')
        bounds = NumericBounds(type_=int, range_=(0, 100))
        result = interpreter.interpret_numerical_parameter(assignment,
                                                           bounds=bounds)
        assert result == 100
        assert mock_rparams.TEST_PARAM == 100

    def test_interpret_numerical_parameter_int_at_bound_fail(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("100", 'TEST_PARAM')
        bounds = NumericBounds(type_=int, range_=(0, 100),
                               accept_limits=(False, False))
        with pytest.raises(ParameterError):
            interpreter.interpret_numerical_parameter(assignment,
                                                      bounds=bounds)

    def test_interpret_numerical_parameter_float_out_of_range(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("-0.5", 'TEST_PARAM')
        bounds = NumericBounds(type_=float, range_=(0, 1))
        with pytest.raises(ParameterError):
            interpreter.interpret_numerical_parameter(assignment,
                                                      bounds=bounds)

    def test_interpret_numerical_parameter_float_at_bound_ok(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("-0.5", 'TEST_PARAM')
        bounds = NumericBounds(type_=float, range_=(-0.5, 1))
        result = interpreter.interpret_numerical_parameter(assignment,
                                                           bounds=bounds)
        assert result == -0.5
        assert mock_rparams.TEST_PARAM == -0.5

    def test_interpret_numerical_parameter_float_at_bound_fail(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("-0.5", 'TEST_PARAM')
        bounds = NumericBounds(type_=float, range_=(-0.5, 1),
                               accept_limits=(False, False))
        with pytest.raises(ParameterError):
            interpreter.interpret_numerical_parameter(assignment,
                                                      bounds=bounds)

    def test_interpret_numerical_parameter_float_out_of_range_event_modulo(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("2.5", 'TEST_PARAM')
        bounds = NumericBounds(type_=float, range_=(0, 2),
                               out_of_range_event='modulo')
        result = interpreter.interpret_numerical_parameter(assignment,
                                                           bounds=bounds)
        assert result == pytest.approx(0.5)
        assert mock_rparams.TEST_PARAM == pytest.approx(0.5)

    def test_interpret_numerical_parameter_int_out_of_range(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("-5", 'TEST_PARAM')
        bounds = NumericBounds(type_=int, range_=(0, 10))
        with pytest.raises(ParameterError):
            interpreter.interpret_numerical_parameter(assignment,
                                                      bounds=bounds)

    def test_interpret_numerical_parameter_int_out_of_range_event_modulo(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("12", 'TEST_PARAM')
        bounds = NumericBounds(type_=int,
                               range_=(0, 10),
                               out_of_range_event='modulo')
        result = interpreter.interpret_numerical_parameter(assignment,
                                                           bounds=bounds)
        assert result == 2
        assert mock_rparams.TEST_PARAM == 2


class TestSymmetryBulk:
    def test_interpret_symmetry_bulk_mirror(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("p2 m[0 1]", "SYMMETRY_BULK")
        interpreter.interpret_symmetry_bulk(assignment)
        assert mock_rparams.SYMMETRY_BULK == {
            'mirror': {(0, 1)},
            'rotation': set(),
            'group': 'p2'
        }

    def test_interpret_symmetry_cm_rot4(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("cm[1 1] r4", "SYMMETRY_BULK")
        interpreter.interpret_symmetry_bulk(assignment)
        assert mock_rparams.SYMMETRY_BULK == {
            'mirror': set(),
            'rotation': {4},
            'group': 'cm[1 1]'
        }

    def test_interpret_symmetry_bulk_pmm(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("pmm", "SYMMETRY_BULK")
        interpreter.interpret_symmetry_bulk(assignment)
        assert mock_rparams.SYMMETRY_BULK == {
            'mirror': set(),
            'rotation': set(),
            'group': 'pmm',
        }

    def test_interpret_symmetry_bulk_invalid_syntax(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("invalid_syntax", "SYMMETRY_BULK")
        with pytest.raises(ParameterValueError):
            interpreter.interpret_symmetry_bulk(assignment)

    def test_interpret_symmetry_bulk_multiple_groups(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("pmm pg", "SYMMETRY_BULK")
        with pytest.raises(ParameterValueError):
            interpreter.interpret_symmetry_bulk(assignment)

    def test_interpret_symmetry_bulk_missing_group(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("", "SYMMETRY_BULK")
        with pytest.raises(ParameterValueError):
            interpreter.interpret_symmetry_bulk(assignment)


class TestBulkRepeat():
    param = 'BULK_REPEAT'
    rpars = Rparams()

    invalid_inputs = ['text', 'y(1.2)', 'z(abc)', '[]']
    valid_inputs = ['c(1.2)', 'z(0.9)']

    def test_interpret_bulk_repeat_invalid(self, slab_ag100):
        slab = slab_ag100
        for val in self.invalid_inputs:
            with pytest.raises(ParameterError):
                interpreter = ParameterInterpreter(self.rpars)
                interpreter.slab = slab
                interpreter.interpret_bulk_repeat(Assignment(val, 'BULK_REPEAT'))


    def test_interpret_bulk_repeat_float(self, slab_ag100):
        val = '1.428'
        interpreter = ParameterInterpreter(self.rpars)
        interpreter.slab = slab_ag100
        interpreter.interpret_bulk_repeat(Assignment(val, 'BULK_REPEAT'))
        assert self.rpars.BULK_REPEAT == pytest.approx(1.428, rel=1e-4)

    def test_interpret_bulk_repeat_c(self, slab_ag100):
        val = 'c(0.1)'
        interpreter = ParameterInterpreter(self.rpars)
        interpreter.slab = slab_ag100
        interpreter.interpret_bulk_repeat(Assignment(val, 'BULK_REPEAT'))
        assert self.rpars.BULK_REPEAT == pytest.approx(2.03646, rel=1e-4)

    def test_interpret_bulk_repeat_z(self, slab_ag100):
        val = 'c(0.1)'
        interpreter = ParameterInterpreter(self.rpars)
        interpreter.slab = slab_ag100
        interpreter.interpret_bulk_repeat(Assignment(val, 'BULK_REPEAT'))
        assert self.rpars.BULK_REPEAT == pytest.approx(2.0364, rel=1e-4)

    def test_interpret_bulk_repeat_vector(self, slab_ag100):
        val = '[1.0 2.0 3.0]'
        interpreter = ParameterInterpreter(self.rpars)
        interpreter.slab = slab_ag100
        interpreter.interpret_bulk_repeat(Assignment(val, 'BULK_REPEAT'))
        assert self.rpars.BULK_REPEAT == pytest.approx([1.0, 2.0, 3.0], rel=1e-4)


class TestFortranComp():
    # TODO: make use of new intepreter class
    param = 'FORTRAN_COMP'
    rpars = Rparams()
    interpreter = ParameterInterpreter(rpars)

    def test_fortran_comp_default_intel(self):
        assignment = Assignment('ifort', 'FORTRAN_COMP')
        self.interpreter.interpret_fortran_comp(assignment, skip_check=True)
        assert 'ifort -O2 -I/opt/intel/mkl/include' in self.rpars.FORTRAN_COMP[0]
        assert '-L/opt/intel/mkl/lib/intel64' in self.rpars.FORTRAN_COMP[1]

    def test_fortran_comp_default_gnu(self):
        assignment = Assignment('gfortran', 'FORTRAN_COMP')
        self.interpreter.interpret_fortran_comp(assignment, skip_check=True)
        assert 'gfortran' in self.rpars.FORTRAN_COMP[0]
        assert '-llapack -lpthread -lblas' in self.rpars.FORTRAN_COMP[1]

    def test_fortran_comp_default_intel_mpi(self):
        assignment = Assignment('mpiifort', 'FORTRAN_COMP', flags_str='mpi')
        self.interpreter.interpret_fortran_comp(assignment, skip_check=True)
        assert 'mpiifort' in self.rpars.FORTRAN_COMP_MPI[0]  # no other flags by default

    def test_fortran_comp_default_gnu_mpi(self):
        assignment = Assignment('mpifort', 'FORTRAN_COMP', flags_str='mpi')
        self.interpreter.interpret_fortran_comp(assignment, skip_check=True)
        assert 'mpifort -Ofast -no-pie' in self.rpars.FORTRAN_COMP_MPI[0]  # no other flags by default

    def test_fortran_comp_custom_without_flag(self):
        assignment = Assignment("'ifort -O3 -march=native'", 'FORTRAN_COMP')
        self.interpreter.interpret_fortran_comp(assignment, skip_check=True)
        assert 'ifort -O3 -march=native' in self.rpars.FORTRAN_COMP[0]

    def test_fortran_comp_custom_post_flag(self):
        assignment = Assignment("'-L/opt/intel/mkl/lib/intel64'",
                                'FORTRAN_COMP', flags_str='post')
        self.interpreter.interpret_fortran_comp(assignment, skip_check=True)
        assert '-L/opt/intel/mkl/lib/intel64' in self.rpars.FORTRAN_COMP[1]

    def test_fortran_comp_custom_mpi_flag(self):
        assignment = Assignment("'mpifort -fallow-argument-mismatch'",
                                'FORTRAN_COMP', flags_str='mpi')
        self.interpreter.interpret_fortran_comp(assignment, skip_check=True)
        assert 'mpifort -fallow-argument-mismatch' in self.rpars.FORTRAN_COMP_MPI[0]


class TestV0Real():
    rpars = Rparams()
    interpreter = ParameterInterpreter(rpars)

    def test_interpret_v0_real_rundgren_type(self):
        assignment = Assignment("rundgren 1.0 2.0 3.0 4.0", "V0_REAL")
        self.interpreter.interpret_v0_real(assignment)
        assert self.rpars.V0_REAL == [1.0, 2.0, 3.0, 4.0]

    def test_interpret_v0_real_other_type(self,):
        assignment = Assignment("EE", "V0_REAL")
        self.interpreter.interpret_v0_real(assignment)
        assert self.rpars.V0_REAL == "EEV+workfn"

    def test_interpret_v0_real_invalid_rundgren_constants(self):
        assignment = Assignment("rundgren 1.0 2.0", "V0_REAL")
        with pytest.raises(ParameterError):
            self.interpreter.interpret_v0_real(assignment)

    def test_interpret_v0_real_value_error(self):
        assignment = Assignment("rundgren 1.0 2.0 3.0 four", "V0_REAL")
        with pytest.raises(ParameterError):
            self.interpreter.interpret_v0_real(assignment)


class TestTheoEnergies:
    def test_interpret_theo_energies_single_value_default(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("_", "THEO_ENERGIES")
        interpreter.interpret_theo_energies(assignment)
        assert mock_rparams.THEO_ENERGIES == pytest.approx(mock_rparams.get_default("THEO_ENERGIES"))

    def test_interpret_theo_energies_single_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("1.0", "THEO_ENERGIES")
        interpreter.interpret_theo_energies(assignment)
        assert mock_rparams.THEO_ENERGIES == pytest.approx([1.0, 1.0, 1.0])

    def test_interpret_theo_energies_three_values_all_default(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("_ _ _", "THEO_ENERGIES")
        interpreter.interpret_theo_energies(assignment)
        assert mock_rparams.THEO_ENERGIES == pytest.approx(mock_rparams.get_default("THEO_ENERGIES"))

    def test_interpret_theo_energies_three_values_partial_default(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("1.0 _ 2.0", "THEO_ENERGIES")
        interpreter.interpret_theo_energies(assignment)
        assert mock_rparams.THEO_ENERGIES == pytest.approx([1.0, mock_rparams.get_default("THEO_ENERGIES")[1], 2.0])

    def test_interpret_theo_energies_three_values_positive_floats(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("1.0 2.0 0.5", "THEO_ENERGIES")
        interpreter.interpret_theo_energies(assignment)
        assert mock_rparams.THEO_ENERGIES == pytest.approx([1.0, 2.0, 0.5])

    def test_interpret_theo_energies_two_defaults(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("1.0 _ _", "THEO_ENERGIES")
        interpreter.interpret_theo_energies(assignment)
        assert mock_rparams.THEO_ENERGIES == pytest.approx([1.0, mock_rparams.get_default("THEO_ENERGIES")[1], mock_rparams.get_default("THEO_ENERGIES")[2]])

    def test_interpret_theo_energies_three_values_negative_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("1.0 -2.0 0.5", "THEO_ENERGIES")
        with pytest.raises(ParameterRangeError):
            interpreter.interpret_theo_energies(assignment)

    def test_interpret_theo_energies_three_values_invalid_range(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("2.0 1.0 0.5", "THEO_ENERGIES")
        with pytest.raises(ParameterValueError):
            interpreter.interpret_theo_energies(assignment)

    def test_interpret_theo_energies_three_values_range_correction(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("1.0 2.5 0.5", "THEO_ENERGIES")
        interpreter.interpret_theo_energies(assignment)
        assert mock_rparams.THEO_ENERGIES == pytest.approx([1.0, 2.5, 0.5])

    def test_interpret_theo_energies_three_values_out_of_range_start(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("-0.5 2.0 0.5", "THEO_ENERGIES")
        with pytest.raises(ParameterRangeError):
            interpreter.interpret_theo_energies(assignment)


class TestNumericalParamsExamples:
    def test_interpret_n_bulk_layers_valid(self):
        rpars = Rparams()
        interpreter = ParameterInterpreter(rpars)
        assignment = Assignment("1", "N_BULK_LAYERS")
        interpreter.interpret_n_bulk_layers(assignment)
        assert rpars.N_BULK_LAYERS == 1

    def test_interpret_n_bulk_layers_invalid(self):
        # N_BULK_LAYERS must be 1 or 2
        rpars = Rparams()
        interpreter = ParameterInterpreter(rpars)
        assignment = Assignment("3", "N_BULK_LAYERS")
        with pytest.raises(ParameterError):
            interpreter.interpret_n_bulk_layers(assignment)

    def test_interpret_t_debye(self):
        rpars = Rparams()
        interpreter = ParameterInterpreter(rpars)
        assignment = Assignment("300.0", "T_DEBYE")
        interpreter.interpret_t_debye(assignment)
        assert rpars.T_DEBYE == pytest.approx(300.0)


class TestLogLevel:
    def test_interpret_log_debug_true(self):
        rpars = Rparams()
        interpreter = ParameterInterpreter(rpars)
        assignment = Assignment("true", "LOG_LEVEL")
        interpreter.interpret_log_level(assignment)
        assert rpars.LOG_LEVEL <= 10

    def test_interpret_log_debug_false(self):
        rpars = Rparams()
        interpreter = ParameterInterpreter(rpars)
        assignment = Assignment("F", "LOG_LEVEL")
        interpreter.interpret_log_level(assignment)
        assert rpars.LOG_LEVEL >= 20

    def test_interpret_log_debug_int(self):
        rpars = Rparams()
        interpreter = ParameterInterpreter(rpars)
        assignment = Assignment("3", "LOG_LEVEL")
        interpreter.interpret_log_level(assignment)
        assert rpars.LOG_LEVEL == 3


class TestFilamentWF:
    def test_interpret_filament_wf_lab6(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("LaB6", "FILAMENT_WF")
        interpreter.interpret_filament_wf(assignment)
        assert mock_rparams.FILAMENT_WF == pytest.approx(2.65)

    def test_interpret_filament_wf_custom(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("1.0", "FILAMENT_WF")
        interpreter.interpret_filament_wf(assignment)
        assert mock_rparams.FILAMENT_WF == pytest.approx(1.0)

    def test_interpret_filament_wf_invalid(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("invalid", "FILAMENT_WF")
        with pytest.raises(ParameterFloatConversionError):
            interpreter.interpret_filament_wf(assignment)

    def test_interpret_filament_wf_flag_invalid(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("1.5", "FILAMENT_WF", flags_str='test')
        with pytest.raises(ParameterUnknownFlagError):
            interpreter.interpret_filament_wf(assignment)


class TestAverageBeams:
    def test_interpret_average_beams_off(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("off", "AVERAGE_BEAMS")
        interpreter.interpret_average_beams(assignment)
        assert mock_rparams.AVERAGE_BEAMS is False

    def test_interpret_average_beams_all(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("all", "AVERAGE_BEAMS")
        interpreter.interpret_average_beams(assignment)
        assert mock_rparams.AVERAGE_BEAMS == (0.0, 0.0)

    def test_interpret_average_beams_custom_angles(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("45.0 60.0", "AVERAGE_BEAMS")
        interpreter.interpret_average_beams(assignment)
        assert mock_rparams.AVERAGE_BEAMS == (45.0, 60.0)

    def test_interpret_average_beams_invalid_input(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("invalid input", "AVERAGE_BEAMS")
        with pytest.raises(ParameterError):
            interpreter.interpret_average_beams(assignment)


class TestBeamIncidence:
    def test_interpret_beam_incidence_custom_angles(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("45.0 60.0", "BEAM_INCIDENCE")
        interpreter.interpret_beam_incidence(assignment)
        assert mock_rparams.THETA == 45.0
        assert mock_rparams.PHI == 60.0

    def test_interpret_beam_incidence_invalid_input(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("invalid input", "BEAM_INCIDENCE")
        with pytest.raises(ParameterError):
            interpreter.interpret_beam_incidence(assignment)


class TestDomain:
    def test_interpret_domain_existing_path(self, mock_rparams, tmp_path):
        interpreter = ParameterInterpreter(mock_rparams)
        domain_path = tmp_path / "domain1"
        domain_path.mkdir()
        assignment = Assignment(values_str=str(domain_path),
                                parameter="DOMAINS",
                                flags_str="domain1")
        interpreter.interpret_domain(assignment)
        assert mock_rparams.DOMAINS == [("domain1", str(domain_path))]

    def test_interpret_domain_new_path(self, mock_rparams, tmp_path):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(str(tmp_path), "DOMAINS")
        interpreter.interpret_domain(assignment)
        assert mock_rparams.DOMAINS == [("1", str(tmp_path))]

    def test_interpret_domain_zip_file(self, mock_rparams, tmp_path):
        zip_file = tmp_path / "domain.zip"
        zip_file.touch()
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(str(zip_file), "DOMAINS")
        interpreter.interpret_domain(assignment)
        assert mock_rparams.DOMAINS == [("1", str(zip_file))]

    def test_interpret_domain_invalid_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("invalid_path", "DOMAINS")
        with pytest.raises(ParameterError):
            interpreter.interpret_domain(assignment)


class TestDomainStep:
    def test_interpret_domain_step_valid_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("10", "DOMAIN_STEP")
        interpreter.interpret_domain_step(assignment)
        assert mock_rparams.DOMAIN_STEP == 10

    def test_interpret_domain_step_invalid_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("200", "DOMAIN_STEP")
        with pytest.raises(ParameterRangeError):
            interpreter.interpret_domain_step(assignment)

    def test_interpret_domain_step_non_integer_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.5", "DOMAIN_STEP")
        with pytest.raises(ParameterError):
            interpreter.interpret_domain_step(assignment)


class TestSuperlattice:
    def test_interpret_superlattice_matrix_notation(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(values_str="2 0, 0 2",
                                parameter="SUPERLATTICE",
                                flags_str="M")
        interpreter.interpret_superlattice(assignment)
        assert mock_rparams.SUPERLATTICE == pytest.approx(np.array([[2, 0], [0, 2,]]))
        assert mock_rparams.superlattice_defined is True

    def test_interpret_superlattice_woods_notation(self, mock_rparams, slab_ag100):
        interpreter = ParameterInterpreter(mock_rparams)
        interpreter.slab = slab_ag100
        assignment = Assignment("p(2x1)", "SUPERLATTICE")
        interpreter.interpret_superlattice(assignment)
        assert mock_rparams.SUPERLATTICE == pytest.approx(np.array([[2, 0], [0, 1,]]))
        assert mock_rparams.superlattice_defined is True

    def test_interpret_superlattice_no_slab(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("c(2x2)", "SUPERLATTICE")
        # should raise because slab is None
        with pytest.raises(ParameterError):
            interpreter.interpret_superlattice(assignment)


class TestTensorOutput:
    def test_interpret_tensor_output_single_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("False", "TENSOR_OUTPUT")
        interpreter.interpret_tensor_output(assignment)
        assert mock_rparams.TENSOR_OUTPUT == [0]

    def test_interpret_tensor_output_multiple_values(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0 1 1 0", "TENSOR_OUTPUT")
        interpreter.interpret_tensor_output(assignment)
        assert mock_rparams.TENSOR_OUTPUT == [0, 1, 1, 0]

    def test_interpret_tensor_output_repeated_values(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("2*1 2*0 1", "TENSOR_OUTPUT")
        interpreter.interpret_tensor_output(assignment)
        assert mock_rparams.TENSOR_OUTPUT == [1, 1, 0, 0, 1]

    def test_interpret_tensor_output_invalid_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("2", "TENSOR_OUTPUT")
        with pytest.raises(ParameterParseError):
            interpreter.interpret_tensor_output(assignment)


class TestElementRename:
    def test_interpret_element_rename_valid_assignment(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(parameter="ELEMENT_RENAME",
                                flags_str="X",
                                values_str="H")
        interpreter.interpret_element_rename(assignment)
        assert mock_rparams.ELEMENT_RENAME == {"X": "H"}

    def test_interpret_element_rename_invalid_element(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(parameter="ELEMENT_RENAME",
                                flags_str="Uk",
                                values_str="Op")
        with pytest.raises(ParameterError):
            interpreter.interpret_element_rename(assignment)


class TestSymmetryFix:
    def test_interpret_symmetry_fix_auto(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(parameter="SYMMETRY_FIX",
                                values_str="t")
        interpreter.interpret_symmetry_fix(assignment)
        assert mock_rparams.SYMMETRY_FIX == ''

    def test_interpret_symmetry_fix_p1(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(parameter="SYMMETRY_FIX",
                                values_str="f")
        interpreter.interpret_symmetry_fix(assignment)
        assert mock_rparams.SYMMETRY_FIX == "p1"

    def test_interpret_symmetry_fix_invalid_group(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(parameter="SYMMETRY_FIX",
                                values_str="invalid")
        with pytest.raises(ParameterParseError):
            interpreter.interpret_symmetry_fix(assignment)

    def test_interpret_symmetry_fix_direction_required(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(parameter="SYMMETRY_FIX",
                                values_str="cm")
        with pytest.raises(ParameterParseError):
            interpreter.interpret_symmetry_fix(assignment)

    def test_interpret_symmetry_fix_custom_group_with_direction(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(parameter="SYMMETRY_FIX",
                                values_str="cm[1 1]")
        interpreter.interpret_symmetry_fix(assignment)
        assert mock_rparams.SYMMETRY_FIX == "cm[1 1]"

    def test_interpret_symmetry_fix_invalid_direction(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(parameter="SYMMETRY_FIX",
                                values_str="pmt [0 x]")
        with pytest.raises(ParameterError):
            interpreter.interpret_symmetry_fix(assignment)


class TestIVShiftRange:
    def test_interpret_iv_shift_range_valid_range(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.0 1.0 0.25", "IV_SHIFT_RANGE")
        interpreter.interpret_iv_shift_range(assignment)
        assert mock_rparams.IV_SHIFT_RANGE == pytest.approx([0.0, 1.0, 0.25])

    def test_interpret_iv_shift_range_invalid_number_of_inputs(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.0", "IV_SHIFT_RANGE")
        with pytest.raises(ParameterNumberOfInputsError):
            interpreter.interpret_iv_shift_range(assignment)

    def test_interpret_iv_shift_range_invalid_float_conversion(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.0 2.0 ()", "IV_SHIFT_RANGE")
        with pytest.raises(ParameterFloatConversionError):
            interpreter.interpret_iv_shift_range(assignment)

    def test_interpret_iv_shift_range_invalid_parse(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.0 2.0 1.a", "IV_SHIFT_RANGE")
        with pytest.raises(ParameterParseError):
            interpreter.interpret_iv_shift_range(assignment)

    def test_interpret_iv_shift_range_end_energy_less_than_start_energy(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("1.0 0.0 0.1", "IV_SHIFT_RANGE")
        with pytest.raises(ParameterError):
            interpreter.interpret_iv_shift_range(assignment)

    def test_interpret_iv_shift_range_invalid_step(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.0 1.0 -0.1", "IV_SHIFT_RANGE")
        with pytest.raises(ParameterError):
            interpreter.interpret_iv_shift_range(assignment)


class TestOptimize:
    def test_interpret_optimize_valid_flag_and_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(parameter= "OPTIMIZE",
                                flags_str="v0i", values_str="step 0.1")
        interpreter.interpret_optimize(assignment)
        assert mock_rparams.OPTIMIZE['which'] == 'v0i'
        assert mock_rparams.OPTIMIZE['step'] == pytest.approx(0.1)

    def test_interpret_optimize_valid_flag_and_value_multiple(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(
            parameter= "OPTIMIZE",
            flags_str="theta",
            values_str="step 0.1, convergence 1e-6, minpoints 10"
            )
        interpreter.interpret_optimize(assignment)
        assert mock_rparams.OPTIMIZE['which'] == 'theta'
        assert mock_rparams.OPTIMIZE['step'] == pytest.approx(0.1)
        assert mock_rparams.OPTIMIZE['minpoints'] == pytest.approx(10)
        assert mock_rparams.OPTIMIZE['convergence'] == pytest.approx(1e-6)

    def test_interpret_optimize_invalid_flag(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("invalid 0.1", "OPTIMIZE")
        with pytest.raises(ParameterUnknownFlagError):
            interpreter.interpret_optimize(assignment)

    def test_interpret_optimize_invalid_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("step not-a-number", "OPTIMIZE")
        with pytest.raises(ParameterError):
            interpreter.interpret_optimize(assignment)

    def test_interpret_optimize_invalid_flag(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("step 0.1", "OPTIMIZE", flags_str="invalid")
        with pytest.raises(ParameterError):
            interpreter.interpret_optimize(assignment)


class TestPhaseshiftEps:
    def test_interpret_phaseshift_eps_valid_float(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.1", "PHASESHIFT_EPS")
        interpreter.interpret_phaseshift_eps(assignment)
        assert mock_rparams.PHASESHIFT_EPS == 0.1

    def test_interpret_phaseshift_eps_invalid(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("invalid", "PHASESHIFT_EPS")
        with pytest.raises(ParameterError):
            interpreter.interpret_phaseshift_eps(assignment)

    def test_interpret_phaseshift_eps_default_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("fine", "PHASESHIFT_EPS")
        interpreter.interpret_phaseshift_eps(assignment)
        assert mock_rparams.PHASESHIFT_EPS == 0.01

    def test_interpret_phaseshift_eps_out_of_range(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("-1.0", "PHASESHIFT_EPS")
        with pytest.raises(ParameterRangeError):
            interpreter.interpret_phaseshift_eps(assignment)

        assignment = Assignment("1.5", "PHASESHIFT_EPS")
        with pytest.raises(ParameterRangeError):
            interpreter.interpret_phaseshift_eps(assignment)


class TestLayerCuts:
    def test_interpret_layer_cuts_simple_values(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.1 0.2 0.3", "LAYER_CUTS")
        interpreter.interpret_layer_cuts(assignment)
        assert mock_rparams.LAYER_CUTS == ["0.1", "0.2", "0.3"]

    def test_interpret_layer_cuts_less_than_greater_than(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("< 0.1 > 0.2", "LAYER_CUTS")
        with pytest.raises(ParameterParseError):
            interpreter.interpret_layer_cuts(assignment)

    def test_interpret_layer_cuts_dz_cutoff(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("dz(1.2)", "LAYER_CUTS")
        interpreter.interpret_layer_cuts(assignment)
        assert mock_rparams.LAYER_CUTS == ["dz(1.2)"]

    def test_interpret_layer_cuts_dc_cutoff(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.15 0.3 < dz(1.2)  ", "LAYER_CUTS")
        interpreter.interpret_layer_cuts(assignment)
        assert mock_rparams.LAYER_CUTS == ["0.15", "0.3", "<", "dz(1.2)"]

    def test_interpret_layer_cuts_invalid_function_format(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("dz(0.1) dc(0.2) invalid(0.3)", "LAYER_CUTS")
        with pytest.raises(ParameterParseError):
            interpreter.interpret_layer_cuts(assignment)

    def test_interpret_layer_cuts_invalid_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.1 invalid 0.3", "LAYER_CUTS")
        with pytest.raises(ParameterParseError):
            interpreter.interpret_layer_cuts(assignment)

    def test_interpret_layer_cuts_two_dz(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("dz(1.0) < 2.0 < dz(0.5) < 4.0", "LAYER_CUTS")
        interpreter.interpret_layer_cuts(assignment)
        assert mock_rparams.LAYER_CUTS == ["dz(1.0)", "<", "2.0", "<", 
                                           "dz(0.5)", "<", "4.0"]

    def test_interpret_layer_cuts_dz_invalid(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.5 1.0 < dz(abcd) < 4.0", "LAYER_CUTS")
        with pytest.raises(ParameterError):
            interpreter.interpret_layer_cuts(assignment)



class TestRun:
    def test_interpret_run_single_section(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("1", "RUN")
        interpreter.interpret_run(assignment)
        assert mock_rparams.RUN == [0, 1]

    def test_interpret_run_multiple_sections(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("1 2 3", "RUN")
        interpreter.interpret_run(assignment)
        assert mock_rparams.RUN == [0, 1, 2, 3]

    def test_interpret_run_range(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("1-3", "RUN")
        interpreter.interpret_run(assignment)
        assert mock_rparams.RUN == [0, 1, 2, 3]

    def test_interpret_run_invalid_section(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("invalid", "RUN")
        with pytest.raises(ParameterValueError):
            interpreter.interpret_run(assignment)

    def test_interpret_run_empty(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("", "RUN")
        with pytest.raises(ParameterError):
            interpreter.interpret_run(assignment)

    def test_interpret_run_invalid_syntax_underscore(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("1 _ 2", "RUN")
        with pytest.raises(ParameterValueError):
            interpreter.interpret_run(assignment)

    def test_interpret_run_invalid_syntax_section(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("1, x", "RUN")
        with pytest.raises(ParameterValueError):
            interpreter.interpret_run(assignment)

class TestSymmetryEps:
    def test_interpret_symmetry_eps_single_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.1", "SYMMETRY_EPS")
        interpreter.interpret_symmetry_eps(assignment)
        assert mock_rparams.SYMMETRY_EPS == 0.1

    def test_interpret_symmetry_eps_multiple_values(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.1 0.2", "SYMMETRY_EPS")
        interpreter.interpret_symmetry_eps(assignment)
        assert mock_rparams.SYMMETRY_EPS == 0.1
        assert mock_rparams.SYMMETRY_EPS_Z == 0.2

    def test_interpret_symmetry_eps_invalid_number_of_inputs(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.1 0.2 0.3", "SYMMETRY_EPS")
        with pytest.raises(ParameterNumberOfInputsError):
            interpreter.interpret_symmetry_eps(assignment)


class TestSearchBeams():
    def test_interpret_search_beams_zero(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment('0', parameter='SEARCH_BEAMS')
        interpreter.interpret_search_beams(assignment)
        assert mock_rparams.SEARCH_BEAMS == 0

    def test_interpret_search_beams_alternative_zero(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment('A', parameter='SEARCH_BEAMS')
        interpreter.interpret_search_beams(assignment)
        assert mock_rparams.SEARCH_BEAMS == 0

    def test_interpret_search_beams_one(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment('1', parameter='SEARCH_BEAMS')
        interpreter.interpret_search_beams(assignment)
        assert mock_rparams.SEARCH_BEAMS == 1

    def test_interpret_search_beams_alternative_one(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment('I', parameter='SEARCH_BEAMS')
        interpreter.interpret_search_beams(assignment)
        assert mock_rparams.SEARCH_BEAMS == 1

    def test_interpret_search_beams_two(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment('2', parameter='SEARCH_BEAMS')
        interpreter.interpret_search_beams(assignment)
        assert mock_rparams.SEARCH_BEAMS == 2

    def test_interpret_search_beams_alternative_two(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment('F', parameter='SEARCH_BEAMS')
        interpreter.interpret_search_beams(assignment)
        assert mock_rparams.SEARCH_BEAMS == 2

    def test_interpret_search_beams_invalid_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment('3', parameter='SEARCH_BEAMS')
        with pytest.raises(ParameterValueError):
            interpreter.interpret_search_beams(assignment)

    def test_interpret_search_beams_multiple_values(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment('0 1', parameter='SEARCH_BEAMS')
        with pytest.raises(ParameterError):
            interpreter.interpret_search_beams(assignment)


class TestSearchCull:
    def test_interpret_search_cull_single_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.5", "SEARCH_CULL")
        interpreter.interpret_search_cull(assignment)
        assert mock_rparams.SEARCH_CULL == 0.5

    def test_interpret_search_cull_integer_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("3", "SEARCH_CULL")
        interpreter.interpret_search_cull(assignment)
        assert mock_rparams.SEARCH_CULL == 3

    def test_interpret_search_cull_float_bigger_one(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("1.5", "SEARCH_CULL")
        with pytest.raises(ParameterError):
            interpreter.interpret_search_cull(assignment)


    def test_interpret_search_cull_greater_than_one_integer(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("4.0", "SEARCH_CULL")
        interpreter.interpret_search_cull(assignment)
        assert mock_rparams.SEARCH_CULL == pytest.approx(4)

    def test_interpret_search_cull_negative_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("-0.5", "SEARCH_CULL")
        with pytest.raises(ParameterError):
            interpreter.interpret_search_cull(assignment)

    def test_interpret_search_cull_invalid_number_of_inputs(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.5 clone 1.0", "SEARCH_CULL")
        with pytest.raises(ParameterNumberOfInputsError):
            interpreter.interpret_search_cull(assignment)

    def test_interpret_search_cull_invalid_search_cull_type(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.5 test", "SEARCH_CULL")
        with pytest.raises(ParameterValueError):
            interpreter.interpret_search_cull(assignment)

class TestSearchConvergence:
    def test_interpret_search_convergence_gaussian(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(parameter="SEARCH_CONVERGENCE",
                                flags_str="gaussian",
                                values_str="0.01 0.9")
        interpreter.interpret_search_convergence(assignment)
        assert mock_rparams.GAUSSIAN_WIDTH == 0.01
        assert mock_rparams.GAUSSIAN_WIDTH_SCALING == 0.9

    def test_interpret_search_convergence_gaussian_no_scaling(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(parameter="SEARCH_CONVERGENCE",
                                flags_str="gaussian",
                                values_str="0.01")
        interpreter.interpret_search_convergence(assignment)
        assert mock_rparams.GAUSSIAN_WIDTH == 0.01
        assert mock_rparams.GAUSSIAN_WIDTH_SCALING == 0.5

    def test_interpret_search_convergence_invalid_flag(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(parameter="SEARCH_CONVERGENCE",
                                flags_str="test",
                                values_str=" 0.01 0.9")
        with pytest.raises(ParameterUnknownFlagError):
            interpreter.interpret_search_convergence(assignment)

    def test_interpret_search_convergence_invalid_number_of_inputs(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(parameter="SEARCH_CONVERGENCE",
                                flags_str="gaussian",
                                values_str="0.01 0.9 0.5")
        with pytest.raises(ParameterNumberOfInputsError):
            interpreter.interpret_search_convergence(assignment)

    def test_interpret_search_convergence_invalid_scaling_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(parameter="SEARCH_CONVERGENCE",
                                flags_str="gaussian",
                                values_str="0.01 -0.5")
        with pytest.raises(ParameterError):
            interpreter.interpret_search_convergence(assignment)

