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
                                                  interpretPARAMETERS,
                                                  ParameterInterpreter,
                                                  Assignment)
from viperleed.tleedmlib.files.poscar import readPOSCAR
from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tleedmlib.files.parameter_errors import (
    ParameterError, ParameterValueError, ParameterParseError,
    ParameterIntConversionError, ParameterFloatConversionError,
    ParameterBooleanConversionError, ParameterNotRecognizedError,
    ParameterNumberOfInputsError, ParameterRangeError,
    ParameterUnknownFlagError, ParameterNeedsFlagError
    )

_FIXTURES_PATH = Path('tests/fixtures/')

@pytest.fixture()
def ag100_parameters_example():
    # read Ag(100) POSCAR and PARAMETERS files
    slab = readPOSCAR(_FIXTURES_PATH / 'Ag(100)' / 'initialization' / 'POSCAR')
    rpars = readPARAMETERS(_FIXTURES_PATH / 'Ag(100)' / 'initialization' / 'PARAMETERS')
    # interpret PARAMETERS file
    interpretPARAMETERS(rpars, slab)
    return (rpars, slab)

@pytest.fixture()
def slab_ag100():
    # read Ag(100) POSCAR
    return readPOSCAR(_FIXTURES_PATH / 'POSCARs' / 'POSCAR_Ag(100)')

@pytest.fixture()
def slab_ir100_2x1_o():
    # read Ir(100)-(2x1)-O POSCAR
    return readPOSCAR(_FIXTURES_PATH / 'POSCARs' / 'POSCAR_Ir(100)-(2x1)-O')

@pytest.fixture()
def ir100_2x1_o_parameters_example(slab_ir100_2x1_o):
    slab = slab_ir100_2x1_o
    rpars = readPARAMETERS(_FIXTURES_PATH / 'parameters' / 'PARAMETERS_Ir(100)-(2x1)-O')
    interpretPARAMETERS(rpars, slab)
    return (rpars, slab)

@pytest.fixture()
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

# _interpret_intpol_deg()
class TestIntpolDeg:
    rpars = Rparams()
    def test__interpret_intpol_deg_valid(self):
        for val in self.rpars.get_limits('INTPOL_DEG'):
            interpreter = ParameterInterpreter(self.rpars, slab=None)
            interpreter._interpret_intpol_deg(Assignment(val))
            assert self.rpars.INTPOL_DEG == int(val)

    def test__interpret_intpol_def_invalid(self):
        incompatible_values = ['1', 'text']
        for val in incompatible_values:
            with pytest.raises(ParameterError):
                        interpreter = ParameterInterpreter(self.rpars, slab=None)
                        interpreter._interpret_intpol_deg(Assignment(val))

# _interpret_bulk_repeat()
class TestInterpretBulkRepeat():
    param = 'BULK_REPEAT'
    rpars = Rparams()

    invalid_inputs = ['text', 'y(1.2)', 'z(abc)', '[]']
    valid_inputs = ['c(1.2)', 'z(0.9)']

    def test__interpret_bulk_repeat_invalid(self, slab_ag100):
        slab = slab_ag100
        for val in self.invalid_inputs:
            with pytest.raises(ParameterError):
                interpreter = ParameterInterpreter(self.rpars, slab=None)
                interpreter._interpret_bulk_repeat(Assignment(val))


    def test__interpret_bulk_repeat_float(self, slab_ag100):
        val = '1.428'
        interpreter = ParameterInterpreter(self.rpars, slab=slab_ag100)
        interpreter._interpret_bulk_repeat(Assignment(val))
        assert self.rpars.BULK_REPEAT == pytest.approx(1.428, rel=1e-4)

    def test__interpret_bulk_repeat_c(self, slab_ag100):
        val = 'c(0.1)'
        interpreter = ParameterInterpreter(self.rpars, slab=slab_ag100)
        interpreter._interpret_bulk_repeat(Assignment(val))
        assert self.rpars.BULK_REPEAT == pytest.approx(2.03646, rel=1e-4)

    def test__interpret_bulk_repeat_z(self, slab_ag100):
        val = 'c(0.1)'
        interpreter = ParameterInterpreter(self.rpars, slab=slab_ag100)
        interpreter._interpret_bulk_repeat(Assignment(val))
        assert self.rpars.BULK_REPEAT == pytest.approx(2.0364, rel=1e-4)

    def test__interpret_bulk_repeat_vector(self, slab_ag100):
        val = '[1.0 2.0 3.0]'
        interpreter = ParameterInterpreter(self.rpars, slab=slab_ag100)
        interpreter._interpret_bulk_repeat(Assignment(val))
        assert self.rpars.BULK_REPEAT == pytest.approx([1.0, 2.0, 3.0], rel=1e-4)

# _interpret_fortran_comp()
class TestInterpretFortranComp():
    # TODO: make use of new intepreter class
    param = 'FORTRAN_COMP'
    rpars = Rparams()
    interpreter = ParameterInterpreter(rpars, slab=None)

    def test__fortran_comp_default_intel(self):
        assignment = Assignment('ifort')
        self.interpreter._interpret_fortran_comp(assignment, skip_check=True)
        assert 'ifort -O2 -I/opt/intel/mkl/include' in self.rpars.FORTRAN_COMP[0]
        assert '-L/opt/intel/mkl/lib/intel64' in self.rpars.FORTRAN_COMP[1]

    def test__fortran_comp_default_gnu(self):
        assignment = Assignment('gfortran')
        self.interpreter._interpret_fortran_comp(assignment, skip_check=True)
        assert 'gfortran' in self.rpars.FORTRAN_COMP[0]
        assert '-llapack -lpthread -lblas' in self.rpars.FORTRAN_COMP[1]

    def test__fortran_comp_default_intel_mpi(self):
        assignment = Assignment('mpiifort', flags ='mpi')
        self.interpreter._interpret_fortran_comp(assignment, skip_check=True)
        assert 'mpiifort' in self.rpars.FORTRAN_COMP_MPI[0]  # no other flags by default

    def test__fortran_comp_default_gnu_mpi(self):
        assignment = Assignment('mpifort', flags='mpi')
        self.interpreter._interpret_fortran_comp(assignment, skip_check=True)
        assert 'mpifort -Ofast -no-pie' in self.rpars.FORTRAN_COMP_MPI[0]  # no other flags by default

    def test__fortran_comp_custom_without_flag(self):
        assignment = Assignment("'ifort -O3 -march=native'")
        self.interpreter._interpret_fortran_comp(assignment, skip_check=True)
        assert 'ifort -O3 -march=native' in self.rpars.FORTRAN_COMP[0]

    def test__fortran_comp_custom_post_flag(self):
        assignment = Assignment("'-L/opt/intel/mkl/lib/intel64'", flags='post')
        self.interpreter._interpret_fortran_comp(assignment, skip_check=True)
        assert '-L/opt/intel/mkl/lib/intel64' in self.rpars.FORTRAN_COMP[1]

    def test__fortran_comp_custom_mpi_flag(self):
        assignment = Assignment("'mpifort -fallow-argument-mismatch'", flags='mpi')
        self.interpreter._interpret_fortran_comp(assignment, skip_check=True)
        assert 'mpifort -fallow-argument-mismatch' in self.rpars.FORTRAN_COMP_MPI[0]

# _interpret_v0_real()
class TestV0Real():
    rpars = Rparams()
    interpreter = ParameterInterpreter(rpars, slab=None)

    def test__interpret_v0_real_rundgren_type(self):
        assignment = Assignment("rundgren 1.0 2.0 3.0 4.0")
        self.interpreter._interpret_v0_real(assignment)
        assert self.rpars.V0_REAL == [1.0, 2.0, 3.0, 4.0]

    def test__interpret_v0_real_other_type(self,):
        assignment = Assignment("EE")
        self.interpreter._interpret_v0_real(assignment)
        assert self.rpars.V0_REAL == "EEV+workfn"

    def test__interpret_v0_real_invalid_rundgren_constants(self):
        assignment = Assignment("rundgren 1.0 2.0")
        with pytest.raises(ParameterError):
            self.interpreter._interpret_v0_real(assignment)

    def test__interpret_v0_real_value_error(self):
        assignment = Assignment("rundgren 1.0 2.0 3.0 four")
        with pytest.raises(ParameterError):
            self.interpreter._interpret_v0_real(assignment)

class TestTheoEnergies:
    rpars = Rparams()

    def test__interpret_theo_energies_single_value(self):
        interpreter = ParameterInterpreter(self.rpars)
        assignment = Assignment("1.0")
        interpreter._interpret_theo_energies(assignment)
        assert self.rpars.THEO_ENERGIES == [1.0, 1.0, 1]

    def test__interpret_theo_energies_multiple_values(self):
        interpreter = ParameterInterpreter(self.rpars)
        assignment = Assignment("1.0 _ 2.0")
        interpreter._interpret_theo_energies(assignment)
        assert self.rpars.THEO_ENERGIES == [1.0, -1, 2.0]

    def test__interpret_theo_energies_invalid_float(self):
        interpreter = ParameterInterpreter(self.rpars)
        assignment = Assignment("invalid")
        with pytest.raises(ParameterFloatConversionError):
            interpreter._interpret_theo_energies(assignment)

    def test__interpret_theo_energies_invalid_positive_value(self):
        interpreter = ParameterInterpreter(self.rpars)
        assignment = Assignment("1.0 -2.0 3.0")
        with pytest.raises(ParameterError):
            interpreter._interpret_theo_energies(assignment)

    def test__interpret_theo_energies_invalid_number_of_values(self):
        interpreter = ParameterInterpreter(self.rpars)
        assignment = Assignment("1.0 2.0")
        with pytest.raises(ParameterNumberOfInputsError):
            interpreter._interpret_theo_energies(assignment)

class TestNumericalParamsExamples:
    def test__interpret_n_bulk_layers_valid(self):
        rpars = Rparams()
        interpreter = ParameterInterpreter(rpars, slab=None)
        assignment = Assignment("1")
        assert rpars.N_BULK_LAYERS == 1

    def test__interpret_n_bulk_layers_invalid(self):
        # N_BULK_LAYERS must be 1 or 2
        rpars = Rparams()
        interpreter = ParameterInterpreter(rpars, slab=None)
        assignment = Assignment("3")
        with pytest.raises(ParameterError):
            interpreter._interpret_n_bulk_layers(assignment)


    def test__interpret_t_debye(self):
        rpars = Rparams()
        interpreter = ParameterInterpreter(rpars, slab=None)
        assignment = Assignment("300.0")
        interpreter._interpret_t_debye(assignment)
        assert rpars.T_DEBYE == pytest.approx(300.0)

class TestBoolParamsExamples:
    def test__interpret_log_debug_true(self):
        rpars = Rparams()
        interpreter = ParameterInterpreter(rpars, slab=None)
        assignment = Assignment("true")
        interpreter._interpret_log_debug(assignment)
        assert rpars.LOG_DEBUG is True

    def test__interpret_log_debug_false(self):
        rpars = Rparams()
        interpreter = ParameterInterpreter(rpars, slab=None)
        assignment = Assignment("F")
        interpreter._interpret_log_debug(assignment)
        assert rpars.LOG_DEBUG is False

class TestFilamentWF:
    def test__interpret_filament_wf_lab6(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams, slab=None)
        assignment = Assignment("LaB6")
        interpreter._interpret_filament_wf(assignment)
        assert mock_rparams.FILAMENT_WF == pytest.approx(2.65)

    def test__interpret_filament_wf_custom(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams, slab=None)
        assignment = Assignment("1.0")
        interpreter._interpret_filament_wf(assignment)
        assert mock_rparams.FILAMENT_WF == pytest.approx(1.0)

    def test__interpret_filament_wf_invalid(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams, slab=None)
        assignment = Assignment("invalid")
        with pytest.raises(ParameterFloatConversionError):
            interpreter._interpret_filament_wf(assignment)

    def test__interpret_filament_wf_flag_invalid(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams, slab=None)
        assignment = Assignment("1.5", flags='test')
        with pytest.raises(ParameterUnknownFlagError):
            interpreter._interpret_filament_wf(assignment)

class TestAverageBeams:
    def test__interpret_average_beams_off(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("off")
        interpreter._interpret_average_beams(assignment)
        assert mock_rparams.AVERAGE_BEAMS is False

    def test__interpret_average_beams_all(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("all")
        interpreter._interpret_average_beams(assignment)
        assert mock_rparams.AVERAGE_BEAMS == (0.0, 0.0)

    def test__interpret_average_beams_custom_angles(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("45.0 60.0")
        interpreter._interpret_average_beams(assignment)
        assert mock_rparams.AVERAGE_BEAMS == (45.0, 60.0)

    def test__interpret_average_beams_invalid_input(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("invalid input")
        with pytest.raises(ParameterError):
            interpreter._interpret_average_beams(assignment)

class TestBeamIncidence:
    def test__interpret_beam_incidence_custom_angles(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("45.0 60.0")
        interpreter._interpret_beam_incidence(assignment)
        assert mock_rparams.THETA == 45.0
        assert mock_rparams.PHI == 60.0

    def test__interpret_beam_incidence_invalid_input(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("invalid input")
        with pytest.raises(ParameterError):
            interpreter._interpret_beam_incidence(assignment)

class TestDomain:
    def test__interpret_domain_existing_path(self, mock_rparams, tmp_path):
        interpreter = ParameterInterpreter(mock_rparams)
        domain_path = tmp_path / "domain1"
        domain_path.mkdir()
        assignment = Assignment(flags="domain1", right_side=str(domain_path))
        interpreter._interpret_domain(assignment)
        assert mock_rparams.DOMAINS == [("domain1", str(domain_path))]

    def test__interpret_domain_new_path(self, mock_rparams, tmp_path):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(str(tmp_path))
        interpreter._interpret_domain(assignment)
        assert mock_rparams.DOMAINS == [("1", str(tmp_path))]

    def test__interpret_domain_zip_file(self, mock_rparams, tmp_path):
        zip_file = tmp_path / "domain.zip"
        zip_file.touch()
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(str(zip_file))
        interpreter._interpret_domain(assignment)
        assert mock_rparams.DOMAINS == [("1", str(zip_file))]

    def test__interpret_domain_invalid_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("invalid_path")
        with pytest.raises(ParameterError):
            interpreter._interpret_domain(assignment)

class TestDomainStep:
    def test__interpret_domain_step_valid_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("10")
        interpreter._interpret_domain_step(assignment)
        assert mock_rparams.DOMAIN_STEP == 10

    def test__interpret_domain_step_invalid_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("200")
        with pytest.raises(ParameterRangeError):
            interpreter._interpret_domain_step(assignment)

    def test__interpret_domain_step_non_integer_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0.5")
        with pytest.raises(ParameterError):
            interpreter._interpret_domain_step(assignment)

    def test__interpret_domain_step_non_divisible_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("25")
        with pytest.raises(ParameterError):
            interpreter._interpret_domain_step(assignment)

class TestSuperlattice:
    def test__interpret_superlattice_matrix_notation(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment(flags="M", right_side="2 0, 0 2")
        interpreter._interpret_superlattice(assignment)
        assert mock_rparams.SUPERLATTICE == pytest.approx(np.array([[2, 0], [0, 2,]]))
        assert mock_rparams.superlattice_defined is True

    def test__interpret_superlattice_woods_notation(self, mock_rparams, slab_ag100):
        interpreter = ParameterInterpreter(mock_rparams, slab=slab_ag100)
        assignment = Assignment("p(2x1)")
        interpreter._interpret_superlattice(assignment)
        assert mock_rparams.SUPERLATTICE == pytest.approx(np.array([[2, 0], [0, 1,]]))
        assert mock_rparams.superlattice_defined is True

    def test__interpret_superlattice_no_slab(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams, slab=None)
        assignment = Assignment("c(2x2)")
        # should raise because slab is None
        with pytest.raises(ParameterError):
            interpreter._interpret_superlattice(assignment)

class TestTensorOutput:
    def test__interpret_tensor_output_single_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("False")
        interpreter._interpret_tensor_output(assignment)
        assert mock_rparams.TENSOR_OUTPUT == [0]

    def test__interpret_tensor_output_multiple_values(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("0 1 1 0")
        interpreter._interpret_tensor_output(assignment)
        assert mock_rparams.TENSOR_OUTPUT == [0, 1, 1, 0]

    def test__interpret_tensor_output_repeated_values(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("2*1 2*0 1")
        interpreter._interpret_tensor_output(assignment)
        assert mock_rparams.TENSOR_OUTPUT == [1, 1, 0, 0, 1]

    def test__interpret_tensor_output_invalid_value(self, mock_rparams):
        interpreter = ParameterInterpreter(mock_rparams)
        assignment = Assignment("2")
        with pytest.raises(ParameterParseError):
            interpreter._interpret_tensor_output(assignment)
