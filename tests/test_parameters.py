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
from viperleed.tleedmlib.files.parameter_errors import ParameterError

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

class TestAg100Parameters():
    def test_read_parameters_for_ag100(self):
        # just check that readPARAMETERS does not crash; not interpreted yet
        filename = 'tests/fixtures/parameters/PARAMETERS_1'
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


# unit tests for parameter interpretation

# _interpret_intpol_deg()
def test_interpret_intpol_deg():
    rpars = Rparams()
    for val in rpars.get_limits('INTPOL_DEG'):
        interpreter = ParameterInterpreter(rpars, slab=None)
        interpreter._interpret_intpol_deg(Assignment(value=val))
        assert rpars.INTPOL_DEG == int(val)
    incompatible_values = ['1', 'text']
    for val in incompatible_values:
        with pytest.raises(ParameterError):
                    interpreter = ParameterInterpreter(rpars, slab=None)
                    interpreter._interpret_intpol_deg(Assignment(value=val))

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
                interpreter._interpret_bulk_repeat(Assignment(value=val, right_side=val))


    def test__interpret_bulk_repeat_float(self, slab_ag100):
        val = '1.428'
        interpreter = ParameterInterpreter(self.rpars, slab=slab_ag100)
        interpreter._interpret_bulk_repeat(Assignment(value=val, right_side=val))
        assert self.rpars.BULK_REPEAT == pytest.approx(1.428, rel=1e-4)

    def test__interpret_bulk_repeat_c(self, slab_ag100):
        val = 'c(0.1)'
        interpreter = ParameterInterpreter(self.rpars, slab=slab_ag100)
        interpreter._interpret_bulk_repeat(Assignment(value=val, right_side=val))
        assert self.rpars.BULK_REPEAT == pytest.approx(2.03646, rel=1e-4)

    def test__interpret_bulk_repeat_z(self, slab_ag100):
        val = 'c(0.1)'
        interpreter = ParameterInterpreter(self.rpars, slab=slab_ag100)
        interpreter._interpret_bulk_repeat(Assignment(value=val, right_side=val))
        assert self.rpars.BULK_REPEAT == pytest.approx(2.0364, rel=1e-4)

    def test__interpret_bulk_repeat_vector(self, slab_ag100):
        val = '[1.0 2.0 3.0]'
        interpreter = ParameterInterpreter(self.rpars, slab=slab_ag100)
        interpreter._interpret_bulk_repeat(Assignment(value=val, right_side=val))
        assert self.rpars.BULK_REPEAT == pytest.approx([1.0, 2.0, 3.0], rel=1e-4)

# _interpret_fortran_comp()
class TestInterpretFortranComp():
    # TODO: make use of new intepreter class
    param = 'FORTRAN_COMP'
    rpars = Rparams()

    def test_fortran_comp_default_intel(self):
        flags, other_values, right_side = [], [], ''
        value = 'ifort'
        parameters._interpret_fortran_comp(self.rpars, self.param, flags, right_side, value, other_values, skip_check=True)
        assert 'ifort -O2 -I/opt/intel/mkl/include' in self.rpars.FORTRAN_COMP[0]
        assert '-L/opt/intel/mkl/lib/intel64' in self.rpars.FORTRAN_COMP[1]

    def test_fortran_comp_default_gnu(self):
        flags, other_values, right_side = [], [], ''
        value = 'gfortran'
        parameters._interpret_fortran_comp(self.rpars, self.param, flags, right_side, value, other_values, skip_check=True)
        assert 'gfortran' in self.rpars.FORTRAN_COMP[0]
        assert '-llapack -lpthread -lblas' in self.rpars.FORTRAN_COMP[1]
    
    def test_fortran_comp_default_intel_mpi(self):
        flags, other_values, right_side = ['mpi'], [], ''
        value = 'mpiifort'
        parameters._interpret_fortran_comp(self.rpars, self.param, flags, right_side, value, other_values, skip_check=True)
        assert 'mpiifort' in self.rpars.FORTRAN_COMP_MPI[0]  # no other flags by default

    def test_fortran_comp_default_gnu_mpi(self):
        flags, other_values, right_side = ['mpi'], [], ''
        value = 'mpifort'
        parameters._interpret_fortran_comp(self.rpars, self.param, flags, right_side, value, other_values, skip_check=True)
        assert 'mpifort -Ofast -no-pie' in self.rpars.FORTRAN_COMP_MPI[0]  # no other flags by default

    def test_fortran_comp_custom_without_flag(self):
        flags, other_values, right_side = [], [],  "'ifort -O3 -march=native'"
        value = "'ifort -O3 -march=native'"
        parameters._interpret_fortran_comp(self.rpars, self.param, flags, right_side, value, other_values, skip_check=True)
        assert 'ifort -O3 -march=native' in self.rpars.FORTRAN_COMP[0]

    def test_fortran_comp_custom_post_flag(self):
        flags, other_values, right_side = ['post'], [], "'-L/opt/intel/mkl/lib/intel64'"
        value = "'-L/opt/intel/mkl/lib/intel64'"
        parameters._interpret_fortran_comp(self.rpars, self.param, flags, right_side, value, other_values, skip_check=True)
        assert '-L/opt/intel/mkl/lib/intel64' in self.rpars.FORTRAN_COMP[1]

    def test_fortran_comp_custom_mpi_flag(self):
        flags, other_values, right_side = ['mpi'], [], "'mpifort -fallow-argument-mismatch'"
        value = "'mpifort -fallow-argument-mismatch'"
        parameters._interpret_fortran_comp(self.rpars, self.param, flags, right_side, value, other_values, skip_check=True)
        assert 'mpifort -fallow-argument-mismatch' in self.rpars.FORTRAN_COMP_MPI[0]

    # test that error is raised if unsupported default compiler is used
    def test_fortran_wrong_compiler(self):
        flags, other_values, right_side = [], [], ''
        value = 'clang'
        with pytest.raises(ParameterError):
            parameters._interpret_fortran_comp(self.rpars, self.param, flags, right_side, value, other_values, skip_check=True)

