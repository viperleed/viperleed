"""Tests for module viperleed.tleedmlib.files.parameters.

Created on 2023-06-09

@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)
"""

from pathlib import Path
import sys

import numpy as np
import pytest
from pytest_cases import parametrize_with_cases

VPR_PATH = str(Path(__file__).resolve().parents[3])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Will be fixed in installable version
from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tleedmlib.files import parameters
from viperleed.tleedmlib.files.parameters import errors as err
from viperleed.tleedmlib.files.parameters._utils import Assignment
from viperleed.tleedmlib.files.parameters._utils import NumericBounds as Bounds

from .case_parameters import case_parameters_slab
# pylint: enable=wrong-import-position


class TestSlabParameters:
    """Tests for successful read/interpretation of slab-related parameters."""

    @staticmethod
    def check_sitedef_consistent(site_def, expected, subtests):
        """Ensure site_def has the expected values."""
        # The tricky bit is that we currently store site_def values
        # as lists, which includes order. Here we do an unordered
        # comparison.
        for element, sites in expected.items():
            for tag, atom_nrs in sites.items():
                with subtests.test(f'SITE_DEF[{element}][{tag}]'):
                    assert set(site_def[element][tag]) == set(atom_nrs)

    @parametrize_with_cases('args', cases=case_parameters_slab)
    def test_parameters_interpreted(self, args, subtests):
        """Check that parameters have been interpreted correctly."""
        slab, rpars, info = args
        interpreter = parameters.ParameterInterpreter(rpars)
        interpreter.interpret(slab)
        expected = info.parameters.expected
        for name, value in expected.items():
            attr = getattr(rpars, name)
            if name == 'SITE_DEF':
                self.check_sitedef_consistent(attr, value, subtests)
                continue
            with subtests.test(name):
                if isinstance(value, np.ndarray):
                    assert attr == pytest.approx(value)
                else:
                    assert attr == value


# TODO: be more specific than ParameterError!
# TODO: would be nice to reduce a bit the verbosity to use (maybe)
# a metaclass in order to implement test_interpret_(in)valid once,
# and parametrize it with the valid/invalid of each class below.
# Would it be cleaner with CASES?


@pytest.fixture(name='interpreter')
def fixture_interpreter():
    """Return a fresh ParameterInterpreter."""
    return parameters.ParameterInterpreter(Rparams())


class _TestInterpretBase:
    """Base class for parameter-interpretation tests."""
    param = None
    rpars_attr = None

    def assignment(self, value_str, **kwargs):
        """Return an Assignment object for self.param."""
        return Assignment(value_str, self.param, **kwargs)

    def interpret(self, interpreter, value_str, **kwargs):
        """Interpret a value for self.param."""
        method = getattr(interpreter, f'interpret_{self.param.lower()}')
        return method(self.assignment(value_str, **kwargs))

    def check_assigned(self, interpreter, value_str, expected, **kwargs):
        """Assert interpretation is successful."""
        self.interpret(interpreter, value_str, **kwargs)
        attr = getattr(interpreter.rpars, self.rpars_attr or self.param)
        assert attr == pytest.approx(expected)

    def check_raises(self, interpreter, value_str, exc, **kwargs):
        """Assert that an attempt to interpret raises an exc."""
        with pytest.raises(exc):
            self.interpret(interpreter, value_str, **kwargs)


class TestNumericalParameter(_TestInterpretBase):
    """Tests for numeric-parameter interpretation."""

    param = 'TEST_PARAM'
    valid = {  # value, bounds, expected
        'float': ('0.01', Bounds(), 0.01),
        'int': ('100', Bounds(type_=int), 100),
        'int at bound': ('100', Bounds(type_=int, range_=(0, 100)), 100),
        'float at bound': ('-0.5', Bounds(range_=(-0.5, 1)), -0.5),
        'float modulo': ('2.5',
                         Bounds(range_=(0, 2), out_of_range_event='modulo'),
                         0.5),
        'int modulo': ('12',
                       Bounds(type_=int, range_=(0, 10),
                              out_of_range_event='modulo'),
                       2),
        }
    invalid = {  # value, bounds
        'int at bound': ('100',
                         Bounds(type_=int, range_=(0, 100),
                                accept_limits=(False, False)),),
        'float at bound': ('-0.5',
                           Bounds(range_=(-0.5, 1),
                                  accept_limits=(False, False))),
        'float range': ('-0.5', Bounds(range_=(0, 1))),
        'int range': ('-5', Bounds(type_=int, range_=(0, 10))),
        }

    # pylint: disable=too-many-arguments
    # 6/5 seems OK here, considering that two are fixtures.
    @pytest.mark.parametrize('val,bounds,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, bounds, expect, interpreter, subtests):
        """Ensure valid values are returned and assigned."""
        result = interpreter.interpret_numerical_parameter(
            self.assignment(val),
            bounds=bounds
            )
        with subtests.test('result'):
            assert result == pytest.approx(expect)
        with subtests.test('attribute'):
            assert interpreter.rpars.TEST_PARAM == pytest.approx(expect)
    # pylint: enable=too-many-arguments

    @pytest.mark.parametrize('val,bounds', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, bounds, interpreter):
        """Ensure exceptions are raised for invalid combinations."""
        with pytest.raises(err.ParameterError):
            interpreter.interpret_numerical_parameter(self.assignment(val),
                                                      bounds=bounds)


class TestSimpleParamsExamples:
    """Tests for a selection of simple numeric PARAMETERS."""

    def test_interpret_n_bulk_layers_valid(self, interpreter):
        """Check assignment of valid N_BULK_LAYERS."""
        assignment = Assignment('1', 'N_BULK_LAYERS')
        interpreter.interpret_n_bulk_layers(assignment)
        assert interpreter.rpars.N_BULK_LAYERS == 1

    def test_interpret_n_bulk_layers_invalid(self, interpreter):
        """Check exceptions are raised for invalid N_BULK_LAYERS."""
        # N_BULK_LAYERS must be 1 or 2
        assignment = Assignment('3', 'N_BULK_LAYERS')
        with pytest.raises(err.ParameterError):
            interpreter.interpret_n_bulk_layers(assignment)

    def test_interpret_t_debye(self, interpreter):
        """Check assignment of valid T_DEBYE."""
        assignment = Assignment('300.0', 'T_DEBYE')
        interpreter.interpret_t_debye(assignment)
        assert interpreter.rpars.T_DEBYE == pytest.approx(300.0)


class TestAverageBeams(_TestInterpretBase):
    """Tests for interpreting AVERAGE_BEAMS."""

    param = 'AVERAGE_BEAMS'
    valid = {'off': ('off', False),
             'all': ('all', (0.0, 0.0)),
             'custom': ('45.0 60', (45.0, 60.0)),}
    invalid = ('invalid input',)

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid AVERAGE_BEAMS."""
        self.check_assigned(interpreter, val, expect)

    @pytest.mark.parametrize('val', invalid, ids=invalid)
    def test_interpret_invalid(self, val, interpreter):
        """Ensure invalid AVERAGE_BEAMS raises exceptions."""
        self.check_raises(interpreter, val, err.ParameterError)


class TestBeamIncidence(_TestInterpretBase):
    """Tests for interpreting BEAM_INCIDENCE."""

    param = 'BEAM_INCIDENCE'
    valid = {'custom': ('45.0 60', (45.0, 60.0)),}
    invalid = ('invalid input',)

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check interpretation of valid BEAM_INCIDENCE into THETA and PHI."""
        self.interpret(interpreter, val)
        rpars = interpreter.rpars
        assert (rpars.THETA, rpars.PHI) == expect

    @pytest.mark.parametrize('val', invalid, ids=invalid)
    def test_interpret_invalid(self, val, interpreter):
        """Ensure invalid BEAM_INCIDENCE raises exceptions."""
        self.check_raises(interpreter, val, err.ParameterError)


class TestBulkRepeat(_TestInterpretBase):
    """Tests for interpreting BULK_REPEAT."""

    param = 'BULK_REPEAT'
    valid = {  # (input, expected)
        'float': ('1.428', 1.428),
        'c_small': ('c(0.1)', 2.0364),
        'c_large': ('c(1.2)', 12 * 2.0364),
        'z': ('z(0.9)', 0.9),
        'vector': ('[1.0 2.0 3.0]', [1.0, 2.0, 3.0])
        }
    invalid = ('text', 'y(1.2)', 'z(abc)', '[]')

    @pytest.fixture(name='ag100_interpreter')
    def fixture_ag100_interpreter(self, ag100, interpreter):
        """Return a ParameterInterpreter for a Ag(100) slab."""
        interpreter.slab, *_ = ag100
        return interpreter

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, ag100_interpreter):
        """Check correct interpretation of valid BULK_REPEAT."""
        self.interpret(ag100_interpreter, val)
        rpars = ag100_interpreter.rpars
        assert rpars.BULK_REPEAT == pytest.approx(expect, rel=1e-4)

    @pytest.mark.parametrize('val', invalid)
    def test_interpret_invalid(self, val, ag100_interpreter):
        """Ensure invalid BULK_REPEAT raises exceptions."""
        self.check_raises(ag100_interpreter, val, err.ParameterError)


class TestDomain(_TestInterpretBase):
    """Tests for interpreting DOMAIN."""

    param = 'DOMAIN'

    def test_interpret_path_with_flag(self, interpreter, tmp_path):
        """Test correct interpretation of a path with a domain name."""
        domain_path = tmp_path / 'domain1'
        domain_path.mkdir()
        domain_path = str(domain_path)
        self.interpret(interpreter, domain_path, flags_str='domain1')
        assert interpreter.rpars.DOMAINS == [('domain1', domain_path)]

    def test_interpret_path_no_flag(self, interpreter, tmp_path):
        """Test correct interpretation of a path without a domain name."""
        tmp_path = str(tmp_path)
        self.interpret(interpreter, tmp_path)
        assert interpreter.rpars.DOMAINS == [('1', tmp_path)]

    def test_interpret_zip_file(self, interpreter, tmp_path):
        """Test correct interpretation of a zip file."""
        zip_file = tmp_path / 'domain.zip'
        zip_file.touch()
        zip_file = str(zip_file)
        self.interpret(interpreter, zip_file)
        assert interpreter.rpars.DOMAINS == [('1', zip_file)]

    def test_interpret_invalid(self, interpreter):
        """Ensure invalid DOMAIN raises exceptions."""
        self.check_raises(interpreter, 'invalid_path', err.ParameterError)


class TestDomainStep(_TestInterpretBase):
    """Tests for interpreting DOMAIN_STEP."""

    param = 'DOMAIN_STEP'
    valid = {'value': ('10', 10),}
    invalid = {'value': ('200', err.ParameterRangeError),
               'non_integer': ('0.5', err.ParameterError),}

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid DOMAIN_STEP."""
        self.check_assigned(interpreter, val, expect)

    @pytest.mark.parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid DOMAIN_STEP raises exceptions."""
        self.check_raises(interpreter, val, exc)


class TestElementRename(_TestInterpretBase):
    """Tests for interpreting ELEMENT_RENAME."""

    param = 'ELEMENT_RENAME'

    def test_interpret_valid(self, interpreter):
        """Check correct interpretation of valid ELEMENT_RENAME."""
        self.interpret(interpreter, 'H', flags_str='X')
        assert interpreter.rpars.ELEMENT_RENAME == {'X': 'H'}

    def test_interpret_invalid_element(self, interpreter):
        """Ensure that an invalid chemical element raises exceptions."""
        self.check_raises(interpreter, 'Op', err.ParameterError,
                          flags_str='Uk')


class TestFilamentWF(_TestInterpretBase):
    """Tests for interpreting FILAMENT_WF."""

    param = 'FILAMENT_WF'
    valid = {'lab6': ('LaB6', 2.65),
             'custom': ('1.0', 1.0),}
    invalid = {
        'invalid_float': ('invalid', '', err.ParameterFloatConversionError),
        'flag': ('1.5', 'test', err.ParameterUnknownFlagError),
        }

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid FILAMENT_WF."""
        self.check_assigned(interpreter, val, expect)

    @pytest.mark.parametrize('val,flag,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, flag, exc, interpreter):
        """Ensure invalid FILAMENT_WF raises exceptions."""
        self.check_raises(interpreter, val, exc, flags_str=flag)


class TestFortranComp(_TestInterpretBase):
    """Tests for interpreting FORTRAN_COMP."""

    param = 'FORTRAN_COMP'
    valid = {  # id: (assign_value, assign_flag, expect_pre, expect_post)
        'default_intel': ('ifort', '',
                          'ifort -O2 -I/opt/intel/mkl/include',
                          '-L/opt/intel/mkl/lib/intel64'),
        'default_gnu': ('gfortran', '', 'gfortran',
                        '-llapack -lpthread -lblas'),
        'default_intel_mpi': ('mpiifort', 'mpi', 'mpiifort', None),
        'default_gnu_mpi': ('mpifort', 'mpi', 'mpifort -Ofast -no-pie', None),
        'custom_no_flag': ('"ifort -O3 -march=native"', '',
                           'ifort -O3 -march=native', None),
        'custom_post_flag': ('"-L/opt/intel/mkl/lib/intel64"', 'post',
                             None, '-L/opt/intel/mkl/lib/intel64'),
        'custom_mpi_flag': ('"mpifort -fallow-argument-mismatch"', 'mpi',
                            'mpifort -fallow-argument-mismatch', None)
        }

    # pylint: disable=too-many-arguments
    # In principle 'pre' and 'post' could be merged into a tuple,
    # but the parametrization above would look even more complex
    @pytest.mark.parametrize('val,flag,pre,post', valid.values(), ids=valid)
    def test_interpret_valid(self, val, flag, pre, post,
                             interpreter, subtests):
        """Check correct interpretation of valid FORTRAN_COMP(_MPI)."""
        assignment = self.assignment(val, flags_str=flag)
        rpars = interpreter.rpars
        interpreter.interpret_fortran_comp(assignment, skip_check=True)
        if 'mpi' in flag:
            compiler = getattr(rpars, self.param + '_MPI')
        else:
            compiler = getattr(rpars, self.param)
        if pre is not None:
            with subtests.test('Check pre'):
                assert pre in compiler[0]
        if post is not None:
            with subtests.test('Check post'):
                assert post in compiler[1]
    # pylint: enable=too-many-arguments


class TestIntpolDeg(_TestInterpretBase):
    """Tests for interpreting INTPOL_DEG."""

    param = 'INTPOL_DEG'
    valid = {v: (v, int(v)) for v in Rparams().get_limits(param)}
    invalid = '1', 'text'

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid INTPOL_DEG."""
        self.check_assigned(interpreter, val, expect)

    @pytest.mark.parametrize('val', invalid, ids=invalid)
    def test_interpret_invalid(self, val, interpreter):
        """Ensure invalid INTPOL_DEG raises exceptions."""
        self.check_raises(interpreter, val, err.ParameterError)


class TestIVShiftRange(_TestInterpretBase):
    """Tests for interpreting IV_SHIFT_RANGE."""

    param = 'IV_SHIFT_RANGE'
    valid = {'range': ('0.0 1.0 0.25', [0.0, 1.0, 0.25]),}
    invalid = {'nr_inputs': ('0.0', err.ParameterNumberOfInputsError),
               'float': ('0.0 2.0 ()', err.ParameterFloatConversionError),
               'parse': ('0.0 2.0 1.a', err.ParameterParseError),
               'out-of-range': ('1.0 0.0 0.1', err.ParameterError),
               'step': ('0.0 1.0 -0.1', err.ParameterError),}

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid IV_SHIFT_RANGE."""
        self.check_assigned(interpreter, val, expect)

    @pytest.mark.parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid IV_SHIFT_RANGE raises exceptions."""
        self.check_raises(interpreter, val, exc)


class TestLayerCuts(_TestInterpretBase):
    """Tests for interpreting LAYER_CUTS."""

    param = 'LAYER_CUTS'
    valid = {
        'simple': ('0.1 0.2 0.3', ['0.1', '0.2', '0.3']),
        'dz': ('dz(1.2)', ['dz(1.2)']),
        'list and dz': ('0.15 0.3 < dz(1.2)  ',
                        ['0.15', '0.3', '<', 'dz(1.2)']),
        'two dz': ('dz(1.0) < 2.0 < dz(0.5) < 4.0',
                   ['dz(1.0)', '<', '2.0', '<', 'dz(0.5)', '<', '4.0']),
        }
    invalid = {
        'less and greater': '< 0.1 > 0.2',
        'cutoff function': 'dz(0.1) dc(0.2) invalid(0.3)',
        'float': '0.1 invalid 0.3',
        'dz': '0.5 1.0 < dz(abcd) < 4.0'
        }

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid LAYER_CUTS."""
        if val.count('dz') == 2:
            pytest.xfail(reason='Known bug in LAYER_CUTS interpreter')
        self.check_assigned(interpreter, val, expect)

    @pytest.mark.parametrize('val', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, interpreter):
        """Ensure invalid LAYER_CUTS raises exceptions."""
        if 'abcd' in val:
            pytest.xfail(reason='Known bug in LAYER_CUTS interpreter')
        self.check_raises(interpreter, val, err.ParameterParseError)


class TestLogLevel(_TestInterpretBase):
    """Tests for interpreting LOG_LEVEL."""

    param = 'LOG_LEVEL'

    def test_interpret_true(self, interpreter):
        """Check correct setting of LOG_LEVEL to True."""
        self.interpret(interpreter, 'true')
        assert interpreter.rpars.LOG_LEVEL <= 10

    def test_interpret_false(self, interpreter):
        """Check correct setting of LOG_LEVEL to False."""
        self.interpret(interpreter, 'F')
        assert interpreter.rpars.LOG_LEVEL >= 20

    def test_interpret_int(self, interpreter):
        """Check correct interpretation of integer LOG_LEVEL."""
        self.check_assigned(interpreter, '3', 3)


class TestOptimize(_TestInterpretBase):
    """Tests for interpreting OPTIMIZE."""

    param = 'OPTIMIZE'
    valid = {  # value, flag, {key: expected}
        'flag_and_value': ('step 0.1', 'v0i', {'step': 0.1}),
        'multiple_flag_value': (
            'step 0.1, convergence 1e-6, minpoints 10', 'theta',
            {'step': 0.1, 'minpoints': 10, 'convergence': 1e-6}
            ),
        }
    invalid = {
        'missing_quantity': ('step 0.1', '', err.ParameterNeedsFlagError),
        'flag': ('invalid 0.1', 'v0i', err.ParameterUnknownFlagError),
        'value': ('step not-a-number', 'v0i', err.ParameterError),
        'quantity': ('step 0.1', 'invalid', err.ParameterError),
        }

    # pylint: disable=too-many-arguments
    # 6/5 Seems OK here, especially considering that two are fixtures
    @pytest.mark.parametrize('val,flag,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, flag, expect, interpreter, subtests):
        """Check correct interpretation of valid OPTIMIZE."""
        self.interpret(interpreter, val, flags_str=flag)
        rpars = interpreter.rpars
        with subtests.test('which'):
            assert rpars.OPTIMIZE['which'] == flag
        for key, value in expect.items():
            with subtests.test(key):
                assert rpars.OPTIMIZE[key] == pytest.approx(value)
    # pylint: enable=too-many-arguments

    @pytest.mark.parametrize('val,flag,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, flag, exc, interpreter):
        """Ensure invalid OPTIMIZE raises exceptions."""
        self.check_raises(interpreter, val, exc, flags_str=flag)


class TestPhaseshiftEps(_TestInterpretBase):
    """Tests for interpreting PHASESHIFT_EPS."""

    param = 'PHASESHIFT_EPS'
    valid = {'float': ('0.1', 0.1),
             'tag': ('fine', 0.01),}
    invalid = {'float': ('invalid', err.ParameterError),
               'negative': ('-1.0', err.ParameterRangeError),
               'too large': ('1.5', err.ParameterRangeError),}

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid PHASESHIFT_EPS."""
        self.check_assigned(interpreter, val, expect)

    @pytest.mark.parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid PHASESHIFT_EPS raises exceptions."""
        self.check_raises(interpreter, val, exc)


class TestRun(_TestInterpretBase):
    """Tests for interpreting RUN."""

    param = 'RUN'
    valid = {
        'single': ('1', [0, 1]),
        'multiple': ('1 2 3', [0, 1, 2, 3]),
        'range': ('1-3', [0, 1, 2, 3]),
        }
    invalid = {'section': ('invalid', err.ParameterValueError),
               'empty': ('', err.ParameterError),
               'syntax underscore': ('1 _ 2', err.ParameterValueError),
               'syntax section': ('1, x', err.ParameterValueError)}

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid RUN."""
        self.check_assigned(interpreter, val, expect)

    @pytest.mark.parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid RUN raises exceptions."""
        self.check_raises(interpreter, val, exc)


class TestSearchBeams(_TestInterpretBase):
    """Tests for interpreting SEARCH_BEAMS."""

    param = 'SEARCH_BEAMS'
    valid = {'average': ('A', 0), 'alt average': ('0', 0),
             'integer': ('I', 1), 'alt integer': ('1', 1),
             'fractional': ('F', 2), 'alt fractional': ('2', 2)}
    invalid = {'one value': ('3', err.ParameterValueError),
               'more values': ('0 1', err.ParameterError)}

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid SEARCH_BEAMS."""
        self.check_assigned(interpreter, val, expect)

    @pytest.mark.parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid SEARCH_BEAMS raises exceptions."""
        self.check_raises(interpreter, val, exc)


class TestSearchConvergence(_TestInterpretBase):
    """Tests for interpreting SEARCH_CONVERGENCE."""

    param = 'SEARCH_CONVERGENCE'
    valid = {'gaussian, scaling': ('0.01 0.9', (0.01, 0.9)),
             'gaussian, no scaling': ('0.01', (0.01, 0.5)),}
    invalid = {
        'flag': (' 0.01 0.9', 'test', err.ParameterUnknownFlagError),
        'nr_values': ('.01 0.9 0.5', 'gaussian',
                      err.ParameterNumberOfInputsError),
        'scaling': ('0.01 -0.5', 'gaussian', err.ParameterError),
        }

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid SEARCH_CONVERGENCE."""
        self.interpret(interpreter, val, flags_str='gaussian')
        rpars = interpreter.rpars
        assert (rpars.GAUSSIAN_WIDTH, rpars.GAUSSIAN_WIDTH_SCALING) == expect

    @pytest.mark.parametrize('val,flag,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, flag, exc, interpreter):
        """Ensure invalid SEARCH_CONVERGENCE raises exceptions."""
        self.check_raises(interpreter, val, exc, flags_str=flag)


class TestSearchCull(_TestInterpretBase):
    """Tests for interpreting SEARCH_CULL and SEARCH_CULL_TYPE."""

    param = 'SEARCH_CULL'
    valid = {'float': ('0.5', 0.5), 'int': ('3', 3),
             'int-like float': ('4.0', 4),
             'explicit type': ('0.48 clone', (0.48, 'clone'))}
    invalid = {
        'float greater than one': ('1.5', err.ParameterValueError),
        'negative': ('-0.5', err.ParameterValueError),
        'too many': ('0.5 clone 1.0', err.ParameterNumberOfInputsError),
        'invalid cull type': ('0.5 test', err.ParameterValueError),
        'not a float': ('abcd', err.ParameterFloatConversionError),
        }

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid SEARCH_CULL."""
        rpars = interpreter.rpars
        if len(val.split()) > 1:
            expect_cull, expect_type = expect
        else:
            expect_cull = expect
            expect_type = rpars.get_default('SEARCH_CULL_TYPE')
        self.check_assigned(interpreter, val, expect_cull)
        assert rpars.SEARCH_CULL_TYPE == expect_type

    @pytest.mark.parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid SEARCH_CULL raises exceptions."""
        self.check_raises(interpreter, val, exc)


class _TestWoodsOrMatrixParam(_TestInterpretBase):
    """Tests for interpreting a parameter in Woods or matrix notation."""

    def test_interpret_matrix(self, interpreter):
        """Check successful interpretation of a matrix."""
        self.interpret(interpreter, '2 0, 0 2', flags_str='M')
        value = getattr(interpreter.rpars, self.param)
        assert value == pytest.approx(np.array([[2, 0], [0, 2]]))

    def test_interpret_woods(self, interpreter, ag100):
        """Check successful interpretation of a Woods notation."""
        interpreter.slab, *_ = ag100
        self.interpret(interpreter, 'p(2x1)')
        value = getattr(interpreter.rpars, self.param)
        assert value == pytest.approx(np.array([[2, 0], [0, 1]]))

    invalid = {
        'wood no slab': ('c(2x2)', '', err.ParameterError),
        'matrix float': ('a 1, 2 3', 'M', err.ParameterFloatConversionError),
        'matrix rows': ('1 2', 'M', err.ParameterParseError),
        'matrix cols': ('1 2, 3', 'M', err.ParameterParseError),
        }

    @pytest.mark.parametrize('val,flag,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, flag, exc, interpreter):
        """Ensure that a Woods SUPERLATTICE without a slab raises."""
        self.check_raises(interpreter, val, err.ParameterError,
                          flags_str=flag)


class TestSuperlattice(_TestWoodsOrMatrixParam):
    """Tests for interpreting SUPERLATTICE."""

    param = 'SUPERLATTICE'

    def test_interpret_matrix(self, interpreter):
        """Check successful interpretation of a SUPERLATTICE matrix."""
        super().test_interpret_matrix(interpreter)
        rpars = interpreter.rpars
        assert rpars.superlattice_defined

    def test_interpret_woods(self, interpreter, ag100):
        """Check successful interpretation of a Woods SUPERLATTICE."""
        super().test_interpret_woods(interpreter, ag100)
        rpars = interpreter.rpars
        assert rpars.superlattice_defined


class TestSymCellTransform(_TestWoodsOrMatrixParam):
    """Tests for interpreting SYMMETRY_CELL_TRANSFORM."""

    param = 'SYMMETRY_CELL_TRANSFORM'


class TestSymmetryBulk(_TestInterpretBase):
    """Tests for interpreting SYMMETRY_BULK."""

    param = 'SYMMETRY_BULK'
    valid = {
        'mirror': ('p2 m[0 1]',
                   {'mirror': {(0, 1)}, 'rotation': set(), 'group': 'p2'}),
        'cm_rot4': ('cm[1 1] r4',
                    {'mirror': set(), 'rotation': {4}, 'group': 'cm[1 1]'}),
        'pmm': ('pmm',
                {'mirror': set(), 'rotation': set(), 'group': 'pmm'}),
        }
    invalid = {'invalid_syntax': 'invalid_syntax',
               'multiple_groups': 'pmm pg',
               'missing_group': ''}

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid SYMMETRY_BULK."""
        self.check_assigned(interpreter, val, expect)

    @pytest.mark.parametrize('val', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, interpreter):
        """Ensure invalid SYMMETRY_BULK raises exceptions."""
        self.check_raises(interpreter, val, err.ParameterValueError)


class TestSymmetryEps(_TestInterpretBase):
    """Tests for interpreting SYMMETRY_EPS/SYMMETRY_EPS_Z."""

    param = 'SYMMETRY_EPS'

    def test_interpret_single_value(self, interpreter):
        """Check correct interpretation of EPS value."""
        self.check_assigned(interpreter, '0.5', 0.5)
        rpars = interpreter.rpars
        assert rpars.SYMMETRY_EPS == rpars.SYMMETRY_EPS_Z

    def test_interpret_multiple_values(self, interpreter):
        """Check correct interpretation of EPS and EPS_Z values."""
        self.interpret(interpreter, '0.1 0.2')
        rpars = interpreter.rpars
        assert (rpars.SYMMETRY_EPS, rpars.SYMMETRY_EPS_Z) == (0.1, 0.2)

    def test_interpret_invalid_number_of_inputs(self, interpreter):
        """Ensure more than two values raises exceptions."""
        self.check_raises(interpreter, '0.1 0.2 0.3',
                           err.ParameterNumberOfInputsError)

    def test_large_values_log(self, interpreter, caplog, re_match):
        """Check correct interpretation of EPS and EPS_Z values."""
        self.interpret(interpreter, '1.5 1.2')
        assert re_match(r'.*[\s\S]*SYMMETRY_EPS.*[\s\S]*very loose constraint',
                        caplog.text)


class TestSymmetryFix(_TestInterpretBase):
    """Tests for interpreting SYMMETRY_FIX."""

    param = 'SYMMETRY_FIX'
    _default = Rparams().get_default(param)

    valid = {'auto': ('t', _default),
             'p1': ('p1', 'p1'),
             'direction': ('cm[1 1]', 'cm[1 1]'),}
    invalid = {'group': ('invalid', err.ParameterParseError),
               'direction_missing': ('cm', err.ParameterParseError),
               'direction_wrong': ('pmt [0 x]', err.ParameterError),}

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid SYMMETRY_FIX."""
        self.check_assigned(interpreter, val, expect)

    @pytest.mark.parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid SYMMETRY_FIX raises exceptions."""
        self.check_raises(interpreter, val, exc)

    def test_default_restored(self, interpreter):
        """Ensure that the default value is used when 't' is given."""
        self.check_assigned(interpreter, 'f', 'p1')
        self.check_assigned(interpreter, 't', self._default)


class TestTensorOutput(_TestInterpretBase):
    """Tests for interpreting TENSOR_OUTPUT."""

    param = 'TENSOR_OUTPUT'
    valid = {'single': ('False', [0]),
             'multiple': ('0 1 1 0', [0, 1, 1, 0]),
             'repeated': ('2*1 2*0 1', [1, 1, 0, 0, 1]),}
    invalid = ('2', '5*5 0 0 1')

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid TENSOR_OUTPUT."""
        self.check_assigned(interpreter, val, expect)

    @pytest.mark.parametrize('val', invalid, ids=invalid)
    def test_interpret_invalid(self, val, interpreter):
        """Ensure invalid TENSOR_OUTPUT raises exceptions."""
        self.check_raises(interpreter, val, err.ParameterParseError)


class TestTheoEnergies(_TestInterpretBase):
    """Tests for interpreting THEO_ENERGIES."""

    param = 'THEO_ENERGIES'
    _defaults = Rparams().get_default(param)
    valid = {
        'single_value_default': ('_', _defaults),
        'single_value': ('1.0', [1.0, 1.0, 1.0]),
        'three_values_all_default': ('_ _ _', _defaults),
        'three_values_one_default': ('1.0 _ 2.0', [1.0, _defaults[1], 2.0]),
        'three_positive': ('1.0 2.0 0.5', [1.0, 2.0, 0.5]),
        'two_defaults': ('1.0 _ _', [1.0, *_defaults[1:]]),
        'range_correction': ('1.1 2.5 0.5', [1.0, 2.5, 0.5]),
        }
    invalid = {
        'one_negative': ('1.0 -2.0 0.5', err.ParameterRangeError),
        'invalid_range': ('2.0 1.0 0.5', err.ParameterValueError),
        'start_out_of_range': ('-0.5 2.0 0.5', err.ParameterRangeError),
        }

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid THEO_ENERGIES."""
        self.check_assigned(interpreter, val, expect)

    @pytest.mark.parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid THEO_ENERGIES raises exceptions."""
        self.check_raises(interpreter, val, exc)


class TestV0Real(_TestInterpretBase):
    """Tests for interpreting V0_REAL."""

    param = 'V0_REAL'
    valid = {
        'rundgren': ('rundgren 1.0 2.0 3.0 4.0', [1.0, 2.0, 3.0, 4.0]),
        'custom': ('EE', 'EEV+workfn'),
        }
    invalid = {
        'too_few': 'rundgren 1.0 2.0',
        'float': 'rundgren 1.0 2.0 3.0 four',
        }

    @pytest.mark.parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid V0_REAL."""
        self.check_assigned(interpreter, val, expect)

    @pytest.mark.parametrize('val', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, interpreter):
        """Ensure invalid V0_REAL raises exceptions."""
        self.check_raises(interpreter, val, err.ParameterError)
