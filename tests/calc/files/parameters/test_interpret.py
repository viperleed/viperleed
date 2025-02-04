"""Tests for module interpret of viperleed.calc.files.parameters."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-06-09'
__license__ = 'GPLv3+'

import logging

import numpy as np
import pytest
from pytest_cases import parametrize, parametrize_with_cases

from viperleed import __version__
from viperleed.calc.classes.rparams import Rparams
from viperleed.calc.classes.rparams.special.l_max import LMax
from viperleed.calc.classes.rparams.special.layer_cuts import (
    LayerCutToken as Cut,
    LayerCutTokenType as CutType
    )
from viperleed.calc.classes.rparams.special.search_cull import SearchCull
from viperleed.calc.files import parameters
from viperleed.calc.files.parameters import errors as err
from viperleed.calc.files.parameters.checker import ParametersChecker
from viperleed.calc.files.parameters.known_parameters import is_deprecated
from viperleed.calc.files.parameters.utils import Assignment
from viperleed.calc.files.parameters.utils import NumericBounds as Bounds
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.version import Version

from ...poscar_slabs import CasePOSCARSlabs
from .case_parameters import case_parameters_slab


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


# TODO: would be nice to reduce a bit the verbosity to use (maybe)
# a metaclass in order to implement test_interpret_(in)valid once,
# and parametrize it with the valid/invalid of each class below.
# Would it be cleaner with CASES?


@pytest.fixture(name='interpreter')
def fixture_interpreter():
    """Return a fresh ParameterInterpreter."""
    return parameters.ParameterInterpreter(Rparams())


@pytest.fixture(name='ag100_interpreter')
def fixture_ag100_interpreter(ag100, interpreter):
    """Return a ParameterInterpreter for a Ag(100) slab."""
    interpreter.slab, *_ = ag100
    return interpreter


class TestInterpreterBasics:
    """Tests for ParameterInterpreter basic functionality."""

    def test_unknown_param(self, interpreter):
        """Check complaints when an unknown parameter is requested."""
        # pylint: disable=protected-access
        with pytest.raises(err.ParameterNotRecognizedError):
            interpreter._interpret_param('UNKNOWN_PARAMETER', None)

    wrong_alias = {
        'not bool': {'abcd': ('abcd alias',)},
        'overlapping': {True: ('alias',), False: ('alias',)},
        'overlapping when joined': {True: ('false',)},
        }

    @parametrize('aliases', wrong_alias.values(), ids=wrong_alias)
    def test_bool_param_wrong_synonyms(self, aliases, interpreter):
        """Check complaints when an invalid synonym dict is passed."""
        with pytest.raises(ValueError):
            interpreter.interpret_bool_parameter(None, aliases)


class _TestInterpretBase:
    """Base class for parameter-interpretation tests."""
    param = None
    rpars_attr = None

    def assignment(self, value_str, **kwargs):
        """Return an Assignment object for self.param."""
        return Assignment(value_str, self.param, **kwargs)

    def interpret(self, interpreter, value_str, **kwargs):
        """Interpret a value for self.param."""
        if not self.param:
            raise AttributeError('Did you forget to set the '
                                 'class attribute param?')
        method = getattr(interpreter, f'interpret_{self.param.lower()}')
        return method(self.assignment(value_str, **kwargs))

    def check_assigned(self, interpreter, value_str, expected, **kwargs):
        """Assert interpretation is successful."""
        self.interpret(interpreter, value_str, **kwargs)
        attr = self.rpars_value(interpreter)
        assert expected == pytest.approx(attr)

    def check_raises(self, interpreter, value_str, exc, **kwargs):
        """Assert that an attempt to interpret raises an exc."""
        with pytest.raises(exc):
            self.interpret(interpreter, value_str, **kwargs)

    def rpars_value(self, interpreter):
        """Return the value of the attribute in interpreter.rpars."""
        return getattr(interpreter.rpars, self.rpars_attr or self.param)


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

    @parametrize('val,bounds,expect', valid.values(), ids=valid)
    # 6/5 seems OK here, considering that two are fixtures.
    # pylint: disable-next=too-many-arguments
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

    @parametrize('val,bounds', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, bounds, interpreter):
        """Ensure exceptions are raised for invalid combinations."""
        with pytest.raises(err.ParameterRangeError):
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
        with pytest.raises(err.ParameterRangeError):
            interpreter.interpret_n_bulk_layers(assignment)

    def test_interpret_t_debye(self, interpreter):
        """Check assignment of valid T_DEBYE."""
        assignment = Assignment('300.0', 'T_DEBYE')
        interpreter.interpret_t_debye(assignment)
        assert interpreter.rpars.T_DEBYE == pytest.approx(300.0)

    def test_interpret_layer_stack_vertical(self, interpreter):
        """Check assignment of valid LAYER_STACK_VERTICAL."""
        assignment = Assignment('c', 'LAYER_STACK_VERTICAL')
        interpreter.interpret_layer_stack_vertical(assignment)
        assert not interpreter.rpars.LAYER_STACK_VERTICAL


class TestAverageBeams(_TestInterpretBase):
    """Tests for interpreting AVERAGE_BEAMS."""

    param = 'AVERAGE_BEAMS'
    valid = {'off': ('off', False),
             'all': ('all', (0.0, 0.0)),
             'custom': ('45.0 60', (45.0, 60.0)),}
    invalid = {
        'invalid input': ('invalid input', err.ParameterFloatConversionError),
        }

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid AVERAGE_BEAMS."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid AVERAGE_BEAMS raises exceptions."""
        self.check_raises(interpreter, val, exc)


class TestBeamIncidence(_TestInterpretBase):
    """Tests for interpreting BEAM_INCIDENCE."""

    param = 'BEAM_INCIDENCE'
    valid = {
        'custom': ('45.0 60', (45.0, 60.0)),
        'negative theta': ('-30 30', (30.0, 210.0)),
        'comma': ('THETA 19.3, PHI 138', (19.3, 138)),
        'comma swap': ('phi 37, THETA 21', (21.0, 37.0)),
        }
    invalid = {
        'not float': ('invalid input', '', err.ParameterFloatConversionError),
        'flag': ('10 20', 'flag', err.ParameterUnknownFlagError),
        'too many': ('10 20 30', '', err.ParameterNumberOfInputsError),
        'too few': ('10 ', '', err.ParameterNumberOfInputsError),
        'invalid angle': ('invalid 10, THETA 20', '',
                          err.ParameterNumberOfInputsError),
        'angle typo': ('PHIinvalid 10, THETA 20', '',
                       err.ParameterUnknownFlagError),
        'missing angle': ('THETA 10, ', '', err.ParameterNumberOfInputsError),
        'missing angle value': ('THETA 10, PHI', '',
                                err.ParameterNumberOfInputsError),
        'missing values': ('THETA, PHI', '', err.ParameterNumberOfInputsError),
        'too many angle values': ('THETA 10 20, PHI 0', '',
                                  err.ParameterNumberOfInputsError),
        'repeated': ('PHI 130, PHI 28, THETA 30', '',
                     err.ParameterNumberOfInputsError),
        'non float angle': ('PHI abcd, THETA 30', '',
                     err.ParameterFloatConversionError),
        }

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check interpretation of valid BEAM_INCIDENCE into THETA and PHI."""
        self.interpret(interpreter, val)
        rpars = interpreter.rpars
        assert (rpars.THETA, rpars.PHI) == expect

    @parametrize('val,flag,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, flag, exc, interpreter):
        """Ensure invalid BEAM_INCIDENCE raises exceptions."""
        self.check_raises(interpreter, val, exc, flags_str=flag)


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
    invalid = {
        'not a float': ('text', err.ParameterFloatConversionError),
        'invalid direction spec': ('y(1.2)', err.ParameterParseError),
        'z, not float': ('z(abc)', err.ParameterParseError),
        'z, invalid float': ('z(1.3.5)', err.ParameterFloatConversionError),
        'vector not enough items': ('[]', err.ParameterParseError),
        }

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, ag100_interpreter):
        """Check correct interpretation of valid BULK_REPEAT."""
        self.interpret(ag100_interpreter, val)
        rpars = ag100_interpreter.rpars
        assert rpars.BULK_REPEAT == pytest.approx(expect, rel=1e-4)

    @parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, ag100_interpreter):
        """Ensure invalid BULK_REPEAT raises exceptions."""
        self.check_raises(ag100_interpreter, val, exc)

    def test_slab_missing(self, interpreter):
        """Check complaints if a slab is not passed."""
        self.check_raises(interpreter, '1.5', err.ParameterNeedsSlabError)


class TestDomain(_TestInterpretBase):
    """Tests for interpreting DOMAIN."""

    param = 'DOMAIN'

    def test_interpret_path_with_flag(self, interpreter, tmp_path):
        """Test correct interpretation of a path with a domain name."""
        domain_path = tmp_path / 'domain1'
        domain_path.mkdir()
        self.interpret(interpreter, str(domain_path), flags_str='domain1')
        assert interpreter.rpars.DOMAINS == {'domain1': domain_path}

    def test_interpret_path_no_flag(self, interpreter, tmp_path):
        """Test correct interpretation of a path without a domain name."""
        self.interpret(interpreter, str(tmp_path))
        assert interpreter.rpars.DOMAINS == {'1': tmp_path}

    def test_interpret_path_relative_to_cwd(self, interpreter, tmp_path):
        """Test correct interpretation of a path relative to cwd."""
        relative_path = 'domain'
        domain_path = tmp_path / relative_path
        domain_path.mkdir()
        with execute_in_dir(tmp_path):
            self.interpret(interpreter, relative_path)
        assert interpreter.rpars.DOMAINS == {'1': domain_path}

    def test_interpret_path_relative_to_calc(self, interpreter, tmp_path):
        """Test interpretation of a path relative to where calc was started."""
        relative_path = 'domain'
        calc_path = tmp_path / 'calc_was_started_here'
        domain_path = calc_path / relative_path
        domain_path.mkdir(parents=True)
        interpreter.rpars.paths.home = calc_path
        self.interpret(interpreter, relative_path)
        assert interpreter.rpars.DOMAINS == {'1': domain_path}

    def test_interpret_zip_file(self, interpreter, tmp_path):
        """Test correct interpretation of a zip file."""
        zip_file = tmp_path / 'domain.zip'
        zip_file.touch()
        self.interpret(interpreter, str(zip_file))
        assert interpreter.rpars.DOMAINS == {'1': zip_file}
        self.interpret(interpreter, str(zip_file.with_suffix('')))
        assert interpreter.rpars.DOMAINS == {'1': zip_file,
                                             '2': zip_file}

    def test_interpret_invalid(self, interpreter):
        """Ensure invalid DOMAIN raises exceptions."""
        self.check_raises(interpreter, 'invalid_path', err.ParameterValueError)

    def test_duplicate_name(self, interpreter):
        """Ensure that no two domains can have the same name."""
        interpreter.rpars.DOMAINS['domain'] = 'path_to_domain'
        self.check_raises(interpreter, 'path_to_domain',
                          err.ParameterValueError, flags_str='domain')


class TestDomainStep(_TestInterpretBase):
    """Tests for interpreting DOMAIN_STEP."""

    param = 'DOMAIN_STEP'
    valid = {'value': ('10', 10),}
    invalid = {
        'value': ('200', err.ParameterRangeError),
        'out of range low': ('-5', err.ParameterRangeError),
        'non integer': ('5.5', err.ParameterIntConversionError),
        'does not divide 100': ('38', err.ParameterValueError),
        'too many': ('12 13', err.ParameterNumberOfInputsError),
        }

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid DOMAIN_STEP."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid DOMAIN_STEP raises exceptions."""
        self.check_raises(interpreter, val, exc)


class _TestSlabNotEmpty(_TestInterpretBase):
    """Collection of tests for parameters that need a non-empty slab."""

    def test_missing_slab(self, interpreter):
        """Check complaints when no slab is provided."""
        self.check_raises(interpreter,
                          'value_that_wont_be_interpreted',
                          err.ParameterNeedsSlabError,
                          flags_str='Ag')

    def test_empty_slab(self, ag100_interpreter):
        """Check complaints when an empty slab is provided."""
        ag100_interpreter.slab.atlist.clear()
        with pytest.raises(err.ParameterError) as exc:
            self.interpret(ag100_interpreter,
                           'value_that_wont_be_interpreted',
                           flags_str='Ag')
        assert exc.match('.*no atoms')

    def test_no_elements(self, ag100_interpreter):
        """Check complaints when an empty slab is provided."""
        ag100_interpreter.slab.n_per_elem.clear()
        with pytest.raises(err.ParameterError) as exc:
            self.interpret(ag100_interpreter,
                           'value_that_wont_be_interpreted',
                           flags_str='Ag')
        assert exc.match('.*no elements')


class TestElementMix(_TestSlabNotEmpty):
    """Tests for interpreting ELEMENT_MIX."""

    param = 'ELEMENT_MIX'
    valid = {
        'two elements': ('Ag', 'Co Cu', {'Ag': ['Co', 'Cu']}),
        }
    invalid = {
        'too many flags': ('D E', 'Mn Fe Co', err.ParameterUnknownFlagError),
        'not chem elem': ('Ag', 'chem_el_unknown Fe', err.ParameterValueError),
        'no elements': ('Ag', '   ', err.ParameterHasNoValueError),
        'not-poscar elem': ('Mn', 'Co Cu', err.ParameterUnknownFlagError),
        'unknown element': ('X', 'Au Hg', err.ParameterUnknownFlagError),
        'single element': ('Ag', 'Mn', err.ParameterNumberOfInputsError),
        }

    @parametrize('poscar_el,mix,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, poscar_el, mix, expect, ag100_interpreter):
        """Check correct interpretation of valid FILAMENT_WF."""
        self.check_assigned(ag100_interpreter, mix, expect, flags_str=poscar_el)

    @parametrize('poscar_el,mix,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, poscar_el, mix, exc, ag100_interpreter):
        """Ensure invalid FILAMENT_WF raises exceptions."""
        self.check_raises(ag100_interpreter, mix, exc, flags_str=poscar_el)


class TestElementRename(_TestSlabNotEmpty):
    """Tests for interpreting ELEMENT_RENAME."""

    param = 'ELEMENT_RENAME'

    def test_interpret_valid(self, ag100_interpreter):
        """Check correct interpretation of valid ELEMENT_RENAME."""
        self.interpret(ag100_interpreter, 'H', flags_str='Ag')
        assert ag100_interpreter.rpars.ELEMENT_RENAME == {'Ag': 'H'}

    invalid = {
        'not chem elem': ('Ag', 'Op', err.ParameterValueError),
        'no elements': ('Ag', '', err.ParameterHasNoValueError),
        'not-poscar elem': ('Mn', 'Co', err.ParameterUnknownFlagError),
        'too many flags': ('Ag B', 'Mn Fe', err.ParameterUnknownFlagError),
        'no flag': ('', 'C', err.ParameterNeedsFlagError),
        }

    @parametrize('flag,val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, flag, val, exc, ag100_interpreter):
        """Ensure that an invalid chemical element raises exceptions."""
        self.check_raises(ag100_interpreter, val, exc, flags_str=flag)


class TestFilamentWF(_TestInterpretBase):
    """Tests for interpreting FILAMENT_WF."""

    param = 'FILAMENT_WF'
    valid = {'lab6': ('LaB6', 2.65),
             'custom': ('1.0', 1.0),}
    invalid = {
        'invalid_float': ('invalid', '', err.ParameterFloatConversionError),
        'flag': ('1.5', 'test', err.ParameterUnknownFlagError),
        }

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid FILAMENT_WF."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val,flag,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, flag, exc, interpreter):
        """Ensure invalid FILAMENT_WF raises exceptions."""
        self.check_raises(interpreter, val, exc, flags_str=flag)


class TestFortranComp(_TestInterpretBase):
    """Tests for interpreting FORTRAN_COMP."""

    param = 'FORTRAN_COMP'
    valid = {  # id: (assign_flag, assign_value, expect_pre, expect_post)
        'default_intel': ('', 'ifort',
                          'ifort -O2 -I/opt/intel/mkl/include',
                          '-L/opt/intel/mkl/lib/intel64'),
        'default_gnu': ('', 'gfortran', 'gfortran',
                        '-llapack -lpthread -lblas'),
        'default_intel_mpi': ('mpi', 'mpiifort', 'mpiifort', None),
        'default_gnu_mpi': ('mpi', 'mpifort', 'mpifort -Ofast -no-pie', None),
        'custom_no_flag': ('', '"ifort -O3 -march=native"',
                           'ifort -O3 -march=native', None),
        'custom_post_flag': ('post', '"-L/opt/intel/mkl/lib/intel64"',
                             None, '-L/opt/intel/mkl/lib/intel64'),
        'custom_mpi_flag': ('mpi', '"mpifort -fallow-argument-mismatch"',
                            'mpifort -fallow-argument-mismatch', None),
        'custom_mpi_post_flag': ('mpipost', '"-L/opt/intel/mkl/lib/intel64"',
                                 None, '-L/opt/intel/mkl/lib/intel64'),
        'without quotes': ('mpipost', '-L/opt/intel/mkl/lib/intel64',
                           None, '-L/opt/intel/mkl/lib/intel64'),
        }
    invalid = {
        'too many flags': ('f1 f2', 'gfortran', err.ParameterUnknownFlagError),
        'unknown flag': ('invalid', '"ifort"', err.ParameterUnknownFlagError),
        'single quote': ('post', '"-L/opt/intel/mkl/lib/intel64',
                         err.ParameterValueError)
        }

    # About the disable below: In principle 'pre' and 'post' could
    # be merged into a tuple, but the parametrization above would
    # look even more complex
    @parametrize('flag,val,pre,post', valid.values(), ids=valid)
    # pylint: disable-next=too-many-arguments
    def test_interpret_valid(self, flag, val, pre, post,
                             interpreter, subtests):
        """Check correct interpretation of valid FORTRAN_COMP(_MPI)."""
        assignment = self.assignment(val, flags_str=flag)
        rpars = interpreter.rpars
        interpreter.interpret_fortran_comp(assignment)
        finish_interpret = ParametersChecker()
        # pylint: disable=protected-access
        finish_interpret._rpars = rpars
        try:
            finish_interpret._check_and_update_fortran_comp()
        except FileNotFoundError:
            pytest.skip(f'Compiler {val} not available')
        # pylint: enable=protected-access
        attr_name = self.param + ('_MPI' if 'mpi' in flag else '')
        compiler = getattr(rpars,  attr_name)
        if pre is not None:
            with subtests.test('Check pre'):
                assert pre in compiler[0]
        if post is not None:
            with subtests.test('Check post'):
                assert post in compiler[1]

    @parametrize('flag,val,exc', invalid.values(), ids=invalid)
    def test_invalid(self, flag, val, exc, interpreter):
        """Ensure invalid FORTRAN_COMP raises exceptions."""
        self.check_raises(interpreter, val, exc, flags_str=flag)


class TestIntpolDeg(_TestInterpretBase):
    """Tests for interpreting INTPOL_DEG."""

    param = 'INTPOL_DEG'
    valid = {v: (v, int(v)) for v in Rparams().get_limits(param)}
    invalid = '1', 'text'

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid INTPOL_DEG."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val', invalid, ids=invalid)
    def test_interpret_invalid(self, val, interpreter):
        """Ensure invalid INTPOL_DEG raises exceptions."""
        self.check_raises(interpreter, val, err.ParameterValueError)


class TestIVShiftRange(_TestInterpretBase):
    """Tests for interpreting IV_SHIFT_RANGE."""

    param = 'IV_SHIFT_RANGE'
    _defaults = Rparams.get_default('IV_SHIFT_RANGE')
    valid = {
        'range': ('0.0 1.0 0.25', [0.0, 1.0, 0.25]),
        'default bound': ('_ 1.0 0.25', [_defaults.start, 1.0, 0.25]),
        'default step': ('-2.5 1.0 _', [-2.5, 1.0, _defaults.step]),
        'no step': ('-2.5 1.0', [-2.5, 1.0, _defaults.step]),
        'swapped': ('3 -2 -0.5', [-2.0, 3.0, 0.5]),
        'adjusted': ('1.05 2.05 0.2', [1.0, 2.2, 0.2]),
        }
    invalid = {
        'nr_inputs': ('0.0', err.ParameterNumberOfInputsError),
        'float': ('0.0 2.0 ()', err.ParameterFloatConversionError),
        'parse': ('0.0 2.0 1.a', err.ParameterParseError),
        'swapped opposite step': ('1.0 0.0 0.1', err.ParameterValueError),
        'opposite step': ('0.0 1.0 -0.1', err.ParameterValueError),
        }

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid IV_SHIFT_RANGE."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid IV_SHIFT_RANGE raises exceptions."""
        self.check_raises(interpreter, val, exc)

    def test_adjusted_logs(self, interpreter, caplog):
        """Check logging messages when IV_SHIFT_RANGE bounds are modified."""
        val, expect = self.valid['adjusted']
        with caplog.at_level(0):
            self.check_assigned(interpreter, val, expect)
        assert 'integer multiple' in caplog.text


class TestLayerCuts(_TestInterpretBase):
    """Tests for interpreting LAYER_CUTS."""

    param = 'LAYER_CUTS'
    valid = {
        'simple': ('0.1 0.2 0.3', ['0.1', '0.2', '0.3']),
        'dz': ('dz(1.2)', ['dz(1.2)']),
        'list and dz': (
            '0.15 0.3 < dz(1.2)  ',
            ['0.15', '0.3', Cut(CutType.AUTO_DZ, 1.2, 0.3, None)]
            ),
        'two dz': (
            'dz(1.0) < 0.28 < dz(0.5) < 0.6',
            [Cut(CutType.AUTO_DZ, 1.0, None, 0.28), '0.28',
             Cut(CutType.AUTO_DZ, 0.5, 0.28, 0.60), '0.60']
            ),
        'larger than': (
            '0.75 0.3 > dz(1.2)  ',
            [Cut(CutType.AUTO_DZ, 1.2, None, 0.3), '0.3', '0.75']
            ),
        }
    invalid = {
        'less and greater': ('< 0.1 > 0.2', err.ParameterParseError),
        'cutoff function': ('dz(0.1) dc(0.2) invalid(0.3)',
                            err.ParameterParseError),
        'float': ('0.1 invalid 0.3', err.ParameterParseError),
        'dz': ('0.5 1.0 < dz(abcd) < 4.0', err.ParameterParseError),
        'no cut': ('', err.ParameterValueError),
        }

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid LAYER_CUTS."""
        self.interpret(interpreter, val)
        value = self.rpars_value(interpreter)
        assert all(cut == cut_expected
                   for cut, cut_expected in zip(value, expect))

    @parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid LAYER_CUTS raises exceptions."""
        self.check_raises(interpreter, val, exc)


class TestLMax(_TestInterpretBase):
    """Tests for interpreting LMAX."""

    param = 'LMAX'
    valid = {
        'single': ('7', LMax(7)),
        'single with spaces': ('6    ', LMax(6)),
        'range hyphen': ('9-12', LMax(9, 12)),
        'range colon': ('8:18', LMax(8, 18)),
        'range space': ('6 17', LMax(6, 17)),
        'range swapped': ('15-11', LMax(11, 15)),
        'space and hyphen': ('12 - 15', LMax(12, 15)),
        'multi space': ('7     10', LMax(7, 10)),
        'multi space and hyphen': ('9    -         12', LMax(9, 12)),
        'stripped': ('  12 - 9     ', LMax(9, 12)),
        'two hyphens': ('5--16', LMax(5, 16)),
        'three hyphens': ('1---2', LMax(1, 2)),
        'four hyphens': ('5----17', LMax(5, 17)),
        'three hyphens and spaces': ('1 --- 2', LMax(1, 2)),
        }
    invalid = {
        'no value': ('   ', err.ParameterHasNoValueError),
        'too many': ('9 12 18', err.ParameterNumberOfInputsError),
        'not int': ('1.3', err.ParameterValueError),
        'not numeric': ('a', err.ParameterValueError),
        'wrong delimiter': ('1 & 2', err.ParameterValueError),
        'out of range single': ('0', err.ParameterRangeError),
        'out of range min': ('0 16', err.ParameterRangeError),
        'out of range max': ('6 1600', err.ParameterRangeError),
        'negative': ('-4', err.ParameterRangeError),
        'range both negative': ('-6 -8', err.ParameterRangeError),
        'range min negative swapped': ('1 -5', err.ParameterRangeError),
        'range min negative hyphen': ('-2 -     5', err.ParameterRangeError),
        'range max negative hyphens': ('1-- -5', err.ParameterRangeError),
        'ends with delimiter': ('2-', err.ParameterValueError),
        'ends with delimiter strip': ('2 -  ', err.ParameterValueError),
        }

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid LMAX."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid LMAX raises exceptions."""
        self.check_raises(interpreter, val, exc)


class TestLogLevel(_TestInterpretBase):
    """Tests for interpreting LOG_LEVEL."""

    param = 'LOG_LEVEL'

    def test_interpret_true(self, interpreter):
        """Check correct setting of LOG_LEVEL to True."""
        self.interpret(interpreter, 'true')
        assert interpreter.rpars.LOG_LEVEL <= logging.DEBUG

    def test_interpret_false(self, interpreter):
        """Check correct setting of LOG_LEVEL to False."""
        self.interpret(interpreter, 'F')
        assert interpreter.rpars.LOG_LEVEL >= logging.INFO

    def test_interpret_int(self, interpreter):
        """Check correct interpretation of integer LOG_LEVEL."""
        self.check_assigned(interpreter, '3', 3)

    def test_interpret_default(self, interpreter):
        """Check correct interpretation of a default string LOG_LEVEL."""
        self.interpret(interpreter, 'vv')
        assert interpreter.rpars.LOG_LEVEL < logging.DEBUG / 2

    invalid = {
        'too many': ('not a level', '', err.ParameterNumberOfInputsError),
        'too few': ('', '', err.ParameterHasNoValueError),
        'flag': ('true', 'flag', err.ParameterUnknownFlagError),
        'invalid str': ('not_a_level', '', err.ParameterValueError),
        }

    @parametrize('val,flag,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, flag, exc, interpreter):
        """Check correct interpretation of a default string LOG_LEVEL."""
        self.check_raises(interpreter, val, exc, flags_str=flag)


class TestOptimize(_TestInterpretBase):
    """Tests for interpreting OPTIMIZE."""

    param = 'OPTIMIZE'
    valid = {  # value, flag, {key: expected}
        'flag_and_value': ('step 0.1', 'v0i', {'step': 0.1}),
        'no flag at right': ('0.8', 'v0i', {'step': 0.8}),
        'multiple_flag_value': (
            'step 0.1, convergence 1e-6, minpoints 10', 'theta',
            {'step': 0.1, 'minpoints': 10, 'convergence': 1e-6}
            ),
        }
    invalid = {
        'missing_quantity': ('step 0.1', '', err.ParameterNeedsFlagError),
        'flag': ('invalid 0.1', 'v0i', err.ParameterUnknownFlagError),
        'value': ('step not-a-number', 'v0i', err.ParameterValueError),
        'quantity': ('step 0.1', 'invalid', err.ParameterUnknownFlagError),
        'non-float step': ('abcd', 'phi', err.ParameterFloatConversionError),
        'empty value': ('step 0.1, minpoints', 'phi',
                        err.ParameterNumberOfInputsError),
        'empty flag-value pair': ('step 0.1, , minpoints 5', 'phi',
                                  err.ParameterNumberOfInputsError),
        }

    @parametrize('val,flag,expect', valid.values(), ids=valid)
    # 6/5 Seems OK here, especially considering that two are fixtures
    # pylint: disable-next=too-many-arguments
    def test_interpret_valid(self, val, flag, expect, interpreter, subtests):
        """Check correct interpretation of valid OPTIMIZE."""
        rpars = interpreter.rpars
        rpars.RUN.append(6)
        self.interpret(interpreter, val, flags_str=flag)
        with subtests.test('which'):
            assert rpars.OPTIMIZE['which'] == flag
        for key, value in expect.items():
            with subtests.test(key):
                assert rpars.OPTIMIZE[key] == pytest.approx(value)

    @parametrize('val,flag,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, flag, exc, interpreter):
        """Ensure invalid OPTIMIZE raises exceptions."""
        rpars = interpreter.rpars
        rpars.RUN.append(6)
        self.check_raises(interpreter, val, exc, flags_str=flag)

    _dont_interpret = {id_: (value, flag)  # No exception info
                       for id_, (value, flag, _) in invalid.items()}
    _dont_interpret['flag_and_value'] = valid['flag_and_value'][:2]

    @parametrize('value,flag', _dont_interpret.values(), ids=_dont_interpret)
    def test_interpret_superfluous(self, value, flag,
                                   interpreter, caplog, subtests):
        """Ensure that defining OPTIMIZE without RUN=6 complains."""
        self.interpret(interpreter, value, flags_str=flag)
        rpars = interpreter.rpars
        interpreted = self.rpars_value(interpreter)
        with subtests.test('value unchanged'):
            assert interpreted == rpars.get_default(self.param)
        with subtests.test('log warning'):
            assert 'RUN does not include' in caplog.text


class TestParabolaFit(_TestInterpretBase):
    """Tests for interpreting PARABOLA_FIT."""

    param = 'PARABOLA_FIT'
    _defaults = getattr(Rparams(), param)
    valid = {
        'off': ('off', {**_defaults, 'type': 'none'}),
        'localize': ('type linear, localise 0.3',
                     {**_defaults, 'type': 'linear', 'localize': 0.3}),
        'none': ('type none, alpha 3.17',
                 {**_defaults, 'type': 'none', 'alpha': 3.17}),
        }
    invalid = {
        'type': ('type not_a_known_type', err.ParameterValueError),
        'float': ('type linear, alpha abc', err.ParameterFloatConversionError),
        'negative': ('type linear, alpha -1.3', err.ParameterRangeError),
        'flag': ('type linear, who_knows 0.5', err. ParameterValueError),
        'too many': ('type linear abcd', err. ParameterNumberOfInputsError),
        'off+invalid': ('off this_is_wrong, so is_this, and this too',
                        err.ParameterValueError),
        }

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid PARABOLA_FIT."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid PARABOLA_FIT raises exceptions."""
        self.check_raises(interpreter, val, exc)

    @pytest.mark.skipif(not is_deprecated(param, __version__),
                        reason='Not deprecated')
    def test_deprecated(self, interpreter, caplog):
        """Ensure that PARABOLA_FIT is flagged as deprecated."""
        val, *_ = self.valid['off']
        self.interpret(interpreter, val)
        assert 'deprecated' in caplog.text


class TestPhaseshiftEps(_TestInterpretBase):
    """Tests for interpreting PHASESHIFT_EPS."""

    param = 'PHASESHIFT_EPS'
    valid = {'float': ('0.1', 0.1),
             'tag': ('fine', 0.01),}
    invalid = {
        'float': ('invalid', '', err.ParameterFloatConversionError),
        'negative': ('-1.0', '', err.ParameterRangeError),
        'too large': ('1.5', '', err.ParameterRangeError),
        'too many': ('1.3 2.4', '', err.ParameterNumberOfInputsError),
        'flag': ('1.23', 'flag', err.ParameterUnknownFlagError),
        }

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid PHASESHIFT_EPS."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val,flag,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, flag, exc, interpreter):
        """Ensure invalid PHASESHIFT_EPS raises exceptions."""
        self.check_raises(interpreter, val, exc, flags_str=flag)


class TestPlotIV(_TestInterpretBase):
    """Tests for interpreting parameter PLOT_IV."""

    param = 'PLOT_IV'
    valid = {   # flag, value, attribute, expected
        'plot true': ('plot', 'true', 'plot', True),
        'plot true t': ('plot', 't', 'plot', True),
        'plot true one': ('plot', '1', 'plot', True),
        'plot false': ('plot', 'false', 'plot', False),
        'plot false f': ('plot', 'f', 'plot', False),
        'plot false zero': ('plot', '0', 'plot', False),
        'legend none': ('legend', 'none', 'legend', 'none'),
        'legend first': ('legends', 'first', 'legend', 'first'),
        'legend top': ('legend', 'topright', 'legend', 'tr'),
        'colors': ('colors', 'blue red', 'colors', ('blue', 'red')),
        'colors color': ('color', 'blue red', 'colors', ('blue', 'red')),
        'colors colour': ('colour', 'blue red', 'colors', ('blue', 'red')),
        'colors colours': ('colours', 'blue red', 'colors', ('blue', 'red')),
        'colors rgb hex': ('colors', '#FF0000 purple',
                           'colors', ('#FF0000', 'purple')),
        'axes': ('border', 'none', 'axes', 'none'),
        'axes all': ('borders', 'all', 'axes', 'all'),
        'axes bottom': ('axes', 'bottom', 'axes', 'b'),
        'axes less': ('axes', 'less', 'axes', 'lb'),
        'overbar true': ('overbar', 'truthfulvalue', 'overbar', True),
        'overbar false': ('overline', 'falsified', 'overbar', False),
        'perpage single': ('perpage', '3', 'perpage', 3),
        'perpage two': ('layout', '3 12', 'perpage', (3, 12)),
        }
    invalid = {
        'no flag': ('', 'abcd', err.ParameterNeedsFlagError),
        'unknown flag': ('invalid', 'value', err.ParameterUnknownFlagError),
        'invalid axes': ('borders', 'invalid', err.ParameterParseError),
        'invalid legend': ('legends', 'invalid', err.ParameterParseError),
        'invalid overbar': ('overline', 'invalid', err.ParameterParseError),
        'invalid perpage': ('perpage', 'wrong',
                            err.ParameterIntConversionError),
        'invalid perpage two': ('perpage', '5 wrong',
                                err.ParameterIntConversionError),
        'perpage neg': ('perpage', '-5',  err.ParameterRangeError),
        'perpage neg two': ('perpage', '-5 -3',  err.ParameterRangeError),
        'empty axes': ('axes', '', err.ParameterHasNoValueError),
        'empty colors': ('colors', '', err.ParameterHasNoValueError),
        'invalid rgb': ('colors', '255 0 0 green', err.ParameterValueError),
        'empty legend': ('legend', '', err.ParameterHasNoValueError),
        'empty overbar': ('overbar', '', err.ParameterHasNoValueError),
        'empty perpage': ('perpage', '', err.ParameterHasNoValueError),
        'empty plot': ('plot', '', err.ParameterHasNoValueError),
        'two axes': ('axes', 'all none', err.ParameterNumberOfInputsError),
        'two legend': ('legend', 'tr none', err.ParameterNumberOfInputsError),
        'two overbar': ('overbar', 't f', err.ParameterNumberOfInputsError),
        'two plot': ('plot', 'true 0', err.ParameterNumberOfInputsError),
        'perpage x3': ('perpage', '1 2 3', err.ParameterNumberOfInputsError),
        'multi flag': ('plot axes', 'true', err.ParameterUnknownFlagError),
        }

    @parametrize('flag,val,attr,expect', valid.values(), ids=valid)
    # About the disable: Perhaps one could pack the attr and
    # expect together, but seems a bit overkill for 6/5.
    # pylint: disable-next=too-many-arguments
    def test_interpret_valid(self, flag, val, attr, expect, interpreter):
        """Check correct interpretation of valid PLOT_IV."""
        self.interpret(interpreter, val, flags_str=flag)
        assert interpreter.rpars.PLOT_IV[attr] == expect

    @parametrize('flag,val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, flag, val, exc, interpreter):
        """Ensure invalid PLOT_IV raises exceptions."""
        self.check_raises(interpreter, val, exc, flags_str=flag)


class TestRun(_TestInterpretBase):
    """Tests for interpreting RUN."""

    param = 'RUN'
    valid = {
        'single': ('1', [0, 1]),
        'multiple': ('1 2 3', [0, 1, 2, 3]),
        'range': ('1-3', [0, 1, 2, 3]),
        }
    invalid = {'section': ('invalid', err.ParameterValueError),
               'empty': ('', err.ParameterHasNoValueError),
               'syntax underscore': ('1 _ 2', err.ParameterValueError),
               'syntax section': ('1, x', err.ParameterValueError)}

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid RUN."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val,exc', invalid.values(), ids=invalid)
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
               'more values': ('0 1', err.ParameterNumberOfInputsError)}

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid SEARCH_BEAMS."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid SEARCH_BEAMS raises exceptions."""
        self.check_raises(interpreter, val, exc)


class TestSearchConvergence(_TestInterpretBase):
    """Tests for interpreting SEARCH_CONVERGENCE."""

    param = 'SEARCH_CONVERGENCE'
    valid_gauss = {'scaling': ('0.01 0.9', (0.01, 0.9)),
                   'no scaling': ('0.01', (0.01, 0.5)),}
    valid_dgen = {'positive': ('1 1.5', (1, 1.5)),}
    invalid = {
        'gaussian no values': ('', 'gaussian', err.ParameterHasNoValueError),
        'dgen no values': ('', 'dgen', err.ParameterHasNoValueError),
        'no flag': ('0.2 0.5', '', err.ParameterNeedsFlagError),
        'unknown flag': (' 0.01 0.9', 'unknown',
                         err.ParameterUnknownFlagError),
        'too many values': ('.01 0.9 0.5', 'gaussian',
                            err.ParameterNumberOfInputsError),
        'scaling': ('0.01 -0.5', 'gaussian', err.ParameterRangeError),
        'gaussian neg': ('-3.5 0.2', 'gaussian', err.ParameterRangeError),
        'no float': ('a 0.3', 'gaussian', err.ParameterFloatConversionError),
        'dgen neg': ('-1 5', 'dgen dec', err.ParameterRangeError),
        'dgen invalid flag': ('1 1', 'dgen invalid',
                              err.ParameterUnknownFlagError),
        'dgen dec small': ('100 0.1', 'dgen dec', err.ParameterRangeError),
        'gaussian too many flags': ('0.5 0.3', 'gaussian invalid',
                                    err.ParameterUnknownFlagError),
        'dgen too many flags': ('100 3', 'dgen dec invalid',
                                err.ParameterUnknownFlagError),
        }

    @parametrize('val,expect', valid_gauss.values(), ids=valid_gauss)
    def test_interpret_valid_gaussian(self, val, expect, interpreter):
        """Check correct interpretation of SEARCH_CONVERGENCE gaussian."""
        self.interpret(interpreter, val, flags_str='gaussian')
        rpars = interpreter.rpars
        assert (rpars.GAUSSIAN_WIDTH, rpars.GAUSSIAN_WIDTH_SCALING) == expect

    @parametrize('val,expect', valid_dgen.values(), ids=valid_dgen)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of SEARCH_CONVERGENCE gaussian."""
        self.interpret(interpreter, val, flags_str='dgen')
        rpars = interpreter.rpars
        assert (rpars.SEARCH_MAX_DGEN['dec'],
                rpars.SEARCH_MAX_DGEN_SCALING['dec']) == expect

    def test_off(self, interpreter):
        """Check correct interpretation of 'off'."""
        self.interpret(interpreter, 'off')
        assert interpreter.rpars.GAUSSIAN_WIDTH_SCALING == 1.0

    @parametrize('val,flag,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, flag, exc, interpreter):
        """Ensure invalid SEARCH_CONVERGENCE raises exceptions."""
        self.check_raises(interpreter, val, exc, flags_str=flag)

    def test_no_change_in_gaussian_when_updating(self, interpreter):
        """Check that gaussian is unchanged while updating from file."""
        rpars = interpreter.rpars
        rpars.GAUSSIAN_WIDTH = rpars.searchConvInit['gaussian'] = 2.0
        assignment = self.assignment('2 0.5', flags_str='gaussian')
        interpreter.interpret_search_convergence(assignment, is_updating=True)
        assert rpars.GAUSSIAN_WIDTH == 2.0

    def test_gaussian_changed_when_updating(self, interpreter):
        """Check that gaussian and init are changed while updating."""
        rpars = interpreter.rpars
        rpars.GAUSSIAN_WIDTH = rpars.searchConvInit['gaussian'] = 2.0
        assignment = self.assignment('8.5 0.5', flags_str='gaussian')
        interpreter.interpret_search_convergence(assignment, is_updating=True)
        assert rpars.GAUSSIAN_WIDTH == 8.5
        assert rpars.searchConvInit['gaussian'] == 8.5

    def test_no_change_in_dgen_when_updating(self, interpreter):
        """Check that dgen is unchanged while updating from file."""
        rpars = interpreter.rpars
        rpars.SEARCH_MAX_DGEN['dec'] = 3
        rpars.searchConvInit['dgen']['dec'] = 3
        assignment = self.assignment('3 1.5', flags_str='dgen dec')
        interpreter.interpret_search_convergence(assignment, is_updating=True)
        assert rpars.SEARCH_MAX_DGEN['dec'] == 3

    def test_dgen_changed_when_updating(self, interpreter):
        """Check that dgen and init are changed while updating."""
        rpars = interpreter.rpars
        rpars.SEARCH_MAX_DGEN['dec'] = 3
        rpars.searchConvInit['dgen']['dec'] = 3
        assignment = self.assignment('200 1.5', flags_str='dgen dec')
        interpreter.interpret_search_convergence(assignment, is_updating=True)
        assert rpars.SEARCH_MAX_DGEN['dec'] == 200
        assert rpars.searchConvInit['dgen']['dec'] == 200


class TestSearchCull(_TestInterpretBase):
    """Tests for interpreting SEARCH_CULL."""

    param = 'SEARCH_CULL'
    _default = Rparams.get_default(param)
    valid = {'float': ('0.5', SearchCull(0.5, _default.type_)),
             'int': ('3', SearchCull(3, _default.type_)),
             'int-like float': ('4.0', SearchCull(4, _default.type_)),
             'explicit type': ('0.48 clone', SearchCull(0.48, 'clone'))}
    invalid = {
        'float greater than one': ('1.5', err.ParameterValueError),
        'negative': ('-0.5', err.ParameterValueError),
        'too many': ('0.5 clone 1.0', err.ParameterNumberOfInputsError),
        'invalid cull type': ('0.5 test', err.ParameterValueError),
        'not a float': ('abcd', err.ParameterFloatConversionError),
        }

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid SEARCH_CULL."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid SEARCH_CULL raises exceptions."""
        self.check_raises(interpreter, val, exc)


class TestSiteDefInvalid(_TestSlabNotEmpty):
    """Tests of invalid conditions when interpreting SITE_DEF."""

    param = 'SITE_DEF'

    invalid = {  # element, spec, exception
        'wrong element': ('Te', 'surf 1-10', err.ParameterUnknownFlagError),
        'top twice': ('Ag', 'topmost top(2) others top(10)',
                      err.ParameterValueError),
        'too few specs': ('Ag', 'spec 1-3, surf', err.ParameterParseError),
        'unknown specs': ('Ag', 'surf abcd', err.ParameterParseError),
        'empty spec': ('Ag', 'surf 1,, below 3',
                       err.ParameterNumberOfInputsError),
        'no flag': ('', 'surf 12', err.ParameterNeedsFlagError),
        'too many flags': ('Ag another', 'top 12',
                           err.ParameterUnknownFlagError),
        'atoms dont exist': ('Ag', 'top 1500-1520', err.ParameterValueError),
        }

    @parametrize('element,val,exc', invalid.values(), ids=invalid)
    def test_invalid(self, element, val, exc, ag100_interpreter):
        """Check complaints when an empty slab is provided."""
        self.check_raises(ag100_interpreter, val, exc, flags_str=element)

    def test_swapped_elements(self, interpreter):
        """Ensure complaints when selecting atoms of the wrong element."""
        fe3o4, *_ = CasePOSCARSlabs().case_poscar_fe3o4_001_cod()
        interpreter.slab = fe3o4
        self.check_raises(interpreter, 'its_actually_iron 1',
                          err.ParameterValueError, flags_str='O')


class _TestWoodsOrMatrixParam(_TestInterpretBase):
    """Tests for interpreting a parameter in Woods or matrix notation."""

    def test_interpret_matrix(self, interpreter):
        """Check successful interpretation of a matrix."""
        self.check_assigned(interpreter, '2 0, 0 2',
                            np.array([[2, 0], [0, 2]]), flags_str='M')

    def test_interpret_woods(self, ag100_interpreter):
        """Check successful interpretation of a Woods notation."""
        self.check_assigned(ag100_interpreter, 'p(2x1)',
                            np.array([[2, 0], [0, 1]]))

    invalid = {
        'wood no slab': ('c(2x2)', '', err.ParameterNeedsSlabError),
        'matrix float': ('a 1, 2 3', 'M', err.ParameterFloatConversionError),
        'matrix rows': ('1 2', 'M', err.ParameterParseError),
        'matrix cols': ('1 2, 3', 'M', err.ParameterParseError),
        }

    @parametrize('val,flag,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, flag, exc, interpreter):
        """Ensure that a Woods SUPERLATTICE without a slab raises."""
        self.check_raises(interpreter, val, exc, flags_str=flag)


class TestSuperlattice(_TestWoodsOrMatrixParam):
    """Tests for interpreting SUPERLATTICE."""

    param = 'SUPERLATTICE'

    def test_interpret_matrix(self, interpreter):
        """Check successful interpretation of a SUPERLATTICE matrix."""
        super().test_interpret_matrix(interpreter)
        rpars = interpreter.rpars
        assert rpars.superlattice_defined

    def test_interpret_woods(self, ag100_interpreter):
        """Check successful interpretation of a Woods SUPERLATTICE."""
        super().test_interpret_woods(ag100_interpreter)
        rpars = ag100_interpreter.rpars
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
        'mirror first neg': (
            'p2 m[-2 1]',
            {'mirror': {(2, -1)}, 'rotation': set(), 'group': 'p2'}
            ),
        'mirror second neg': (
            'p2 m[0 -1]',
            {'mirror': {(0, 1)}, 'rotation': set(), 'group': 'p2'}
            ),
        'cm_rot4': ('cm[1 1] r4',
                    {'mirror': set(), 'rotation': {4}, 'group': 'cm[1 1]'}),
        'pmm': ('pmm',
                {'mirror': set(), 'rotation': set(), 'group': 'pmm'}),
        }
    invalid = {
        'invalid_syntax': 'invalid_syntax',
        'multiple_groups': 'pmm pg',
        'missing_group': '',
        'wrong glide': 'pmm m[3 -2]',
        'wrong screw': 'p4g r(12)',
        }

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid SYMMETRY_BULK."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, interpreter):
        """Ensure invalid SYMMETRY_BULK raises exceptions."""
        self.check_raises(interpreter, val, err.ParameterValueError)


class TestSymmetryEps(_TestInterpretBase):
    """Tests for interpreting SYMMETRY_EPS."""

    param = 'SYMMETRY_EPS'

    def test_interpret_single_value(self, interpreter):
        """Check correct interpretation of EPS value."""
        self.check_assigned(interpreter, '0.5', 0.5)
        rpars = interpreter.rpars
        assert rpars.SYMMETRY_EPS == rpars.SYMMETRY_EPS.z

    def test_interpret_multiple_values(self, interpreter):
        """Check correct interpretation of EPS and EPS.z values."""
        self.interpret(interpreter, '0.1 0.2')
        eps = interpreter.rpars.SYMMETRY_EPS
        assert (float(eps), eps.z) == (0.1, 0.2)

    def test_interpret_invalid_number_of_inputs(self, interpreter):
        """Ensure more than two values raises exceptions."""
        self.check_raises(interpreter, '0.1 0.2 0.3',
                           err.ParameterNumberOfInputsError)

    def test_large_values_log(self, interpreter, caplog, re_match):
        """Check correct interpretation of EPS and EPS.z values."""
        self.interpret(interpreter, '1.5 1.2')
        assert re_match(r'.*[\s\S]*SYMMETRY_EPS.*[\s\S]*very loose constraint',
                        caplog.text)


class TestSymmetryFix(_TestInterpretBase):
    """Tests for interpreting SYMMETRY_FIX."""

    param = 'SYMMETRY_FIX'
    _default = Rparams.get_default(param)

    valid = {'auto': ('t', _default),
             'p1': ('p1', 'p1'),
             'direction': ('cm[1 1]', 'cm[1 1]'),}
    invalid = {'group': ('invalid', err.ParameterParseError),
               'direction_missing': ('cm', err.ParameterParseError),
               'direction_wrong': ('pmt [0 x]', err.ParameterParseError),}

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid SYMMETRY_FIX."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val,exc', invalid.values(), ids=invalid)
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

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid TENSOR_OUTPUT."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val', invalid, ids=invalid)
    def test_interpret_invalid(self, val, interpreter):
        """Ensure invalid TENSOR_OUTPUT raises exceptions."""
        self.check_raises(interpreter, val, err.ParameterParseError)


class TestTheoEnergies(_TestInterpretBase):
    """Tests for interpreting THEO_ENERGIES."""

    param = 'THEO_ENERGIES'
    _defaults = Rparams.get_default(param)
    valid = {
        'single_value_default': ('_', _defaults),
        'single_value': ('5.0', [5.0, 5.0, 1.0]),
        'three_values_all_default': ('_ _ _', _defaults),
        'three_values_one_default': ('1.0 _ 2.0', [1.0, _defaults.stop, 2.0]),
        'three_positive': ('1.0 2.0 0.5', [1.0, 2.0, 0.5]),
        'two_defaults': ('1.0 _ _', [1.0, _defaults.stop, _defaults.step]),
        'range_correction': ('1.1 2.5 0.5', [1.0, 2.5, 0.5]),
        'start rounded negative': ('0.1 2.5 0.5', [0.5, 2.5, 0.5]),
        'start rounded zero': ('0.3 2.5 0.5', [0.5, 2.5, 0.5]),
        }
    invalid = {
        'one_negative': ('1.0 -2.0 0.5', err.ParameterValueError),
        'invalid_range': ('2.0 1.0 0.5', err.ParameterValueError),
        'start_out_of_range': ('-0.5 2.0 0.5', err.ParameterValueError),
        'zero step': ('1.0 2.0 0.0', err.ParameterValueError),
        'too few': ('1.0 2.0 ', err.ParameterNumberOfInputsError),
        'too many': ('1.0 2.0 0.3 9', err.ParameterNumberOfInputsError),
        'nan value': ('nan 2.0 0.3', err.ParameterParseError),
        'nan step': ('0.1 2.0 nan', err.ParameterParseError),
        'inf value': ('inf 2.0 0.3', err.ParameterParseError),
        'inf step': ('0.3 2.5 inf', err.ParameterParseError),
        'float(inf) step': ('0.3 2.5 float("inf")', err.ParameterParseError),
        }

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid THEO_ENERGIES."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid THEO_ENERGIES raises exceptions."""
        self.check_raises(interpreter, val, exc)


class TestTLVersion(_TestInterpretBase):
    """Test for interpreting TL_VERSION"""
    param = 'TL_VERSION'
    valid = {
        '2.0.0': ('v2.0.0', Version(2, 0, 0)),
        'old format': ('1.76', Version(1, 7, 6)),
        'with v': ('v2.0.0', Version(2, 0, 0)),
    }
    invalid = {
        'non version string': ('abc', err.ParameterConversionError)
    }

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid TL_VERSION."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid TL_VERSION raises exceptions."""
        self.check_raises(interpreter, val, exc)


class TestV0Real(_TestInterpretBase):
    """Tests for interpreting V0_REAL."""

    param = 'V0_REAL'
    valid = {
        'rundgren': ('rundgren 1.0 2.0 3.0 4.0', [1.0, 2.0, 3.0, 4.0]),
        'custom': ('EE', 'EEV+workfn'),
        }
    invalid = {
        'too_few': ('rundgren 1.0 2.0', err.ParameterNumberOfInputsError),
        'float': ('rundgren 1 2 3.0 four', err.ParameterFloatConversionError),
        }

    @parametrize('val,expect', valid.values(), ids=valid)
    def test_interpret_valid(self, val, expect, interpreter):
        """Check correct interpretation of valid V0_REAL."""
        self.check_assigned(interpreter, val, expect)

    @parametrize('val,exc', invalid.values(), ids=invalid)
    def test_interpret_invalid(self, val, exc, interpreter):
        """Ensure invalid V0_REAL raises exceptions."""
        self.check_raises(interpreter, val, exc)
