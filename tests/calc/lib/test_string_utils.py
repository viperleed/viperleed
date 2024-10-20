"""Tests for module string_utils of viperleed.calc.lib."""


__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-06-13'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.calc.lib.string_utils import parent_name
from viperleed.calc.lib.string_utils import range_to_str
from viperleed.calc.lib.string_utils import read_int_line
from viperleed.calc.lib.string_utils import read_int_range
from viperleed.calc.lib.string_utils import read_vector
from viperleed.calc.lib.string_utils import rsplit_once
from viperleed.calc.lib.string_utils import split_string_range
from viperleed.calc.lib.string_utils import strip_comments
from viperleed.calc.lib.string_utils import to_snake_case


class TestParentName:
    """Tests for the parent_name function."""

    _valid = {
        'remove None': ('a.b.c', None, 'a.b'),
        'empty remove': ('a.b.c', '', 'a.b'),
        'tail two': ('a.b.c', 'b.c', 'a'),
        'tail two dotted': ('a.b.c', '.b.c', 'a'),
        'middle': ('a.b.c', 'b', 'a'),
        'middle dotted': ('a.b.c', '.b', 'a'),
        'four': ('a.b.c.d', 'c.d', 'a.b'),
        'not in name': ('a.b.c.d', 'e', 'a.b.c.d'),
        'not dotted': ('a', '', 'a'),
        'not dotted, not in name': ('a', 'b', 'a'),
        }

    @parametrize('string,remove,expect', _valid.values(), ids=_valid)
    def test_valid(self, string, remove, expect):
        """Check simple behavior with a custom remove kwarg."""
        assert parent_name(string, remove=remove) == expect

    _invalid = {
        'non-string dotted_name': ((['a.b.c'], None), AttributeError),
        'non-string remove': (('a.b.c', 123), AttributeError),
        }

    @parametrize('args,exc', _invalid.values(), ids=_invalid)
    def test_raises(self, args, exc):
        """Check complaints for invalid types."""
        with pytest.raises(exc):
            parent_name(*args)


class TestRangeToStr:
    """Tests for the range_to_str function."""

    _valid = {
        'one element': ([1], '1'),
        'empty': (tuple(), ''),
        'generator': ((i for i in range(10) if i%4), '1-3, 5-7, 9'),
        'sorted list': ([1, 2, 3, 5, 6, 8], '1-3, 5-6, 8'),
        'unsorted list': ([3, 1, 2, 8, 5, 6], '1-3, 5-6, 8'),
        'unsorted duplicate list': ([8, 3, 1, 2, 5, 5, 8, 5, 6],
                                    '1-3, 5-6, 8'),
        'one-element group': ([1, 2, 4, 6, 7, 8], '1-2, 4, 6-8'),
        }
    _fail = {  # No known failing conditions so far
        }

    @parametrize('integers,expect', _valid.values(), ids=_valid)
    def test_valid(self, integers, expect):
        """Check expected outcome."""
        assert range_to_str(integers) == expect

    @parametrize('integers,expect', _fail.values(), ids=_fail)
    @pytest.mark.xfail(reason='Current implementation does not support this')
    def test_fail(self, integers, expect):
        """Check conditions that we do not handle correctly."""
        self.test_valid(integers, expect)

    _raises = {
        'string item': ([1, '2', 3], {}, TypeError),
        'string sequence': ('831255856', {}, TypeError),
        'invalid separator': ([1,2,3], {'sep': None}, AttributeError),
        }

    @parametrize('iterable,kwargs,exc', _raises.values(), ids=_raises)
    def test_raises(self, iterable, kwargs, exc):
        """Check complaints when using an invalid iterable."""
        with pytest.raises(exc):
            range_to_str(iterable, **kwargs)


class TestReadIntLine:
    """Tests for the read_int_line function."""

    _valid = {  # args, expected result
        ('123456', 3): (123, 456),
        ('123456   \n\t', 3): (123, 456),
        ('123456', 2): (12, 34, 56),
        ('', 3): (),
        ('12345', 3): (123, 45),
        ('12 34 5', 3): (12, 34, 5),
        }

    @parametrize('args,expect', _valid.items())
    def test_valid(self, args, expect):
        """Check expected outcome."""
        assert read_int_line(*args) == expect

    def test_raises(self):
        """Check complaints if items are not integers."""
        with pytest.raises(ValueError):
            read_int_line('1 2 3 a 5', 2)


class TestReadIntRange:
    """Tests for the read_int_range function."""

    _valid = {
        '1-3 5 7-9': [1, 2, 3, 5, 7, 8, 9],
        '1 3 5': [1, 3, 5],
        '1-3  5  7-9': [1, 2, 3, 5, 7, 8, 9],
        '1-3 2 3 5-8': [1, 2, 3, 5, 6, 7, 8],
        '1-3 2 3 5:8': [1, 2, 3, 5, 6, 7, 8],
        '1-1': [1],
        }
    _fail = {  # string, expected, and actual results
        'not a valid range': ('1-3-5', [], [1, 2, 3, 4, 5]),
        }

    @parametrize('string,expect', _valid.items())
    def test_valid(self, string, expect):
        """Check expected outcome."""
        assert read_int_range(string) == expect

    @parametrize('string,what_should_be,what_is', _fail.values(), ids=_fail)
    def test_fail(self, string, what_should_be, what_is):
        """Check currently incorrect results."""
        result = read_int_range(string)
        assert result == what_is
        assert result != what_should_be

    _invalid = {
        'not an int': '1-3 a 5',
        'not an int bound': '1-3 a-12 5',
        'partial range': '1- 3-6 9',
        'empty range': '\n\t  ',
        }

    @parametrize(string=_invalid.values(), ids=_invalid)
    def test_invalid(self, string):
        """Check expected outcome."""
        with pytest.raises(ValueError):
            read_int_range(string)


class TestReadVector:
    """Tests for the read_vector function."""

    _valid = {
        'xyz': (('xyz[1   2 3 ]',), [1, 2, 3]),
        'abc': (('abc[1 2 3]', np.eye(3)), [1, 2, 3]),
        'abc, non-simple ucell': (
            ('abc[1 1 1]',
             np.array([[1, 2, 3], [0, 1, 4], [0, 0, 1]])),
            [1, 3, 8],
            )
        }

    @parametrize('args,expect', _valid.values(), ids=_valid)
    def test_valid(self, args, expect):
        """Check expected outcome for acceptable arguments."""
        result = read_vector(*args)
        assert result == pytest.approx(expect)

    _invalid = {
        '[ ]': ValueError,         # No items
        'xyz[1 2]': ValueError,    # Not enough items
        'abc[1 2 3]': TypeError,   # Missing ucell
        'xyz[a b c]': ValueError,  # Non-numeric items
        'xyz[1 2 ?]': ValueError,  # Non-numeric items
        'pxyz[1 2 3]': ValueError,      # Invalid prefix
        'xyz[1 2.a 3]': ValueError,     # Non-numeric items
        'xyz[f1 f2 f3]': ValueError,    # Non-numeric items
        'xyz[0.1.2 3 4]': ValueError,   # Non-numeric items
        'xyz 1.0 2.0 3.0': ValueError,  # No brackets
        # TODO: The next one currently passes. Do we want to complain?
        # 'xyz[1.0 2.0 3.0] extra': ValueError,
        }

    @parametrize('string,exc', _invalid.items(), ids=_invalid)
    def test_invalid(self, string, exc):
        """Check outcome for invalid arguments."""
        with pytest.raises(exc):
            read_vector(string)


class TestRsplitOnce:
    """Tests for the rsplit_once function."""

    _valid = {
        'empty string': ('', '.', ['']),
        'multi sep': ('a.b.c.d.e', '.', ['a.b.c.d', 'e']),
        'multi sep, consecutive': ('a.a....a', '.', ['a.a...', 'a']),
        'multi sep, same chars': ('a.a.a.a.a.a.a', '.', ['a.a.a.a.a.a', 'a']),
        'no sep': ('abcdefghi', '.', ['abcdefghi']),
        'one sep': ('abc.def', '.', ['abc', 'def']),
        'two sep': ('abc.def.ghi', '.', ['abc.def', 'ghi']),
        'sep at end': ('abcdef.', '.', ['abcdef', '']),
        'sep at start': ('.abcdef', '.', ['', 'abcdef']),
        'sep long': ('abcXXdefXXghi', 'XX', ['abcXXdef', 'ghi']),
        'sep only': ('.', '.', ['', '']),
        'sep unicode': ('abc.def.ghi.jkl', '。', ['abc.def.ghi.jkl']),
        'special chars': ('a$b^c.d!e@f#g', '#', ['a$b^c.d!e@f', 'g']),
        }

    @parametrize('string,sep,expect', _valid.values(), ids=_valid)
    def test_valid(self, string, sep, expect):
        """Check correct outcome with valid arguments."""
        assert rsplit_once(string, sep) == expect

    _raises = {
        'not a string': (None, ',', TypeError),
        'not a string sep': ('abc,def', 1, TypeError),
        'empty sep': ('abc,def', '', ValueError),
        }

    @parametrize('string,sep,exc', _raises.values(), ids=_raises)
    def test_raises(self, string, sep, exc):
        """Check complaints for invalid inputs."""
        with pytest.raises(exc):
            rsplit_once(string, sep)


class TestSplitStringRange:
    """Tests for the split_string_range function."""

    _valid = {
        '1-5': ('1', '5'),
        '2-2': ('2', '2'),
        '1 -   5': ('1 ', '   5'),
        '1:5': ('1', '5'),
        'start:stop': ('start', 'stop'),
        }
    _fail = { # Should these all raise some form of syntax error?
        'partial range': ('1-', ('1', '')),
        'multi dashes': ('1--5', ('1', '5')),
        'multi colons': ('1::5', ('1', '5')),
        'three items increasing': ('1:5:10', ('1', '10')),
        'three items non monotonic': ('1:5:2', ('1', '2')),
        'dash only': ('-', ('', '')),
        'colon only': (':', ('', '')),
        }

    @parametrize('string,expect', _valid.items())
    def test_valid(self, string, expect):
        """Check expected outcome."""
        assert split_string_range(string) == expect

    @parametrize('string,current_result', _fail.values(), ids=_fail)
    def test_fail(self, string, current_result):
        """Check currently incorrect results."""
        result = split_string_range(string)
        assert result == current_result

    _raises = {
        'empty': '',
        'only spaces': '   ',
        'no separator': '1',
        }

    @parametrize(string=_raises.values(), ids=_raises)
    def test_raises(self, string):
        """Check complaints when using an invalid iterable."""
        with pytest.raises(ValueError):
            split_string_range(string)


class TestStripComments:
    """Tests for the strip_comments function."""

    _valid = {
        ('code!and comment with exclamation mark',): 'code',
        ('code #and comment with hash',): 'code',
        ('code %   comment with percent',): 'code',
        ('  code with spaces   !comment',): 'code with spaces',
        ('  code with spaces   !comment', False): '  code with spaces   ',
        ('no comments',): 'no comments',
        (' code!with_more_than_one_comment#character',): 'code',
        ('#comment only ! with multiple % comment characters',): '',
        ('!comment only # with multiple % comment characters',): '',
        ('%comment only ! with multiple # comment characters',): '',
        }

    @parametrize('args,expect', _valid.items(), ids=(str(k) for k in _valid))
    def test_valid(self, args, expect):
        """Check expected outcome of stripping comments."""
        assert strip_comments(*args) == expect

    _raises = {
        'list': ['test'],
        'tuple': tuple('test'),
        'None': None,
        }

    @parametrize(invalid_line=_raises.values(), ids=_raises)
    def test_raises(self, invalid_line):
        """Check complaints when an invalid type is given."""
        with pytest.raises(TypeError):
            strip_comments(invalid_line)


class TestToSnakeCase:
    """Tests for the to_snake_case function."""

    _valid = {
        'PascalCase': 'pascal_case',
        'camelCaseExample': 'camel_case_example',
        'Word': 'word',
        'Test3': 'test3',
        'lowercase': 'lowercase',
        'already_snake_case': 'already_snake_case',
        'A': 'a',
        'CamelC': 'camel_c',
        'HTTPResponse': 'http_response',
        'CamelCase35': 'camel_case35',
        'HTTP2Response': 'http2_response',
        'HTTP25ResponseForNumber9': 'http25_response_for_number9',
        'Version2Update3': 'version2_update3',
        '': '',
        'URLProcessor': 'url_processor',
        'HTTPRequestHandler': 'http_request_handler',
        'HTTP2REQUESTHANDLER': 'http2requesthandler',  # Not super
        '_SOME_MODULE_CONSTANT': '_some_module_constant',
        'Update2024Version3': 'update2024_version3',
        }

    @parametrize('name,expect', _valid.items(), ids=_valid)
    def test_valid(self, name, expect):
        """Check correct conversion to_snake_case."""
        assert to_snake_case(name) == expect

    _raises = {
        'mixed_Case': NotImplementedError,
        'Pascal__Case': NotImplementedError,
        'notPascal_case': NotImplementedError,
        'PascalCase_with_error': NotImplementedError,
        'mixedCaseExample_withError': NotImplementedError,
        'Pascal$Case': ValueError,
        'Pascal@Case': ValueError,
        'camelCase_WithPascalCase': NotImplementedError,
        None: TypeError,
        (1, 2, 3): TypeError,
        }

    @parametrize('name,exc', _raises.items())
    def test_raises(self, name, exc):
        """Check complaints for invalid inputs."""
        with pytest.raises(exc):
            to_snake_case(name)
