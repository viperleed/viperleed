"""Tests for module string_field of viperleed.calc.bookkeeper.history.entry."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-09-03'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.history.entry.enums import FieldTag
from viperleed.calc.bookkeeper.history.entry.field import EmptyField
from viperleed.calc.bookkeeper.history.entry.field import MissingField
from viperleed.calc.bookkeeper.history.entry.string_field import FolderField
from viperleed.calc.bookkeeper.history.entry.string_field import JobNameField
from viperleed.calc.bookkeeper.history.entry.string_field import StringField
from viperleed.calc.bookkeeper.history.errors import FixableSyntaxError

from .....helpers import not_raises
from .test_field import _TestFieldUtils


class TestFolderField(_TestFieldUtils):
    """Tests for the tag-bearing FolderField class."""

    test_cls = FolderField

    @fixture(name='folder_field')
    def fixture_folder_field(self, make_field_factory):
        """Return a factory of FolderField instances."""
        return make_field_factory(self.test_cls)

    def test_class_attrs(self):
        """Check the expected class-level attributes of FolderField."""
        field = self.test_cls
        assert field.is_mandatory
        assert field.tag is FieldTag.FOLDER

    _init = {
        'empty': (
            '',
            {'value': EmptyField,
             'job_name': None, 'is_missing': False,
             'is_empty': True, 'was_understood': False},
            ),
        'invalid': (
            'invalid_format',
            {'value': 'invalid_format',
             'job_name': None, 'is_missing': False,
             'is_empty': False, 'was_understood': False},
            ),
        'job': (
            't123.r456_789101-121314_jobname',
            {'value': 't123.r456_789101-121314_jobname',
             'job_name': 'jobname', 'is_missing': False,
             'was_understood': True, 'needs_fixing': False},
            ),
        'job, hyphen': (
            't123.r456_230101-123456_job-1-23-b_moved-123456',
            {'value': 't123.r456_230101-123456_job-1-23-b_moved-123456',
             'job_name': None, 'was_understood': False},
            ),
        'job, invalid chars': (
            't123.r456_789101-121314_!@#$',
            {'value': 't123.r456_789101-121314_!@#$',
             'job_name': None, 'was_understood': False},
            ),
        'job, white spaces': (
            't123.r456_789101-121314_job name  ',
            {'value': 't123.r456_789101-121314_job name',
             'job_name': 'job name', 'was_understood': True},
            ),
        'leading space': (
            '   t123.r456_789101-121314_my_job',
            {'value': 't123.r456_789101-121314_my_job', 'job_name': 'my_job'},
            ),
        'long job': (
            't123.r456_230101-123456_' + 'a'*1000,
            {'value': 't123.r456_230101-123456_' + 'a'*1000,
             'job_name': 'a'*1000},
            ),
        'missing bits': (
            't9r9',
            {'value': 't9r9', 'was_understood': False},
            ),
        'moved once, missing time, no job': (
            't123.r456_230101-123456_moved-123456',
            {'value': 't123.r456_230101-123456_moved-123456',
             'job_name': None, 'was_understood': False},
            ),
        'multiple jobs': (
            't123.r456_789101-121314_job1_job2',
            {'value': 't123.r456_789101-121314_job1_job2',
             'job_name': 'job1_job2'},
            ),
        'no job': (
            't123.r456_789101-121314',
            {'value': 't123.r456_789101-121314',
             'job_name': None, 'is_missing': False,
             'was_understood': True, 'needs_fixing': False},
            ),
        'only optional parts': (
            '_moved-230101-123456_job',
            {'value': '_moved-230101-123456_job', 'was_understood': False},
            ),
        'suffix, moved': (
            't123.r456_moved-789101-121314',
            {'value': 't123.r456_moved-789101-121314',
             'job_name': None, 'was_understood': True},
            ),
        'suffix, moved, job': (
            't123.r456_moved-789101-121314_job',
            {'value': 't123.r456_moved-789101-121314_job',
             'job_name': 'job', 'was_understood': True},
            ),
        'suffix, moved, time conflict': (
            't123.r456_moved-789101-121314_moved-123456-789012',
            {'value': 't123.r456_moved-789101-121314_moved-123456-789012',
             'job_name': None, 'was_understood': True},
            ),
        'suffix, moved, time conflict, job': (
            't123.r456_moved-789101-121314_job_moved-123456-789012',
            {'value': 't123.r456_moved-789101-121314_job_moved-123456-789012',
             'job_name': 'job', 'was_understood': True},
            ),
        'suffix, search': (
            't123.r456.999_789101-121314_job',
            {'value': 't123.r456.999_789101-121314_job',
             'job_name': 'job', 'is_missing': False,
             'was_understood': True, 'needs_fixing': False},
            ),
        'suffix, time conflict': (
            't123.r456_789101-121314_moved-123456-789012',
            {'value': 't123.r456_789101-121314_moved-123456-789012',
             'job_name': None, 'was_understood': True},
            ),
        'suffix, time conflict, job': (
            't123.r456_789101-121314_job_moved-123456-789012',
            {'value': 't123.r456_789101-121314_job_moved-123456-789012',
             'job_name': 'job', 'was_understood': True},
            ),
        'short tensor and run, zeros': (
            't1.r1_000000-000000',
            {'value': 't1.r1_000000-000000', 'was_understood': True},
            ),
        'trailing space': (
            't123.r456_789101-121314_my_job ',
            {'value': 't123.r456_789101-121314_my_job',
             'job_name': 'my_job'},
            ),
        'wrong date-time': (
            't123.r456_23-01-01_123456',
            {'value': 't123.r456_23-01-01_123456',
             'was_understood': False},
            ),
        }

    @parametrize('value,attrs', _init.values(), ids=_init)
    def test_init(self, value, attrs, folder_field):
        """Check initialization."""
        self.check_attrs(folder_field, attrs, value)

    _other_job = (
        'no one has this',
        JobNameField('no one has this'),
        )

    @parametrize('value,attrs', _init.values(), ids=_init)
    @parametrize(other_job=_other_job)
    def test_check_has_job_name_different(self, value, attrs, other_job,
                                          folder_field):
        """Check complaints for inconsistent job names."""
        field = folder_field(value)
        context = (pytest.raises if attrs.get('was_understood', True)
                   else not_raises)
        with context(FixableSyntaxError) as exc_info:
            field.check_has_job_name(other_job)
        if not exc_info:
            return
        fixed_value = exc_info.value.fixed_value
        try:
            other_job_str = other_job.value
        except AttributeError:
            other_job_str = other_job
        assert f'_{other_job_str}' in fixed_value
        assert field.needs_fixing
        fixed_field = folder_field(fixed_value)
        assert fixed_field.job_name == other_job_str
        if field.job_name:
            assert field.job_name not in fixed_value
        assert str(None) not in fixed_value

    @parametrize('value,_', _init.values(), ids=_init)
    def test_check_has_job_name_consistent(self, value, _, folder_field):
        """Check there are no complaints for consistent job names."""
        field = folder_field(value)
        if not field.job_name:
            return
        with not_raises(FixableSyntaxError):
            field.check_has_job_name(field.job_name)

    _invalid_job = {
        'not a string': [],
        'none': None,
        'empty': '',
        'missing field': JobNameField(),
        'hyphens': 'job-cannot-have-hyphens',
        'special chars': 'job&can$only@contain@alphanumeric=,_,and spaces',
        }

    @parametrize(invalid_job=_invalid_job.values(), ids=_invalid_job)
    def test_check_has_job_name_raises(self, invalid_job, folder_field):
        """Check complaints when an invalid job_name is used for checking."""
        field = folder_field('t123.r456_000007-000008')
        with pytest.raises(ValueError):
            field.check_has_job_name(invalid_job)


class TestJobNameField(_TestFieldUtils):
    """Tests for the tag-bearing JobNameField class."""

    test_cls = JobNameField

    @fixture(name='job_field')
    def fixture_job_field(self, make_field_factory):
        """Return a factory of JobNameField instances."""
        return make_field_factory(self.test_cls)

    def test_class_attrs(self):
        """Check the expected class-level attributes of JobNameField."""
        field = self.test_cls
        assert not field.is_mandatory
        assert field.tag is FieldTag.JOB_NAME

    _init = {
        'empty': (
            '',
            {'value': EmptyField, 'is_missing': False, 'is_empty': True,
             'was_understood': False},
            ),
        'hyphens': (
            'job-name123',
            {'value': 'job-name123', 'was_understood': False},
            ),
        'invalid type': (
            [],
            {'value': [], 'is_missing': False, 'is_empty': False,
             'was_understood': False},
            ),
        'leading spaces': (
            '   jobname',
            {'value': 'jobname'},
            ),
        'long': (
            'a' * 257,
            {'value': 'a' * 257},
            ),
        'special characters': (
            '$myjob@123#!',
            {'value': '$myjob@123#!', 'was_understood': False},
            ),
        'valid': (
            'jobname',
            {'value': 'jobname', 'is_missing': False, 'is_empty': False,
             'was_understood': True},
            ),
        }

    @parametrize('value,attrs', _init.values(), ids=_init)
    def test_init(self, value, attrs, job_field):
        """Check initialization."""
        self.check_attrs(job_field, attrs, value)


class TestStringField(_TestFieldUtils):
    """Tests for the tag-less StringField class."""

    test_cls = StringField

    @fixture(name='string_field')
    def fixture_string_field(self, make_concrete_field_instance):
        """Return a concrete StringField subclass."""
        return make_concrete_field_instance(self.test_cls)

    _init = {
        'empty': (
            '',
            {'value': EmptyField, '_value_str': EmptyField,
             'is_missing': False, 'is_empty': True, 'was_understood': False},
            ),
        'hyphens': (
            'job-name123',
            {'value': 'job-name123', 'was_understood': True},
            ),
        'invalid type': (
            123,
            {'value': 123, 'was_understood': False},
            ),
        'leading and trailing spaces': (
            '   test   ',
            {'value': 'test', 'was_understood': True},
            ),
        'leading spaces': (
            ' \t    Leading spaces',
            {'value': 'Leading spaces', 'was_understood': True},
            ),
        'long string': (
            'a' * 1000,
            {'value': 'a' * 1000, 'was_understood': True},
            ),
        'missing': (
            MissingField,
            {'is_missing': True, '_value_str': None, 'was_understood': True},
            ),
        'no spaces': (
            'NoSpacesAnywhere',
            {'value': 'NoSpacesAnywhere', 'was_understood': True},
            ),
        'none is empty': (
            None,
            {'value': EmptyField, '_value_str': EmptyField,
             'was_understood': False},
            ),
        'not empty': (
            'Test string',
            {'value': 'Test string', 'is_missing': False, 'is_empty': False,
             'was_understood': True},
            ),
        'special characters': (
            'test@#^&*!',
            {'value': 'test@#^&*!', 'was_understood': True},
            ),
        'trailing spaces': (
            'Trailing spaces\t    \t ',
            {'value': 'Trailing spaces', 'was_understood': True},
            ),
        'white space only': (  # Treated like empty
            '   ',
            {'value': EmptyField, '_value_str': EmptyField, 'is_empty': True,
             'is_missing': False, 'was_understood': False},
            ),
        }

    @parametrize('value,attrs', _init.values(), ids=_init)
    def test_init(self, value, attrs, string_field):
        """Check initialization."""
        field = self.check_attrs(string_field, attrs, value)
        str_expect = field.value if not field.is_missing else None
        # pylint: disable-next=protected-access           # OK in tests
        assert field._value_str is str_expect
