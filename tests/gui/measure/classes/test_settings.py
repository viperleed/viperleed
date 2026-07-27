"""Tests for module settings of viperleed.gui.measure.classes."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-09-30'
__license__ = 'GPLv3+'

from configparser import ConfigParser
from configparser import NoOptionError
from configparser import NoSectionError
from pathlib import Path

from pytest_cases import fixture
from pytest_cases import parametrize
import pytest

from viperleed.gui.measure.classes.settings import AliasConfigParser
from viperleed.gui.measure.classes.settings import SystemSettings
from viperleed.gui.measure.classes.settings import ViPErLEEDSettings
from viperleed.gui.measure.classes.settings import ensure_aliases_exist
from viperleed.gui.measure.classes.settings import get_aliases_path
from viperleed.gui.measure.classes.settings import interpolate_config_path

_MODULE = 'viperleed.gui.measure.classes.settings'


def test_get_aliases_path(mocker):
    """Test the get_aliases_path function."""
    fake_path = Path('aliases.ini')
    fake_qs = mocker.Mock()
    fake_qs.fileName.return_value = str(fake_path)
    mock_get = mocker.patch(f'{_MODULE}.get_qsettings', return_value=fake_qs)
    result = get_aliases_path()
    assert result == fake_path.resolve()
    mock_get.assert_called_once_with('aliases')


class TestEnsureAliasesExist:
    """Tests for the ensure_aliases_exist function."""

    def test_merges_installed_and_user(self, tmp_path, mocker):
        """Check that old aliases are overwritten and other aliases stay."""
        user_aliases = tmp_path / 'aliases.ini'
        user_aliases.write_text('[Foo]\nstays_user=stays1\nchanged=user')
        installed_aliases = ConfigParser()
        installed_aliases.read_string(
            '[Foo]\nstays_installed=stays2\nchanged=installed'
            )
        defaults = tmp_path / '_defaults'
        defaults.mkdir()
        with (defaults/'_aliases.ini').open('w') as installed_ini:
            installed_aliases.write(installed_ini)
        mocker.patch(f'{_MODULE}.get_aliases_path', return_value=user_aliases)
        mocker.patch(f'{_MODULE}.SRC_ALIASES_PATH', defaults / '_aliases.ini')
        ensure_aliases_exist()
        merged_aliases = ConfigParser()
        merged_aliases.read(user_aliases)
        # pylint: disable=magic-value-comparison
        assert merged_aliases['Foo']['stays_user'] == 'stays1'
        assert merged_aliases['Foo']['stays_installed'] == 'stays2'
        assert merged_aliases['Foo']['changed'] == 'installed'
        # pylint: enable=magic-value-comparison

    def test_new_aliases_written(self, tmp_path, mocker):
        """Check that a user aliases.ini file is always created."""
        user_aliases = tmp_path / 'aliases.ini'
        defaults = tmp_path / '_defaults'
        mocker.patch(f'{_MODULE}.get_aliases_path', return_value=user_aliases)
        mocker.patch(f'{_MODULE}.SRC_ALIASES_PATH', defaults / '_aliases.ini')
        ensure_aliases_exist()
        assert user_aliases.is_file()


class TestInterpolateConfigPath:
    """Tests for the interpolate_config_path function."""

    _no_replace = {
        'no sys path': (None, ('__CONFIG__/foo',)),
        'no __CONFIG__': ('/base', ('nothing/to/replace',)),
        'multiple __CONFIG__': ('/base', ('__CONFIG__/__CONFIG__/bar',)),
        '__CONFIG__ not at start': ('/base', ('other/__CONFIG__/file',)),
        'not a path': ('/base', (3, tuple(), {1: 1})),
        }

    @parametrize('sys_path,files', _no_replace.values(), ids=_no_replace)
    def test_no_replacement(self, sys_path, files, mocker):
        """Check situations that cause no path-segment replacements."""
        fake_qs = mocker.Mock()
        fake_qs.value.return_value = sys_path
        mock_get = mocker.patch(f'{_MODULE}.get_qsettings',
                                return_value=fake_qs)
        no_replacements = list(files)
        interpolate_config_path(no_replacements)
        assert no_replacements == list(files)
        mock_get.assert_called_once_with('measurement')

    def test_success(self, mocker):
        """Check successful interpolation of a path."""
        fake_qs = mocker.Mock()
        fake_qs.value.return_value = '/base'
        mocker.patch(f'{_MODULE}.get_qsettings', return_value=fake_qs)
        # TODO: fails on WindowsPath, because os.fspath gives '\\'
        # separators but we only check for a single '/'
        # filenames = ['__CONFIG__/bar', Path('__CONFIG__/bar')]
        filenames = ['__CONFIG__/bar']
        interpolate_config_path(filenames)
        assert all(str(f).startswith('/base') for f in filenames)


class TestAliasConfigParserNoAlias:
    """Tests for AliasConfigParser sections/options without aliases.

    These tests ensure that AliasConfigParser behaves as a normal
    ConfigParser when accessing sections/options that have no
    corresponding alias.
    """

    @fixture(autouse=True)
    def mock_aliases(self, tmp_path, mocker):
        """Create an example user-aliases file."""
        aliases = tmp_path/'aliases.ini'
        aliases.write_text('''
[SomeClass]
Foo/new_opt=('Foo/old_opt',)
''')
        mocker.patch(f'{_MODULE}.get_aliases_path', return_value=aliases)

    @parametrize(user_fallback=(None, ''))
    def test_get_empty_fallback(self, user_fallback):
        """Check retrieval of an empty fallback without aliases."""
        parser = AliasConfigParser(cls_name='SomeClass')
        value = parser.get('sec', 'opt', fallback=user_fallback)
        assert value is user_fallback

    def test_get_existing_option(self):
        """Check correct retrieval of an existing option value."""
        parser = AliasConfigParser(cls_name='SomeClass')
        explicit_value = 'value'
        parser.read_dict({'Foo': {'opt': explicit_value}})
        assert parser['Foo']['opt'] == explicit_value
        assert parser.get('Foo', 'opt') == explicit_value

    def test_get_fallback(self):
        """Check fallback retrieval for a missing option."""
        parser = AliasConfigParser(cls_name='SomeClass')
        parser.read_dict({'Foo': {}})
        fallback = object()
        assert parser.get('Foo', 'opt', fallback=fallback) is fallback

    def test_missing_option_raises(self):
        """Check complaints for a missing non-aliased option."""
        parser = AliasConfigParser(cls_name='SomeClass')
        parser.read_dict({'Foo': {}})
        with pytest.raises(NoOptionError):
            parser.get('Foo', 'opt')

    def test_missing_section_raises(self):
        """Check complaints for a missing non-aliased section."""
        parser = AliasConfigParser(cls_name='SomeClass')
        parser.read_dict({'Foo': {}})
        with pytest.raises(NoSectionError):
            parser.get('Missing', 'opt')


class TestAliasConfigParser:
    """Tests for AliasConfigParser when aliases are present."""

    @fixture(autouse=True)
    def mock_aliases(self, tmp_path, mocker):
        """Create an example user-aliases file."""
        aliases = tmp_path/'aliases.ini'
        aliases.write_text('''
[WithAliases]
new_sections = ('new_section', 'another_new_section')
new_section/new_option=('oldsection/option','even_older/old_option')
fallback_values = (('A/opt', 'fb'), ('A/opt2', 'pfb'),)

[Foo]
new_sections = ('foo',)
foo/opt = ('old/old',)

[HasParent]
parent_aliases = ('WithAliases', )
fallback_values = (('A/opt2', 'cfb'),)
''')
        mocker.patch(f'{_MODULE}.get_aliases_path', return_value=aliases)

    def test_fallback_from_aliases(self):
        """Check retrieval of a fallback_value from the aliases."""
        parser = AliasConfigParser(cls_name='WithAliases')
        parser.read_dict({'A': {'opt': ''}})
        for option, expect in (('opt', 'fb'), ('opt2', 'pfb')):
            assert parser.get('A', option) == expect
            assert parser['A'][option] == expect

    def test_fallback_from_aliases_with_parent(self):
        """Check retrieval of a fallback_value with parent aliases."""
        parser = AliasConfigParser(cls_name='HasParent')
        parser.read_dict({'A': {'opt': ''}})
        for option, expect in (('opt', 'fb'), ('opt2', 'cfb')):
            assert parser.get('A', option) == expect
            assert parser['A'][option] == expect

    def test_get_from_alias_and_empty_fallback(self):
        """Ensure retrieval of an old-named section/option from new ones."""
        parser = AliasConfigParser(cls_name='Foo')
        old_value = 'aliasval'
        parser.read_dict({'old': {'old': old_value}})
        assert parser.get('foo', 'opt') == old_value
        assert parser['foo']['opt'] == old_value

    def test_iter_aliases(self):
        """Check expected iteration of known aliases."""
        parser = AliasConfigParser(cls_name='WithAliases')
        # pylint: disable-next=protected-access
        aliases = list(parser._iter_aliases('new_section', 'new_option'))
        expect = [['oldsection', 'option'], ['even_older', 'old_option']]
        assert aliases == expect

    def test_multiple_old_files_with_alias_overwrite_dict(self):
        """Ensure aliases persist when multiple files are read."""
        parser = AliasConfigParser(cls_name='WithAliases')

        # first old file
        expect_first = 'first'
        parser.read_dict({'oldsection': {'option': expect_first}})
        assert parser.get('new_section', 'new_option') == expect_first
        assert parser['new_section']['new_option'] == expect_first

        # second old file, expected to overwrite with "second"
        expect_second = 'second'
        parser.read_dict({'even_older': {'old_option': expect_second}})

        assert parser.get('new_section', 'new_option') == expect_second
        assert parser['new_section']['new_option'] == expect_second

    def test_multiple_old_files_with_alias_overwrite_string(self):
        """Ensure aliases persist when multiple files are read."""
        parser = AliasConfigParser(cls_name='WithAliases')

        # first old file
        expect_first = 'first'
        parser.read_string(f'[oldsection]\noption={expect_first}')
        assert parser.get('new_section', 'new_option') == expect_first
        assert parser['new_section']['new_option'] == expect_first

        # second old file, expected to overwrite with "second"
        expect_second = 'second'
        parser.read_string(f'[even_older]\nold_option={expect_second}')

        assert parser.get('new_section', 'new_option') == expect_second
        assert parser['new_section']['new_option'] == expect_second

    def test_new_sections_added(self):
        """Check addition of sections with new names."""
        parser = AliasConfigParser(cls_name='WithAliases')
        added_section = 'new_section'
        assert added_section in parser.sections()

    @parametrize(cls_name=('IHaveNoAliases', '', None))
    def test_no_aliases(self, cls_name):
        """Check emptiness of aliases when none exist."""
        parser = AliasConfigParser(cls_name=cls_name)
        assert not parser._aliases      # pylint: disable=protected-access
        assert not parser._fallbacks    # pylint: disable=protected-access

    def test_old_alias_section_removal(self):
        """Check removal of emptied alias sections."""
        parser = AliasConfigParser(cls_name='Foo')
        old_section = 'old'
        parser.read_dict({old_section: {'old': 'aliasval'}})
        assert old_section not in parser.sections()

    def test_parent_aliases(self):
        """Check if parent aliases are applied."""
        parser = AliasConfigParser(cls_name='HasParent')
        expected = 'expected'
        parser.read_string(f'[oldsection]\noption={expected}')
        assert parser['new_section']['new_option'] == expected

# pylint: disable-next=too-few-public-methods
class TestSystemSettings:
    """Tests for SystemSettings."""

    def test_hidden_folder_and_settings_creation(self, tmp_path, mocker):
        """Test whether the hidden folder and settings were created."""
        fake_path = tmp_path / 'ViPErLEED' / 'Measurement.ini'
        # Patch detected folder path to be in tmp_path.
        mocker.patch('PyQt5.QtCore.QSettings.fileName',
                     return_value=str(fake_path))
        # Patch QSettings.allKeys to force creation of settings folder.
        mocker.patch('PyQt5.QtCore.QSettings.allKeys', return_value=None)
        sys_settings = SystemSettings()
        # pylint: disable-next=protected-access
        settings_path = Path(sys_settings._sys_qsettings.fileName()).resolve()
        assert settings_path.parent.is_dir()
        assert settings_path.is_file()


class TestCheckMandatorySettings:
    """Tests for _check_mandatory_settings behavior."""

    @fixture
    def settings(self, tmp_path, mocker):
        """Create a mock SystemSettings instance pointing to a tmp file."""
        fake_path = tmp_path / 'settings.ini'
        fake_qs = mocker.Mock()
        fake_qs.fileName.return_value = str(fake_path)
        fake_qs.childGroups.return_value = []
        fake_qs.childKeys.return_value = []
        fake_qs.allKeys.return_value = []
        fake_qs.value.return_value = ''

        mocker.patch(f'{_MODULE}.get_qsettings', return_value=fake_qs)
        sys_settings = SystemSettings()
        mocker.patch.object(sys_settings, 'update_file')
        return sys_settings

    # pylint: disable=protected-access
    def test_auto_fills_missing_non_null_settings(self, settings):
        """Check that missing settings are auto-created."""
        # Define requirements
        settings._SystemSettings__non_null = [('SecA', 'opt1')]
        settings._SystemSettings__non_mandatory = [('SecB',)]
        settings._SystemSettings__mandatory = []
        settings._check_mandatory_settings()

        # Check auto-creation
        assert settings.has_section('SecA')
        assert settings.has_option('SecA', 'opt1')
        assert settings.get('SecA', 'opt1') == ''   # pylint: disable=C1804
        assert settings.has_section('SecB')
        settings.update_file.assert_called_once()

    def test_raises_runtime_error_on_missing_mandatory(self, settings):
        """Check that RuntimeError is raised when mandatory settings fail."""
        settings._SystemSettings__non_null = []
        settings._SystemSettings__non_mandatory = []
        settings._SystemSettings__mandatory = [('MandatorySec', 'opt')]

        with pytest.raises(RuntimeError, match='MandatorySec/opt'):
            settings._check_mandatory_settings()
    # pylint: enable=protected-access


class TestMissesSettings:
    """Tests for the misses_settings method."""

    def test_all_settings_valid(self):
        """Check that no invalid settings are returned when all exist."""
        parser = ViPErLEEDSettings(cls_name='Foo')
        parser.read_dict({'Sec': {'opt1': 'val1', 'opt2': 'a'}})

        invalid = parser.misses_settings(
            ('Sec',),
            ('Sec', 'opt1'),
            ('Sec', 'opt2', ['a', 'b'])
        )
        assert not invalid

    def test_missing_section(self):
        """Check reporting when a section is missing."""
        parser = ViPErLEEDSettings(cls_name='Foo')
        invalid = parser.misses_settings(('MissingSec',))
        assert invalid == ['MissingSec']

    def test_missing_option(self):
        """Check reporting when an option within a section is missing."""
        parser = ViPErLEEDSettings(cls_name='Foo')
        parser.read_dict({'Sec': {}})
        invalid = parser.misses_settings(('Sec', 'MissingOpt'))
        assert invalid == ['Sec/MissingOpt']

    def test_invalid_admissible_value(self):
        """Check reporting when an option is not in admissible_values."""
        parser = ViPErLEEDSettings(cls_name='Foo')
        parser.read_dict({'Sec': {'opt': 'wrong_val'}})
        invalid = parser.misses_settings(('Sec', 'opt', ['val1', 'val2']))
        assert invalid == ['Sec/opt not one of val1, val2']

    @parametrize('bad_setting', [(), ('a', 'b', 'c', 'd')])
    def test_invalid_setting_format_raises(self, bad_setting):
        """Check that setting tuples with bad length raise TypeError."""
        parser = ViPErLEEDSettings(cls_name='Foo')
        with pytest.raises(TypeError):
            parser.misses_settings(bad_setting)
