"""Tests for _aliases.ini of viperleed.gui.measure._defaults."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-10-02'
__license__ = 'GPLv3+'

import ast
from configparser import ConfigParser

from pytest_cases import fixture

from viperleed.gui.measure.constants import SRC_ALIASES_PATH


@fixture
def alias_config():
    """Fixture that creates a ConfigParser containing the aliases."""
    config = ConfigParser()
    config.read(SRC_ALIASES_PATH)
    return config


class TestDefaultAliases:
    """Tests for the aliases provided in the _defaults."""

    def test_no_multiple_inheritance(self, alias_config):
        """Check that there is no multiple inheritance."""
        # Collect all sections that declare parent_aliases themselves.
        inheriting_sections = (
            section
            for section in alias_config.sections()
            if alias_config.has_option(section, 'parent_aliases')
            )

        # For each section, collect its declared parent_aliases.
        for section in alias_config.sections():
            parents_raw = alias_config.get(section, 'parent_aliases',
                                           fallback=None)
            if not parents_raw:
                continue
            parents = ast.literal_eval(parents_raw)
            # Check if any of the parent_aliases are among
            # the sections that declare parent_aliases.
            for parent in parents:
                assert parent not in inheriting_sections, (
                    f'Invalid multiple inheritance: Section "{section}"'
                    f' tries to inherit from "{parent}", which already '
                    'has its own parent_aliases.'
                    )

    def test_no_capital_letters_in_options_and_values(self, alias_config):
        """Ensure all options and non-parent_aliases values are lowercase."""
        for section in alias_config.sections():
            for option in alias_config.options(section):
                assert option == option.lower(), (
                    f'Option "{option}" in section [{section}] '
                    'contains uppercase letters.'
                    )
                if option == 'parent_aliases':
                    # Skip class names as these contain uppercase strings.
                    continue
                value = alias_config[section][option]
                assert value == value.lower(), (
                    f'Alias "{value}" in [{section}][{option}] '
                    'contains uppercase letters.'
                    )

    def test_no_child_parent_option_overlap(self, alias_config):
        """Ensure that child sections do not redefine parent_aliases."""
        for section in alias_config.sections():
            parents_raw = alias_config.get(section, 'parent_aliases',
                                           fallback=None)
            if not parents_raw:
                continue
            parents = ast.literal_eval(parents_raw)
            # Collect all option names from parents.
            parent_options = set()
            for parent in parents:
                parent_options.update(alias_config.options(parent))

            # Check overlap between child options and parent options.
            child_options = set(alias_config.options(section))
            overlap = parent_options & child_options

            assert not overlap, (
                f'Section [{section}] redefines option(s) already '
                f'defined in its parent_aliases {parents}: {overlap}'
                )
