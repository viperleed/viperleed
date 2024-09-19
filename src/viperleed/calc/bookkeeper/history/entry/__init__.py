"""Package entry of viperleed.calc.bookkeeper.history.

Defines the functionality that pertains to each separate
'block' of a history.info file.

Modules
-------
entry
    Defines the HistoryInfoEntry and PureCommentEntry classes
    representing an entry with and without any known fields,
    respectively.
enums
    Enumeration classes used in multiple spots in the package.
field
    Base classes for handling single and multiple (related)lines
    in an entry.
field_collection
    Classes representing collections of fields.
list_of_int_field
    Classes handling fields whose value is a list of integers.
notes_field
    Classes handling the user-notes field of an entry.
rfactor_field
    Classes handling lines of an entry carrying information
    about R factors.
string_field
    Classes for handling simple string-only fields of entries.
time_field
    Classes handling the TIME field of an entry.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-25'
__license__ = 'GPLv3+'
