.. _log_level:

LOG_LEVEL
=========

The parameter LOG_LEVEL controls the verbosity of the output log.
LOG_LEVEL follows the standard Python logging levels, as defined by the
`logging module <https://docs.python.org/3/library/logging.html>`__.
As such, higher lower values correspond to higher verbosity with any value
below 10 being debug messages. In :ref:`PARAMETERS` LOG_LEVEL can be assigned 
an integer value in the range [0, 50], or one of the following shorthands:

- ``debug``, ``True``: 10
- ``verbose``, ``v``: 5
- ``vverbose``, ``vv``: 1

The ``verbose`` and ``very-verbose`` options may lead to a very larger
amount of output, which can be useful for debugging purposes.

**Default**: LOG_LEVEL = 20 (shows ``INFO`` messages and higher)

**Syntax**:

::

   LOG_LEVEL = 10          ! show DEBUG level messages and higher
   LOG_LEVEL = debug       ! sames as above

**Acceptable values**: :math:`0 \le` LOG_LEVEL :math:`\le 50`

.. note::
    If a calculation is invoked with the ``--verbose`` or ``--very-verbose``
    command-line options, the log level will be overridden to 5 or 1,
    respectively. This can be useful for quick debugging without having
    to change the PARAMETERS file.

.. versionadded:: 0.10

   Prior to version 0.10, this parameter was called LOG_DEBUG and only had
   the options ``True``/``False``. ``LOG_DEBUG = True`` is interpreted as
   an alias for ``LOG_LEVEL = debug``.