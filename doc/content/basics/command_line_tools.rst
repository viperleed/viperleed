.. _command_line_tools:

.. include:: /substitutions.rst

Command Line Tools
==================

.. todo:: viperleed gui

.. todo::
    Look into the sphinx-argparse extension for an easier way to
    automatically produce documentation for command-line tools.

The ViPErLEED Python package (viperleed) provides a number of command-line
tools that can be used to run calculations, invoke the
:ref:`bookkeeper<bookkeeper>`, and run the :ref:`utilities<utilities>`
and :ref:`poscar utilities<poscar_utils>`.


.. _cli_calc:

``viperleed calc``
------------------

``viperleed calc`` (also |calc|) is the main command line tool
for running calculations.

**Usage:**

.. code-block:: console

    viperleed calc [options]

**Options:**

- ``-h, --help``: Show a list of all available options and exit.
- ``--version``: Show version number and exit.
- ``-v, --verbose``: Increase output verbosity.
  Overwrites the parameter :ref:`LOG_LEVEL <log_level>` to the value
  ``verbose``.
- ``-vv, --very-verbose``: Increase output verbosity further.
  Overwrites the parameter :ref:`LOG_LEVEL <log_level>` to the value
  ``vverbose``.
- ``-w, --workdir``: Specify the working directory of the calculation.
  Default is ``./work``.
- ``-t, --tensorleed``: Specify the path to the tensorleed source code.
  If not provided the path is searched in the environment variable
  ``VIPERLEED_TENSORLEED``.
- ``--delete-workdir``: Delete the working directory of the calculation
  after it has finished.
- ``-n, --name``: Set the name of the system that is being calculated.
- ``--no-cont``: Do not automatically call :ref:`bookkeeper<bookkeeper>`
  after the calculation has finished. Progress may be lost if the
  :ref:`bookkeeper<bookkeeper>` is not called manually with ``--cont``
  before the next calculation.
- ``-j, --job-name``: Specify a name for the current run.
  Will be appended to the name of the history folder that is created,
  and is logged in history.info
  Passed along to the :ref:`bookkeeper<bookkeeper>`.
- ``--history-name``: Specify the name of the history folder to be used.
  Default is ``history``.
  Passed along to the :ref:`bookkeeper<bookkeeper>`.
- ``--work-history-name``: Specify the name of the work history folder to
  be used. Default is ``workhistory``.
  Passed along to the :ref:`bookkeeper<bookkeeper>`.


.. _cli_bookkeeper:

``viperleed bookkeeper``
------------------------

The command ``viperleed bookkeeper`` manually invokes the
:ref:`bookkeeper<bookkeeper>`.

The bookkeeper runs automatically runs in *default* mode before
and  in *continuation* mode after a calculation. See the
:ref:`bookkeeper<bookkeeper>` page for details.

The bookkeeper can also be run manually with ``viperleed bookkeeper``.
It can safely be run multiple times.
If no new output is detected, it will simply exit without doing anything.

**Usage:**

.. code-block:: console

    viperleed bookkeeper [options]

**Options:**

- ``-h, --help``: Show a list of all available options and exit.
- ``-c, --cont``: Run in :ref:`continuation mode<bookkeeper>`.
- ``-d, --discard``: Run in :ref:`discard mode<bookkeeper>`.
- ``-j, --job-name``: Specify a name for the current run.
  Will be appended to the name of the history folder that is created,
  and is logged in history.info. Passed along to the
  :ref:`bookkeeper<bookkeeper>`.
- ``--history-name``: Specify the name of the history folder to be used.
  Default is ``history``.
  Passed along to the :ref:`bookkeeper<bookkeeper>`.
- ``--work-history-name``: Specify the name of the work history folder to
  be used. Default is ``workhistory``.
  Passed along to the :ref:`bookkeeper<bookkeeper>`.

.. _cli_util_and_poscar:

``viperleed util`` and ``viperleed poscar``
-------------------------------------------

The commands ``viperleed util`` and ``viperleed poscar`` are used to
invoke the ViPErLEED :ref:`utilities<utilities>` and
:ref:`poscar utilities<poscar_utils>` respectively.
See those pages for details.
