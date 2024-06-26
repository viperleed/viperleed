.. include:: /substitutions.rst

.. _command_line_tools:

Command-line tools
==================


.. todo:: viperleed gui

.. todo::
    Look into the sphinx-argparse extension for an easier way to
    automatically produce documentation for command-line tools.

The ViPErLEED Python package ``viperleed`` provides a number of command-line
tools that can be used to run calculations, invoke the
:ref:`bookkeeper<bookkeeper>`, and run the :ref:`utilities<utilities>`
and :ref:`poscar utilities<poscar_utils>`. The ViPErLEED graphical user
interface can also be invoked from the command line via ``viperleed gui``
(currently, only on systems that have access to a graphical display).

.. _cli_calc:

``viperleed``
-------------

    .. argparse::
        :module: viperleed.cli
        :class:  ViPErLEEDMain
        :func:   parser
        :prog:   viperleed
        :nodefaultconst:
        :nosubcommands:
        :usagesection:
        :showusagemain: viperleed

        command
            See also :ref:`cli_calc`,
            ``viperleed gui``,
            :ref:`cli_bookkeeper`,
            :ref:`utilities`,
            and :ref:`poscar_utils`


.. _cli_calc:

``viperleed calc``
------------------

``viperleed calc`` is the main command-line
tool for running |LEED-IV| calculations.

    .. argparse::
        :module: viperleed.cli
        :class:  ViPErLEEDMain
        :func:   parser
        :prog:   viperleed
        :path:   calc
        :nodefaultconst:
        :usagesection:
        :showusagemain: viperleed calc, viperleed.calc

        -v --verbose
            Overwrites the :ref:`LOG_LEVEL` parameter to the value
            ``verbose``.
        -vv --very-verbose
            Overwrites the :ref:`LOG_LEVEL` parameter to the value
            ``vverbose``.
        --tensorleed -t
            If not provided, the path is searched in the environment variable
            :envvar:`VIPERLEED_TENSORLEED`.
        -j --job-name
            See also :ref:`bookkeeper`.
        --no-cont : @after
            :ref:`bookkeeper` is not called manually with ``--cont``
            before the next calculation.
        --history-name : @after
            See also :ref:`bookkeeper`.
        --work-history-name : @after
            See also :ref:`bookkeeper`.
        -w --work
            Default: ``'./work'``


.. _cli_bookkeeper:

``viperleed bookkeeper``
------------------------

The command ``viperleed bookkeeper`` manually invokes the
:ref:`bookkeeper<bookkeeper>`.


The bookkeeper automatically runs in ``archive`` mode before and in
``clear`` mode after a calculation.
See the :ref:`bookkeeper<bookkeeper>` page for details.


The bookkeeper can also be run manually with ``viperleed bookkeeper``.
It can safely be run multiple times.
If no new output is detected, it will simply exit without doing anything.
For details on the different modes, see the :ref:`bookkeeper<bookkeeper>` page.

**Usage:**

.. code-block:: console

    viperleed bookkeeper [options]

**Options:**

- ``-h, --help``: Show a list of all available options and exit.
- ``-a, --archive``: Run in :ref:`archive mode<bookkeeper>`.
- ``-c, --clear``: Run in :ref:`clear mode<bookkeeper>`.
- ``-d, --discard``: Run in :ref:`discard mode<bookkeeper>`.
- ``-df, --discard-full``: Run in :ref:`discard full mode<bookkeeper>`.
- ``-j, --job-name``: Specify a name for the current run.
  Will be appended to the name of the history folder that is created, and is
  logged in history.info
  Passed along to the :ref:`bookkeeper<bookkeeper>`.
- ``--history-name``: Specify the name of the history folder to be used.
  Default is ``history``.
  Passed along to the :ref:`bookkeeper<bookkeeper>`.
- ``--work-history-name``: Specify the name of the work history folder to be
  used.
  Default is ``workhistory``.

  Passed along to the :ref:`bookkeeper<bookkeeper>`.

.. _cli_util_and_poscar:

``viperleed util`` and ``viperleed poscar``
-------------------------------------------

The commands ``viperleed util`` and ``viperleed poscar`` are used to
invoke the ViPErLEED :ref:`utilities<utilities>` and
:ref:`poscar utilities<poscar_utils>` respectively.
See those pages for details.
