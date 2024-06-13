.. include:: /substitutions.rst

.. _command_line_tools:

Command-line tools
==================

.. todo:: viperleed gui

The ViPErLEED Python package ``viperleed`` provides a number of command-line
tools that can be used to run calculations, invoke the
:ref:`bookkeeper<bookkeeper>`, and run the :ref:`utilities<utilities>`
and :ref:`poscar utilities<poscar_utils>`. The ViPErLEED graphical user
interface can also be invoked from the command line via ``viperleed gui``
(currently, only on systems that have access to a graphical display).


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
            Overwrites the :ref:`LOG_LEVEL <log_level>` parameter to the value
            ``verbose``.
        -vv --very-verbose
            Overwrites the :ref:`LOG_LEVEL <log_level>` parameter to the value
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

The bookkeeper runs automatically runs in *default* mode before
and  in *continuation* mode after a calculation. See the
:ref:`bookkeeper<bookkeeper>` page for details.

The bookkeeper can also be run manually with ``viperleed bookkeeper``.
It can safely be run multiple times.
If no new output is detected, it will simply exit without doing anything.

    .. argparse::
        :module: viperleed.cli
        :class:  ViPErLEEDMain
        :func:   parser
        :prog:   viperleed
        :path:   bookkeeper
        :nodefaultconst:
        :usagesection:
        :showusagemain: viperleed bookkeeper, viperleed.calc.bookkeeper

        -c --cont
            See also :ref:`bookkeeper`.
        -d --discard
            See also :ref:`bookkeeper`.


.. _cli_util_and_poscar:

``viperleed util`` and ``viperleed poscar``
-------------------------------------------

The commands ``viperleed util`` and ``viperleed poscar`` are used to
invoke the ViPErLEED :ref:`utilities<utilities>` and
:ref:`poscar utilities<poscar_utils>`, respectively.
See those pages for details.
