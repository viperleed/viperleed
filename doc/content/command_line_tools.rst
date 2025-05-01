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

.. _cli_viperleed:

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
        -w --work
            Default: ``./work``

.. versionadded:: 0.12.0
    The ``--all-tensors``, ``--delete-workdir``, ``--history-name``,
    ``--job-name``, ``--name``, ``--no-cont``, ``--verbose``,
    ``--very-verbose``, and ``--work-history-name`` arguments. In earlier
    versions, the functionality of ``--all-tensors`` and ``--delete-workdir``
    was available by editing the job script.

.. versionchanged:: 0.12.0
    Renamed the ``--source`` argument to ``--tensorleed``.
    Modified the behavior of the ``--work`` argument. In earlier versions,
    calculations would run in a :file:`work` subfolder of the path specified
    via the ``--work`` argument. Now calculations run in the path given as
    the ``--work`` argument.

.. versionremoved:: 0.13.0
    The ``--no-cont`` argument. The same effect can be obtained by manually
    running :ref:`bookkeeper` in ``--discard`` (or ``--discard-full``) mode
    after |calc|.

.. versionremoved:: 0.13.0
    The ``--job-name``, ``--history-name``, and ``--work-history-name``
    arguments.

.. versionchanged:: 0.13.0
    Removed the ``--delete-workdir`` argument, whose default was ``False``.
    It is replaced by ``--keep-workdir``.
