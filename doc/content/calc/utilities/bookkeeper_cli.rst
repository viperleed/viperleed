.. include:: /substitutions.rst

.. _cli_bookkeeper:

``viperleed bookkeeper``
========================

The command ``viperleed bookkeeper`` manually invokes the
:ref:`bookkeeper<bookkeeper>`.

The bookkeeper runs automatically in ``--clear`` mode before
and in ``--archive`` mode after a calculation. See the
:ref:`bookkeeper<bookkeeper>` page for details.

The |bookkeeper| can also be run manually with ``viperleed bookkeeper``.
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


.. versionchanged:: 0.12.0
    |bookkeeper| now runs automatically in ``--cont`` mode after
    |calc| (in addition to running before |calc|, in "default" mode).

.. versionchanged:: 0.12.0
    The ``--name`` command-line argument was renamed to ``--job-name``.

.. versionadded:: 0.12.0
    ``--history-name`` and ``--work-history-name`` command-line arguments.
    These arguments were removed again in v0.13.0.

.. versionchanged:: 0.13.0
    The behavior of the ``--discard`` mode has been changed. The new equivalent
    mode is ``--discard-full``.
    The ``--cont`` command-line argument was renamed to ``--archive``. The
    behavior was changed to keep both inputs and outputs from a calculation
    in the root folder. What was previously the "default" mode can now be
    executed with ``--clear``. See the :ref:`bookkeeper<bookkeeper>` page
    for more details.

    |bookkeeper| automatically runs (in the same mode) in all domain
    subfolders when executed in the root of a multi-domain calculation.

.. versionremoved:: 0.13.0
    The ``--job-name``, ``--history-name``, and ``--work-history-name``
    command-line arguments.
