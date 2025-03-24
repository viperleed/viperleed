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

