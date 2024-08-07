

.. _cli_bookkeeper:

``viperleed bookkeeper``
========================

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
