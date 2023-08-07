.. _aux_to_exp:

=======================
AUXEXPBEAMS_TO_EXPBEAMS
=======================

``AUXEXPBEAMS_to_EXPBEAMS`` is an interactive Python :ref:`utility<utilities>` that is packed with ViPErLEED to help users working with old TensErLEED input data.

It reads a file containing experimental beam data formatted in the file format expected by TensErLEED (:ref:`see AUXEXPBEAMS<auxexpbeams>`) and converts it to the (more natural) comma-separated :ref:`EXPBEAMS.csv format<expbeams>` used and expected by ViPErLEED.

Usage
=====

``AUXEXPBEAMS_to_EXPBEAMS`` is installed with the ViPErLEED package by default.

It is invoked with the following command:

.. code-block:: bash

    $ viperleed util AUXEXPBEAMS_to_EXPBEAMS