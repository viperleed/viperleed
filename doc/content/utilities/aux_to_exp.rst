.. _aux_to_exp:

=============
beams utility
=============

The beams utility is a Python script that is packed with ViPErLEED to help users
working with experimental input data formatted for other programs.
It supports conversion from the formats used by :term:`TensErLEED` 
(:ref:`see AUXEXPBEAMS<auxexpbeams>`) and
`SATLEED <https://www.icts.hkbu.edu.hk/VanHove_files/leed/leedpack.html>`.
The utility converts the input data to the comma-separate
:ref:`EXPBEAMS.csv format<expbeams>` used and expected by ViPErLEED.


Usage
=====

The beams is installed with the ViPErLEED package by default.

It is invoked with the following commands:

.. code-block:: bash

    $ viperleed util beams from-TensErLEED
    $ viperleed util beams from-SATLEED

For more information on usage, see

.. code-block:: bash

    $ viperleed util beams from-TensErLEED --help
    $ viperleed util beams from-SATLEED --help
