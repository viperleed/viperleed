.. aux_to_exp:

=======================
AUXEXPBEAMS_TO_EXPBEAMS
=======================

``AUXEXPBEAMS_TO_EXPBEAMS.py`` is a short python utility script that is packed with ViPErLEED to help useres working with old TensErLEED input data.

It reads a file containing experimental beam data formated in the file format expected by TensErLEED (:ref:`see AUXEXPBEAMS<auxexpbeams>`) and converts it to the (more natural) :ref:`EXPBEAMS.csv format<expbeams>` used and expected by ViPErLEED.
The :ref:`EXPBEAMS.csv file<expbeams>` will be automatically translated to :ref:`see AUXEXPBEAMS<auxexpbeams>` at run-time when required.

Usage
=====

To use the AUXEXPBEAMS_TO_EXPBEAMS utility, simply call the :term:`Python` script ``AUXEXPBEAMS_TO_EXPBEAMS.py`` from the command line (e.g. ``python AUXEXPBEAMS_TO_EXPBEAMS.py``).
The script will prompt the user for the path to the file to be converted and an output filname.