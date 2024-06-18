.. include:: /substitutions.rst

.. _api_calc_files:

Files and I/O
=============

The |calc| package provides a number of functions for reading, writing and 
interpreting files. In particular, it also provides the functions to read 
the input files (e.g. :ref:`POSCAR`, :ref:`PARAMETERS`, :ref:`VIBROCC`, 
etc.) and translate the information into :term:`TensErLEED` input files.

.. currentmodule:: viperleed.calc.files

.. autosummary::
    :toctree: files
    :recursive:

    beams
    displacements
    iodeltas
    ioerrorcalc
    iofdopt
    iorfactor
    iosearch
    iosuperpos
    ivplot
    parameters
    patterninfo
    poscar
    searchpdf
    vibrocc
