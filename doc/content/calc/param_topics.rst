.. include:: /substitutions.rst

.. _paramtopics:

===========================
List of PARAMETERS by topic
===========================

This page lists parameters for ViPErLEED by grouping them
into categories. Not all parameters are listed. See the
:ref:`main PARAMETERS page<PARAMETERS>` for a complete
alphabetical list and other groupings.

.. note::
   While all parameters have a default value, parameters marked with a
   "**→**" usually require user input for simple but non-trivial
   systems.

|calc| execution
================

Defines what to calculate and where to start.

.. table::
  :width: 100%
  :widths: 30 70

  +----------------------------------+--------------------------------------------------------------------------+
  | Parameter                        | Function                                                                 |
  +==================================+==========================================================================+
  | :ref:`HALTING`                   | Sensitivity to things going wrong, i.e. how easily should ViPErLEED stop |
  +----------------------------------+--------------------------------------------------------------------------+
  | **→** :ref:`RUN`                 | Which parts of ViPErLEED / TensErLEED should be run, in order            |
  +----------------------------------+--------------------------------------------------------------------------+
  | :ref:`STOP`                      | Stop execution of ViPErLEED at next opportunity                          |
  +----------------------------------+--------------------------------------------------------------------------+
  | :ref:`TENSOR_INDEX`              | Which Tensor files to use for the delta calculation and search           |
  +----------------------------------+--------------------------------------------------------------------------+
  | :ref:`TENSOR_OUTPUT`             | Disable Tensor output for some or all layers                             |
  +----------------------------------+--------------------------------------------------------------------------+
  | :ref:`THEO_ENERGIES`             | What energy range to calculate                                           |
  +----------------------------------+--------------------------------------------------------------------------+
  | :ref:`TL_VERSION`                | Which version of TensErLEED to use                                       |
  +----------------------------------+--------------------------------------------------------------------------+
  | :ref:`TL_IGNORE_CHECKSUM`        | Skip calculation of TensErLEED source code checksums                     |
  +----------------------------------+--------------------------------------------------------------------------+

.. _input_structure_settings:

Input structure
===============

.. table::
  :width: 100%
  :widths: 30 70

  +----------------------------------+----------------------------------------------------------------------------------+
  | Parameter                        | Function                                                                         |
  +==================================+==================================================================================+
  | :ref:`BULK_LIKE_BELOW`           | Onset of unrelaxed slab, for automatic detection of minimal bulk and bulk repeat |
  +----------------------------------+----------------------------------------------------------------------------------+
  | :ref:`BULK_REPEAT`               | Thickness of the bulk repeat unit, or a bulk repeat vector                       |
  +----------------------------------+----------------------------------------------------------------------------------+
  | :ref:`LAYER_CUTS`                | How to separate the :ref:`POSCAR` file into layers                               |
  +----------------------------------+----------------------------------------------------------------------------------+
  | **→** :ref:`N_BULK_LAYERS`       | Define how many layers in the :ref:`POSCAR` file represent the bulk              |
  +----------------------------------+----------------------------------------------------------------------------------+
  | **→** :ref:`SUPERLATTICE`        | The relationship between the surface and bulk unit cells                         |
  +----------------------------------+----------------------------------------------------------------------------------+

Elements, vibration amplitudes and element concentrations
=========================================================

.. table::
  :width: 100%
  :widths: 30 70

  +----------------------------------+------------------------------------------------------------------------------------------------------+
  | Parameter                        | Function                                                                                             |
  +==================================+======================================================================================================+
  | :ref:`ELEMENT_MIX`               | Declare that sites in the :ref:`POSCAR` file can be occupied by different chemical elements          |
  +----------------------------------+------------------------------------------------------------------------------------------------------+
  | :ref:`ELEMENT_RENAME`            | Declare that an element in the POSCAR file is actually a different chemical element                  |
  +----------------------------------+------------------------------------------------------------------------------------------------------+
  | **→** :ref:`SITEDEF`             | Define which sites in the :ref:`POSCAR` file are special, i.e. have different vibration amplitude    |
  +----------------------------------+------------------------------------------------------------------------------------------------------+
  | :ref:`T_DEBYE`                   | Debye temperature of the system (only for automatically generating :ref:`VIBROCC`)                   |
  +----------------------------------+------------------------------------------------------------------------------------------------------+
  | :ref:`T_EXPERIMENT`              | Measurement temperature in experiment (only for automatically generating :ref:`VIBROCC`)             |
  +----------------------------------+------------------------------------------------------------------------------------------------------+
  | :ref:`VIBR_AMP_SCALE`            | Scaling factor, only for automatically generating :ref:`VIBROCC`                                     |
  +----------------------------------+------------------------------------------------------------------------------------------------------+

.. _symmetry_settings:

Symmetry determination and control
==================================

.. table::
  :width: 100%
  :widths: 30 70

  +----------------------------------+---------------------------------------------------------------------------------------------------+
  | Parameter                        | Function                                                                                          |
  +==================================+===================================================================================================+
  | :ref:`SYMMETRIZE_INPUT`          | Whether to move atoms in the :ref:`POSCAR` file to perfectly match the symmetry                   |
  +----------------------------------+---------------------------------------------------------------------------------------------------+
  | :ref:`SYMMETRYBULK`              | Manually set the symmetry to be used in beam averaging for the bulk, ignoring automatic detection |
  +----------------------------------+---------------------------------------------------------------------------------------------------+
  | :ref:`SYM_EPS`                   | Error tolerance during symmetry search                                                            |
  +----------------------------------+---------------------------------------------------------------------------------------------------+
  | :ref:`ISYM`                      | Manually set a symmetry, or turn symmetry off                                                     |
  +----------------------------------+---------------------------------------------------------------------------------------------------+
  | :ref:`SYMMETRY_FIND_ORI`         | Whether the symmetry search should look for the highest-symmetry point.                           |
  +----------------------------------+---------------------------------------------------------------------------------------------------+

Experimental setup
==================

.. table::
  :width: 100%
  :widths: 30 70

  +----------------------------------+---------------------------------------------------------------------------+
  | Parameter                        | Function                                                                  |
  +==================================+===========================================================================+
  | :ref:`AVERAGEBEAMS`              | Set beam averaging to assume an incidence other than :ref:`BEAMINCIDENCE` |
  +----------------------------------+---------------------------------------------------------------------------+
  | :ref:`BEAMINCIDENCE`             | Incidence angle and direction of the electron beam in experiment          |
  +----------------------------------+---------------------------------------------------------------------------+
  | :ref:`FILWF`                     | The LEED filament work function                                           |
  +----------------------------------+---------------------------------------------------------------------------+
  | :ref:`SCREEN_APERTURE`           | The aperture of the acceptance cone of the LEED screen                    |
  +----------------------------------+---------------------------------------------------------------------------+

Inner potential
===============

.. table::
  :width: 100%
  :widths: 30 70

  +----------------------------------+-----------------------------------------------------------+
  | Parameter                        | Function                                                  |
  +==================================+===========================================================+
  | :ref:`V0_IMAG`                   | Imaginary part of the inner potential                     |
  +----------------------------------+-----------------------------------------------------------+
  | :ref:`MUFTIN`                    | Real part of the inner potential                          |
  +----------------------------------+-----------------------------------------------------------+
  | :ref:`v0_z_onset`                | How far from the topmost atom the inner potential begins  |
  +----------------------------------+-----------------------------------------------------------+

Computational setup
===================

.. table::
  :width: 100%
  :widths: 30 70

  +----------------------------------+-----------------------------------------------------------+
  | Parameter                        | Function                                                  |
  +==================================+===========================================================+
  | :ref:`FORTRAN_COMP`              | Which Fortran compiler to use, and tags for compiling     |
  +----------------------------------+-----------------------------------------------------------+
  | **→** :ref:`NCORES`              | The number of CPUs to use                                 |
  +----------------------------------+-----------------------------------------------------------+

|R factor|
==========

.. table::
  :width: 100%
  :widths: 30 70

  +----------------------------------+-----------------------------------------------------------------------------------------------------+
  | Parameter                        | Function                                                                                            |
  +==================================+=====================================================================================================+
  | :ref:`IVSHIFTRANGE`              | Range and step size for shifting experimental and theoretical curves during |R-factor| calculation  |
  +----------------------------------+-----------------------------------------------------------------------------------------------------+
  | :ref:`RFACTORTYPE`               | Which definition of the |R factor| to use                                                           |
  +----------------------------------+-----------------------------------------------------------------------------------------------------+
  | :ref:`RFACTORLEGACY`             | Use legacy TensErLEED |R factor|                                                                    |
  +----------------------------------+-----------------------------------------------------------------------------------------------------+
  | :ref:`RFACTORSMOOTH`             | How strongly experimental beams are smoothed                                                        |
  +----------------------------------+-----------------------------------------------------------------------------------------------------+

.. _search_settings:

Search behavior
===============

.. table::
  :width: 100%
  :widths: 30 70

  +----------------------------------+-----------------------------------------------------------------------------------------------------+
  | Parameter                        | Function                                                                                            |
  +==================================+=====================================================================================================+
  | :ref:`MAX_TL_DISPLACEMENT`       | Whether to automatically re-do reference calculations when displacements get too large              |
  +----------------------------------+-----------------------------------------------------------------------------------------------------+
  | :ref:`SEARCHBEAMS`               | Whether to use |R factor| of integer, fractional, or all beams for the search                       |
  +----------------------------------+-----------------------------------------------------------------------------------------------------+
  | :ref:`SEARCH_CONVERGENCE`        | Convergence criteria for the search, and convergence-dependent parameter control                    |
  +----------------------------------+-----------------------------------------------------------------------------------------------------+
  | :ref:`SEARCH_CULL`               | Controls regular culling of worst-performing structures, and what to replace them with              |
  +----------------------------------+-----------------------------------------------------------------------------------------------------+
  | :ref:`SEARCHGENMAX`              | Maximum total number of generations that the search should run for                                  |
  +----------------------------------+-----------------------------------------------------------------------------------------------------+
  | :ref:`SEARCHPOP`                 | Number of trial structures used in parallel during the search                                       |
  +----------------------------------+-----------------------------------------------------------------------------------------------------+
  | :ref:`SEARCHSTART`               | How to initialize the search population                                                             |
  +----------------------------------+-----------------------------------------------------------------------------------------------------+
  | :ref:`OPTIMIZE`                  | Controls behavior of :ref:`full-dynamic optimization<Fdoptimization>` runs                          |
  +----------------------------------+-----------------------------------------------------------------------------------------------------+

Structural domains
==================

.. table::
  :width: 100%
  :widths: 30 70

  +----------------------------------+------------------------------------------------------------------------------------------------------------------------------------------+
  | Parameter                        | Function                                                                                                                                 |
  +==================================+==========================================================================================================================================+
  | :ref:`DOMAIN`                    | Define domains for :ref:`calculations involving multiple coexisting structural domains<domain_calculation>`                              |
  +----------------------------------+------------------------------------------------------------------------------------------------------------------------------------------+
  | :ref:`DOMAIN_STEP`               | Step width for structural domain coverage during search                                                                                  |
  +----------------------------------+------------------------------------------------------------------------------------------------------------------------------------------+
  | :ref:`SYMMETRY_CELL_TRANSFORM`   | Relationship between a supercell and the primitive surface unit cell (only relevant for :ref:`domain calculations<domain_calculation>`)  |
  +----------------------------------+------------------------------------------------------------------------------------------------------------------------------------------+

TensErLEED approximations
=========================

.. table::
  :width: 100%
  :widths: 30 70

  +----------------------------------+-------------------------------------------------------------------------------+
  | Parameter                        | Function                                                                      |
  +==================================+===============================================================================+
  | :ref:`ATTENUATION_EPS`           | Cutoff for beam propagation                                                   |
  +----------------------------------+-------------------------------------------------------------------------------+
  | :ref:`BULKDOUBLEEPS`             | Convergence criterion for bulk thickness in the TensErLEED calculation        |
  +----------------------------------+-------------------------------------------------------------------------------+
  | :ref:`BULKDOUBLEITER`            | Maximum bulk thickness in TensErLEED calculation                              |
  +----------------------------------+-------------------------------------------------------------------------------+
  | :ref:`LMAX`                      | Maximum angular momentum number; usually determined via :ref:`PHASESHIFTMIN`  |
  +----------------------------------+-------------------------------------------------------------------------------+
  | :ref:`PHASESHIFTMIN`             | Cutoff in phaseshifts magnitudes to determine :ref:`LMAX`                     |
  +----------------------------------+-------------------------------------------------------------------------------+

Debugging
=========

.. table::
  :width: 100%
  :widths: 30 70

  +----------------------------------+-------------------------------------------------------------------------------------------------+
  | Parameter                        | Function                                                                                        |
  +==================================+=================================================================================================+
  | :ref:`KEEP_REFCALC_DIRS`         | Toggle to keep the reference calculating execution directories                                  |
  +----------------------------------+-------------------------------------------------------------------------------------------------+
  | :ref:`LAYER_STACK_VERTICAL`      | How to choose layer stacking vectors in the TensErLEED input (debugging functionality only)     |
  +----------------------------------+-------------------------------------------------------------------------------------------------+
  | :ref:`LOG_LEVEL`                 | Set verbosity of the log file                                                                   |
  +----------------------------------+-------------------------------------------------------------------------------------------------+
  | :ref:`LOG_SEARCH`                | Output the search log file (may be very large, mostly for debugging)                            |
  +----------------------------------+-------------------------------------------------------------------------------------------------+
  | :ref:`SUPPRESS_EXE`              | Generate TensErLEED input files, but stop ViPErLEED before executing TensErLEED (for debugging) |
  +----------------------------------+-------------------------------------------------------------------------------------------------+

Output style
============

.. table::
  :width: 100%
  :widths: 30 70

  +----------------------------------+-------------------------------------------------+
  | Parameter                        | Function                                        |
  +==================================+=================================================+
  | :ref:`PLOT_IV`                   | Change appearance of the |R-factor| plot files  |
  +----------------------------------+-------------------------------------------------+


All PARAMETERS
==============

.. toctree::
    :maxdepth: 1
    :glob:

    parameters/*

