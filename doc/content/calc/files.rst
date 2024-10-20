.. include:: /substitutions.rst

.. _files:

=====
Files
=====


.. _list_input_files:

Input files
-----------

ViPErLEED has a number of input and control files, some of which are optional,
depending on the work-segment (:ref:`see above<work-segments>`).

.. table::
    :width: 100%
    :widths: 22 39 39

    +--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
    | File                                 | Function                                          | Required                                                          |
    +======================================+===================================================+===================================================================+
    | :ref:`PARAMETERS`                    | What to execute, and how to interpret other input | Always required                                                   |
    +--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
    | :ref:`POSCAR`                        | Structural information                            | Always required                                                   |
    +--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
    | :ref:`VIBROCC`                       | Vibrational amplitudes and site occupations       | Always required, can be generated automatically                   |
    +--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
    | :ref:`IVBEAMS`                       | Which beams to calculate                          | Always required, can be generated automatically from EXPBEAMS.csv |
    +--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
    | :ref:`DISPLACEMENTS`                 | What to vary during the search                    | Required for delta calculations and search                        |
    +--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
    | :ref:`EXPBEAMS.csv<EXPBEAMS>`        | Experimental |IV| curves                          | Required for |R-factor| calculations and search                   |
    +--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
    | :ref:`PHASESHIFTS`                   | Contains elastic electron scattering phaseshifts  | Generated automatically if needed                                 |
    +--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+

.. toctree::
   :hidden:
   :glob:

   files/input/*


.. _output_files:

Output files
------------

ViPErLEED produces a number of output files, containing the results of
the requested calculations. They are stored in the ``OUT`` subfolder.


.. table::
    :width: 100%
    :widths: 32 42 26

    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | File                                            | Function                                                               | Output by                     |
    +=================================================+========================================================================+===============================+
    | :ref:`POSCAR_OUT<POSCAR>`                       | Best-fit structure                                                     | :ref:`Search<sec_search>`     |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`VIBROCC_OUT<vibrocc>`                     | Best-fit vibrational amplitudes and occupations                        | :ref:`Search<sec_search>`     |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`THEOBEAMS.csv, THEOBEAMS.pdf<THEOBEAMS>`  | Theoretical |IV| curves                                                | :ref:`ref-calc`               |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`FITBEAMS.csv<FITBEAMS>`                   | Theoretical |IV| curves                                                | :ref:`super_pos`              |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`searchprogresspdf`                        | Progress information on current search                                 | :ref:`Search<sec_search>`     |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`searchreportpdf`                          | Summary of several consecutive searches                                | :ref:`Search<sec_search>`     |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`rfactorplots`                             | Plots of experimental and theoretical |IV| curves with |R factor|\ s   | |R-factor| calculation        |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`rfactoranalysis`                          | Same as :ref:`rfactorplots`, with Y-functions and analysis             | |R-factor| calculation        |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`errorspdf_header`                         | |R factor|\ s for one-dimensional parameter variation                  | :ref:`error_calculation`      |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`fdoptimizationdata`                       | |R factor|\ s for variation in full-dynamic optimization               | :ref:`Fdoptimization`         |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`fdoptimizationbeams`                      | |IV| curves generated during full-dynamic optimization                 | :ref:`Fdoptimization`         |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`PatternInfo`                              | Input for :term:`GUI`                                                  | :ref:`Initialization`         |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`Tensors.zip<tensorszip>`                  | Tensor files for Delta-amplitudes calculation                          | :ref:`ref-calc`               |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`Deltas.zip<Deltaszip>`                    | Delta files for search                                                 | :ref:`sec_deltas`             |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | **Raw output files from TensErLEED:**                                                                                                                    |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`R_OUT`                                    | |R factor|\ s for every beam and every inner potential shift           | |R-factor| calculation        |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`controlchem`                              | Search parameter values for the latest generation                      | :ref:`Search<sec_search>`     |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`sdtl`                                     | Search parameter values and |R factor|\ s per generation               | :ref:`Search<sec_search>`     |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`fd_out`                                   | Theoretical |IV| curves                                                | :ref:`ref-calc`               |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+
    | :ref:`superpos-spec_out`                        | Theoretical |IV| curves                                                | :ref:`super_pos`              |
    +-------------------------------------------------+------------------------------------------------------------------------+-------------------------------+

.. toctree::
    :hidden:
    :glob:

    files/output/*
    files/raw_tenserleed/*


.. _supp_files:

Supplementary files
-------------------

ViPErLEED produces supplementary files that are required during execution, that
contain intermediate results or that may be of interest for debugging purposes.
These files are stored in the ``SUPP`` subfolder of the root directory.

+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| File                                                                                         | Function                                                                                  |
+==============================================================================================+===========================================================================================+
| **Input for TensErLEED**                                                                     |                                                                                           |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`BEAMLIST`                                                                              | All relevant beams                                                                        |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`AUXBEAMS`                                                                              | Information on beams used by multiple TensErLEED parts, based on :ref:`IVBEAMS<IVBEAMS>`  |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`AUXEXPBEAMS`                                                                           | Contents of :ref:`EXPBEAMS.csv<EXPBEAMS>`, formatted for TensErLEED                       |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`refcalc-PARAM, muftin.f, AUXGEO, AUXLATGEO, AUXNONSTRUCT, refcalc-FIN<refcalc-input>`  | Input files for the reference calculation                                                 |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`delta-input`                                                                           | Input files for the delta-amplitudes calculations                                         |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`search.steu, search-rf.info, restrict.f, search-PARAM<search-input>`                   | Input files for the search                                                                |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`superpos-PARAM, superpos-CONTRIN<superpos-input>`                                      | Input files for the superpos calculation                                                  |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`rfactor-PARAM, rfactor-WEXPEL<rfactor-input>`                                          | Input files for the |R-factor| calculation                                                |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| **Files for user information and troubleshooting**                                           |                                                                                           |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`POSCAR_bulk, POSCAR_bulk_appended, POSCAR_oricell, POSCAR_vacuum_corrected<POSCAR>`    | Modified POSCAR files to check structure                                                  |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`searchpars.info<searchparsinfo>`                                                       | Information on search parameters used in the search                                       |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`phaseshiftplots`                                                                       | Plots of PHASESHIFTS, generated during initialization                                     |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+


.. todo::
    "manifest" is missing from the above table, but present in the
    toctree below.
    
    delta-input names a generic "PARAM" but I guess this is named differently
    in SUPP. Perhaps would be better to list the files like the others?

.. toctree::
    :glob:

    files/supp/*
