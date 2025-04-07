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

    +-------------------------------+----------------------------+----------------------+
    | File                          | Function                   | Required             |
    +===============================+============================+======================+
    | :ref:`PARAMETERS`             | What to execute, and how \ | Always required      |
    |                               | to interpret other input   |                      |
    +-------------------------------+----------------------------+----------------------+
    | :ref:`POSCAR`                 | Structural information     | Always required      |
    +-------------------------------+----------------------------+----------------------+
    | :ref:`VIBROCC`                | Vibrational amplitudes \   | Always required, \   |
    |                               | and site occupations       | can be generated \   |
    |                               |                            | automatically        |
    +-------------------------------+----------------------------+----------------------+
    | :ref:`IVBEAMS`                | Which beams to calculate   | Always required, \   |
    |                               |                            | can be generated \   |
    |                               |                            | automatically        |
    |                               |                            | from EXPBEAMS.csv    |
    +-------------------------------+----------------------------+----------------------+
    | :ref:`DISPLACEMENTS`          | What to vary during \      | Required for delta \ |
    |                               | the search                 | calculations and \   |
    |                               |                            | search               |
    +-------------------------------+----------------------------+----------------------+
    | :ref:`EXPBEAMS.csv<EXPBEAMS>` | Experimental |IV| curves   | Required for \       |
    |                               |                            | |R-factor| \         |
    |                               |                            | calculation and \    |
    |                               |                            | search               |
    +-------------------------------+----------------------------+----------------------+
    | :ref:`PHASESHIFTS`            | Contains elastic \         | Generated \          |
    |                               | electron scattering \      | automatically \      |
    |                               | phase shifts               | if needed            |
    +-------------------------------+----------------------------+----------------------+

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

    +-------------------------------+---------------------------------+------------------------------+
    | File                          | Function                        | Output by                    |
    +===============================+=================================+==============================+
    | :ref:`PARAMETERS`             | Automatically edited version \  | :ref:`Initialization`        |
    |                               | of the user-given PARAMETERS    |                              |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`POSCAR`                 | Modified input POSCAR \         | :ref:`Initialization`, \     |
    |                               | (\ :ref:`Initialization`) \     | :ref:`Search<sec_search>`, \ |
    |                               | or best-fit structure           | :ref:`Fdoptimization`        |
    |                               | (\ :ref:`Search<sec_search>`, \ |                              |
    |                               | :ref:`Fdoptimization`)          |                              |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`VIBROCC`                | Auto-generated VIBROCC \        | :ref:`Initialization`, \     |
    |                               | (\ :ref:`Initialization`) \     | :ref:`Search<sec_search>`    |
    |                               | or best-fit vibration \         |                              |
    |                               | amplitudes and occupations \    |                              |
    |                               | (\ :ref:`Search<sec_search>`)   |                              |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`THEOBEAMS`.csv/.pdf     | Theoretical |IV| curves         | :ref:`ref-calc`              |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`FITBEAMS.csv<FITBEAMS>` | Theoretical |IV| curves         | :ref:`super_pos`             |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`searchprogresspdf`      | Progress information \          | :ref:`Search<sec_search>`    |
    |                               | on current search               |                              |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`searchreportpdf`        | Summary of several \            | :ref:`Search<sec_search>`    |
    |                               | consecutive searches            |                              |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`rfactorplots`           | Plots of experimental and \     | |R-factor| calculation       |
    |                               | theoretical |IV| curves with \  |                              |
    |                               | |R factor|\ s                   |                              |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`rfactoranalysis`        | Same as :ref:`rfactorplots`, \  | |R-factor| calculation       |
    |                               | with Y-functions and analysis   |                              |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`errorspdf_header`       | |R factor|\ s for \             | :ref:`error_calculation`     |
    |                               | one-dimensional \               |                              |
    |                               | parameter variation             |                              |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`fdoptimizationdata`     | |R factor|\ s for variation \   | :ref:`Fdoptimization`        |
    |                               | in full-dynamic optimization    |                              |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`fdoptimizationbeams`    | |IV| curves generated during \  | :ref:`Fdoptimization`        |
    |                               | full-dynamic optimization       |                              |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`PatternInfo`            | Input for :term:`GUI`           | :ref:`Initialization`        |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`Tensors.zip<tensorszip>`| Tensor files for \              | :ref:`ref-calc`              |
    |                               | delta-amplitude calculation     |                              |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`Deltas.zip<Deltaszip>`  | Delta files for search          | :ref:`sec_deltas`            |
    +-------------------------------+---------------------------------+------------------------------+
    | **Raw output files** \                                                                         |
    | **from TensErLEED:**                                                                           |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`R_OUT`                  | |R factor|\ s for every \       | |R-factor| calculation       |
    |                               | beam and every \                |                              |
    |                               | inner-potential shift           |                              |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`controlchem`            | Search-parameter values for \   | :ref:`Search<sec_search>`    |
    |                               | the latest generation           |                              |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`sdtl`                   | Search parameter values and \   | :ref:`Search<sec_search>`    |
    |                               | |R factor|\ s per generation    |                              |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`fd_out`                 | Theoretical |IV| curves         | :ref:`ref-calc`              |
    +-------------------------------+---------------------------------+------------------------------+
    | :ref:`superpos-spec_out`      | Theoretical |IV| curves         | :ref:`super_pos`             |
    +-------------------------------+---------------------------------+------------------------------+

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

+------------------------------------------+-----------------------------------+
| File                                     | Function                          |
+==========================================+===================================+
| **Input for TensErLEED**                 |                                   |
+------------------------------------------+-----------------------------------+
| :ref:`BEAMLIST`                          | All relevant beams                |
+------------------------------------------+-----------------------------------+
| :ref:`AUXBEAMS`                          | Information on beams used by \    |
|                                          | multiple TensErLEED parts, \      |
|                                          | based on :ref:`IVBEAMS<IVBEAMS>`  |
+------------------------------------------+-----------------------------------+
| :ref:`AUXEXPBEAMS`                       | Contents of \                     |
|                                          | :ref:`EXPBEAMS.csv<EXPBEAMS>`, \  |
|                                          | formatted for TensErLEED          |
+------------------------------------------+-----------------------------------+
| :ref:`refcalc-PARAM<refcalc-input>`, \   | Input files for the \             |
| :ref:`muftin.f<refcalc-input>`, \        | reference calculation             |
| :ref:`AUXGEO<refcalc-input>`, \          |                                   |
| :ref:`AUXLATGEO<refcalc-input>`, \       |                                   |
| :ref:`AUXNONSTRUCT<refcalc-input>`, \    |                                   |
| :ref:`refcalc-FIN<refcalc-input>`        |                                   |
+------------------------------------------+-----------------------------------+
| :ref:`delta-input`                       | Input files for the \             |
|                                          | delta-amplitudes calculation      |
+------------------------------------------+-----------------------------------+
| :ref:`search.steu<search-input>`, \      | Input files for the search        |
| :ref:`search-rf.info<search-input>`, \   |                                   |
| :ref:`restrict.f<search-input>`, \       |                                   |
| :ref:`rsearch-PARAM<search-input>`       |                                   |
+------------------------------------------+-----------------------------------+
| :ref:`superpos-PARAM<superpos-input>`, \ | Input files for the \             |
| :ref:`superpos-CONTRIN<superpos-input>`  | superpos calculation              |
+------------------------------------------+-----------------------------------+
| :ref:`rfactor-PARAM<rfactor-input>`, \   | Input files for the \             |
| :ref:`rfactor-WEXPEL<rfactor-input>`     | |R-factor| calculation            |
+------------------------------------------+-----------------------------------+
| **Files for user information** \         |                                   |
| **and troubleshooting**                  |                                   |
+------------------------------------------+-----------------------------------+
| :ref:`POSCAR_bulk<POSCAR>`, \            | Modified POSCAR files to \        |
| :ref:`POSCAR_bulk_appended<POSCAR>`, \   | check structure                   |
| :ref:`POSCAR_oricell<POSCAR>`, \         |                                   |
| :ref:`POSCAR_vacuum_corrected<POSCAR>`   |                                   |
+------------------------------------------+-----------------------------------+
| :ref:`searchpars.info<searchparsinfo>`   | Information on search \           |
|                                          | parameters used in the search     |
+------------------------------------------+-----------------------------------+
| :ref:`phaseshiftplots`                   | Plots of PHASESHIFTS, generated \ |
|                                          | during initialization             |
+------------------------------------------+-----------------------------------+


.. todo::
    "manifest" is missing from the above table, but present in the
    toctree below.
    
    delta-input names a generic "PARAM" but I guess this is named differently
    in SUPP. Perhaps would be better to list the files like the others?

.. toctree::
    :glob:

    files/supp/*
