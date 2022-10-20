.. _index:

ViPErLEED documentation
=======================

The Vienna Package for TensErLEED (ViPErLEED) is a collection of tools for modern LEED-I(V).

*TODO*: add nice introduction

.. include:: table_of_contents.rst

.. _work-segments:

ViPErLEED work segments
-----------------------

**Note: not sure where to put this, but here seemed fine for now**

ViPErLEED operates using a set of self-contained work segments (see also the :ref:`RUN<RUN>`  parameter). 
The three main segments, following the logic of calculations using the Tensor LEED approximation, are:

#. :ref:`Reference calculation<ref-calc>`: Full-dynamic LEED calculation, 
   which outputs a set of :ref:`theoretical beams<THEOBEAMS>` for a 
   given structure and :ref:`the "Tensors"<Tensorszip>` 
#. :ref:`Delta-Amplitudes generation<sec_deltas>`: Based on the Tensors and a :ref:`set of parameter variations specified by the user<DISPLACEMENTS>`, produces :ref:`"Delta files"<Deltaszip>`, specifying how these changes affect the beams
#. :ref:`Search<sec_search>`: Using the :ref:`Delta files<Deltaszip>`  to vary the theoretical beams, looks for a set of parameters such that the :ref:`R-factor<r-factor_calculation>` between the theoretical beams and :ref:`a given set of experimental beams<EXPBEAMS>`  is minimized.

Which of these segments should be executed must be specified using the :ref:`RUN<RUN>`  parameter, using the segment numbers in the list above. Besides these main three segments, there are also the following minor segments, which during normal ViPErLEED execution will be inserted automatically where appropriate:

-  :ref:`Initialization<initialization>`: Always runs at the beginning; reads and checks input files, runs symmetry search, generates derivative input files if appropriate.
-  :ref:`Superpos calculation<super_pos>`: Automatically runs after the search. Generates a set of theoretical beams based on the Tensor LEED approximation,
-  :ref:`R-factor calculation<r-factor_calculation>`: Automatically runs after the :ref:`reference calculation<ref-calc>`  and superpos segments, if an :ref:`experimental beams file<EXPBEAMS>` is present.
   Calculates the R-factor per beam and for the entire set of beams, and outputs an :ref:`Rfactor_plots pdf file<Rfactorplots>`.

Further specialized segments include:

-  :ref:`Error calculations<error_calculation>`: Based on a given reference structure (i.e. after a reference calculation has been run), calculate one-dimensional error curves for variation of a single parameter. Effectively, this calculates delta amplitudes for variations of a single parameter, and outputs the R-factor for every single configuration along that axis.
-  :ref:`Full-dynamic optimization<fdoptimization>`: Optimize parameters that cannot be varied during the search, like :ref:`BEAM_INCIDENCE<BEAMINCIDENCE>`, :ref:`V0_IMAG<INPOIM>`  or unit cell scaling. This is achieved by performing multiple full-dynamic (i.e. "reference") calculations (without Tensor output). Behaviour is controlled by the :ref:`OPTIMIZE<OPTIMIZE>`  parameter.

The pages listed above cover normal operation, in which the theoretical beams correspond to only one surface structure. If multiple structure coexist on the sample, the same segments need to be executed, but their behaviour is somewhat different, as described here:

-  :ref:`Domain calculations<domain_calculation>`: Reference calculations are run separately for the different domains (if necessary) and Delta-amplitudes are generated independently. The search then combines the optimization of the different structures for the best overall R-factor, compared to only one experimental beam set.


.. _list_input_files:

List of input files
-------------------

ViPErLEED has a number of input and control files, some of which are optional, depending on the work-segment (:ref:`see above<work-segments>`).

+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| File                                 | Function                                          | Required                                                          |
+======================================+===================================================+===================================================================+
| :ref:`PARAMETERS<PARAMETERS>`        | What to execute, and how to interpret other input | Always required                                                   |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`POSCAR<POSCAR>`                | Structural information                            | Always required                                                   |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`VIBROCC<VIBOCCIN>`             | Vibrational amplitudes and site occupations       | Always required, can be generated automatically                   |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`IVBEAMS<IVBEAMS>`              | Which beams to calculate                          | Always required, can be generated automatically from EXPBEAMS.csv |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`DISPLACEMENTS<DISPLACEMENTS>`  | What to vary during the search                    | Required for delta calculations and search                        |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`EXPBEAMS.csv<EXPBEAMS>`        | Experimental I(V) curves                          | Required for R-factor calculations and search                     |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`PHASESHIFTS<PHASESHIFTS>`      | Phase shifts per energy \* site \* L              | Generated automatically if needed                                 |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+

List of output files
--------------------

Output files, stored in the OUT subfolder
'''''''''''''''''''''''''''''''''''''''''

+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| File                                                                 | Function                                                         | Output by                                         |
+======================================================================+==================================================================+===================================================+
| :ref:`POSCAR_OUT<POSCAR>`                                            | Best-fit structure                                               | Search                                            |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`VIBROCC_OUT<VIBOCCIN>`                                         | Best-fit vibrational amplitudes and occupations                  | Search                                            |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`THEOBEAMS.csv, THEOBEAMS.pdf<THEOBEAMS>`                       | Theoretical I(V) curves                                          | Reference calculation                             |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`FITBEAMS.csv<FITBEAMS>`                                        | Theoretical I(V) curves                                          | Superpos calculation                              |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`Search-progress.pdf<searchprogresspdf>`                        | Progress information on current search                           | Search                                            |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`Search-report.pdf<searchreportpdf>`                            | Summary of several consecutive searches                          | Search                                            |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`Rfactor_plots.pdf<rfactorplots>`                               | Plots of experimental and theoretical I(V) curves with R-factors | R-factor calculation                              |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`Rfactor_analysis.pdf<rfactoranalysis>`                         | Same as Rfactor_plots, with Y-functions and analysis             | R-factor calculation                              |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`Errors.csv, Errors.pdf<Errorspdf>`                             | R-factors for one-dimensional parameter variation                | :ref:`Error calculation<error_calculation>`       |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`FD_Optimization.csv, FD_Optimization.pdf<fdoptimizationdata>`  | R-factors for variation in full-dynamic optimization             | :ref:`Full-dynamic optimization<Fdoptimization>`  |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`FD_Optimization_beams.pdf<fdoptimizationbeams>`                | I(V) curves generated during full-dynamic optimization           | :ref:`Full-dynamic optimization<Fdoptimization>`  |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`PatternInfo.tlm<PatternInfo>`                                  | Input for :term:`GUI`                                            | Initialization                                    |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`Tensors.zip<tensorszip>`                                       | Tensor files for Delta-amplitudes calculation                    | Reference calculation                             |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`Deltas.zip<Deltaszip>`                                         | Delta files for search                                           | Delta-amplitudes calculation                      |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| Raw output files from TensErLEED:                                    |                                                                  |                                                   |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`R_OUT<r_out>`                                                  | R-factors for every beam and every inner potential shift         | R-factor calculation                              |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`control.chem<controlchem>`                                     | Search parameter values for the latest generation                | Search                                            |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`SD.TL<sdtl>`                                                   | Search parameter values and R-factors per generation             | Search                                            |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`refcalc-fd.out<fd_out>`                                        | Theoretical I(V) curves                                          | Reference calculation                             |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`superpos-spec.out<superpos-spec_out>`                          | Theoretical I(V) curves                                          | Superpos calculation                              |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+

--------------

Supplementary files, stored in the SUPP subfolder
'''''''''''''''''''''''''''''''''''''''''''''''''

+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| File                                                                                         | Function                                                                                  |
+==============================================================================================+===========================================================================================+
| Input for TensErLEED                                                                         |                                                                                           |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`BEAMLIST<BEAMLIST>`                                                                    | All relevant beams                                                                        |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`AUXBEAMS<AUXBEAMS>`                                                                    | Information on beams used by multiple TensErLEED parts, based on :ref:`IVBEAMS<IVBEAMS>`  |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`AUXEXPBEAMS<AUXEXPBEAMS>`                                                              | Contents of :ref:`EXPBEAMS.csv<EXPBEAMS>`, formatted for TensErLEED                       |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`refcalc-PARAM, muftin.f, AUXGEO, AUXLATGEO, AUXNONSTRUCT, refcalc-FIN<refcalc-input>`  | Input files for the reference calculation                                                 |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`delta-input<delta-input>`                                                              | Input files for the delta-amplitudes calculations                                         |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`search.steu, search-rf.info, restrict.f, search-PARAM<search-input>`                   | Input files for the search                                                                |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`superpos-PARAM, superpos-CONTRIN<superpos-input>`                                      | Input files for the superpos calculation                                                  |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`rfactor-PARAM, rfactor-WEXPEL<rfactor-input>`                                          | Input files for the R-factor calculation                                                  |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| Files for user information and troubleshooting                                               |                                                                                           |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`POSCAR_bulk, POSCAR_bulk_appended, POSCAR_oricell<POSCAR>`                             | Modified POSCAR files to check structure                                                  |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`searchpars.info<searchparsinfo>`                                                       | Information on search parameters used in the search                                       |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`Phaseshifts-plots.pdf<phaseshiftplots>`                                                | Plots of PHASESHIFTS, generated during initialization                                     |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+

.. _utilities:

Utilities
---------

ViPErLEED includes several small utilities for file processing and organization. Proper embedding of these into a shared interface is pending; the following list is a loose summary and may not be complete.

-  :ref:`Bookkeeper<bookkeeper>`: Helper utility for the TensErLEED Manager tleedm. Sorts files from previous runs into a "history" folder, and keeps track in the history.info file.
-  Combine-POSCAR, POSCAR_get_bulk_repeat: Very simple POSCAR editing scripts. **TODO** documentation
-  ModifyPhaseshifts: Simple utility for taking an existing PHASESHIFTS file and duplicating and/or re-arranging blocks. **TODO** usage (how to call it)
-  AUXEXPBEAMS_TO_EXPBEAMS: Transforms :ref:`AUXEXPBEAMS<AUXEXPBEAMS>`  format files (input for TensErLEED) to :ref:`EXPBEAMS.csv<EXPBEAMS>`  format. **TODO** usage (how to call it)
