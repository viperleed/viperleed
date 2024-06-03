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

+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| File                                 | Function                                          | Required                                                          |
+======================================+===================================================+===================================================================+
| :ref:`PARAMETERS<PARAMETERS>`        | What to execute, and how to interpret other input | Always required                                                   |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`POSCAR<POSCAR>`                | Structural information                            | Always required                                                   |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`VIBROCC<vibrocc>`              | Vibrational amplitudes and site occupations       | Always required, can be generated automatically                   |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`IVBEAMS<IVBEAMS>`              | Which beams to calculate                          | Always required, can be generated automatically from EXPBEAMS.csv |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`DISPLACEMENTS<DISPLACEMENTS>`  | What to vary during the search                    | Required for delta calculations and search                        |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`EXPBEAMS.csv<EXPBEAMS>`        | Experimental |IV| curves                          | Required for R-factor calculations and search                     |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`PHASESHIFTS<PHASESHIFTS>`      | Contains elastic electron scattering phaseshifts  | Generated automatically if needed                                 |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+


.. _output_files:

Output files
------------

ViPErLEED produces a number of output files, containing the results of 
the requested calculations. They are stored in the ``OUT`` subfolder.

+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| File                                                                 | Function                                                         | Output by                                         |
+======================================================================+==================================================================+===================================================+
| :ref:`POSCAR_OUT<POSCAR>`                                            | Best-fit structure                                               | Search                                            |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`VIBROCC_OUT<vibrocc>`                                          | Best-fit vibrational amplitudes and occupations                  | Search                                            |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`THEOBEAMS.csv, THEOBEAMS.pdf<THEOBEAMS>`                       | Theoretical |IV| curves                                          | Reference calculation                             |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`FITBEAMS.csv<FITBEAMS>`                                        | Theoretical |IV| curves                                          | Superpos calculation                              |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`Search-progress.pdf<searchprogresspdf>`                        | Progress information on current search                           | Search                                            |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`Search-report.pdf<searchreportpdf>`                            | Summary of several consecutive searches                          | Search                                            |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`Rfactor_plots.pdf<rfactorplots>`                               | Plots of experimental and theoretical |IV| curves with R factors | R-factor calculation                              |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`Rfactor_analysis.pdf<rfactoranalysis>`                         | Same as Rfactor_plots, with Y-functions and analysis             | R-factor calculation                              |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`Errors.csv, Errors.pdf<Errorspdf>`                             | R-factors for one-dimensional parameter variation                | :ref:`Error calculation<error_calculation>`       |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`FD_Optimization.csv, FD_Optimization.pdf<fdoptimizationdata>`  | R-factors for variation in full-dynamic optimization             | :ref:`Full-dynamic optimization<Fdoptimization>`  |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`FD_Optimization_beams.pdf<fdoptimizationbeams>`                | |IV| curves generated during full-dynamic optimization           | :ref:`Full-dynamic optimization<Fdoptimization>`  |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`PatternInfo.tlm<PatternInfo>`                                  | Input for :term:`GUI`                                            | Initialization                                    |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`Tensors.zip<tensorszip>`                                       | Tensor files for Delta-amplitudes calculation                    | Reference calculation                             |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`Deltas.zip<Deltaszip>`                                         | Delta files for search                                           | Delta-amplitudes calculation                      |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| **Raw output files from TensErLEED:**                                                                                                                                                       |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`R_OUT<r_out>`                                                  | R-factors for every beam and every inner potential shift         | R-factor calculation                              |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`control.chem<controlchem>`                                     | Search parameter values for the latest generation                | Search                                            |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`SD.TL<sdtl>`                                                   | Search parameter values and R-factors per generation             | Search                                            |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`refcalc-fd.out<fd_out>`                                        | Theoretical |IV| curves                                          | Reference calculation                             |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :ref:`superpos-spec.out<superpos-spec_out>`                          | Theoretical |IV| curves                                          | Superpos calculation                              |
+----------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+


.. _supp_files:

Supplementary files
-------------------

ViPErLEED produces supplementary files that are required during execution, that 
contain intermediate results or that may be of interest for debugging purposes.
These files are stored in the ``SUPP`` subfolder of the ``work`` directory.

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
| :ref:`POSCAR_bulk, POSCAR_bulk_appended, POSCAR_oricell, POSCAR_vacuum_corrected<POSCAR>`    | Modified POSCAR files to check structure                                                  |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`searchpars.info<searchparsinfo>`                                                       | Information on search parameters used in the search                                       |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+
| :ref:`Phaseshifts-plots.pdf<phaseshiftplots>`                                                | Plots of PHASESHIFTS, generated during initialization                                     |
+----------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------+

.. toctree::
   :glob:

   files/input/*


.. toctree::

   files/output/theobeams
   files/output/fitbeams
   files/output/rfactorplots
   files/output/rfactoranalysis
   files/output/searchreportpdf
   files/output/searchprogresspdf
   files/output/errorspdf
   files/output/fdoptimizationdata
   files/output/fdoptimizationbeams
   files/output/log_files
   files/output/tensorszip
   files/output/deltaszip

.. toctree::
   :glob:

   files/supp/*
   files/raw_tenserleed/*