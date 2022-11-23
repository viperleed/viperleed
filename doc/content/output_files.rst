.. _output_files:

List of output files
--------------------

ViPErLEED produces a number of output files, containing the results of 
the requested calculations. They are stored in the ``OUT`` subfolder.

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
| **Raw output files from TensErLEED:**                                                                                                                                                       |
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

