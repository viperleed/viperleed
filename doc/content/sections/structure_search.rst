.. _sec_search:

================
Structure Search
================

The structure search (refered to as ``search`` in the code) is 
the third part of a :ref:`Tensor LEED<tensor_leed>` calculation as implemented 
in ViPErLEED (:ref:`RUN<run>` = 3).
It must follow a :ref:`refercence calculation<ref-calc>` and a 
:ref:`delta amplitude calculation<sec_deltas>` and requires stored 
:ref:`delta files<deltaszip>` to run.
In almost all cases it is the most computationally expensive part of the
calculation to run. The complexity is dependent on the size of the 
parameter space defined in the :ref:`DISPLACEMENTS file<displacements>`.

During the structure search, the employed optimization algorithm (see ref. :cite:p:`kottckeNewApproachAutomated1997`) samples surface structures in the configuration-space defined in :ref:`DISPLACEMENTS<displacements>`.
Diffraction intensities and a corresponding R-factor are calculated for every structure based on combinations of the pre-computed delta-amplitudes.
The optimzation tries to find the combination of parameters yielding the lowest possible R-factor.

The behaviour of the structure optimization algorithm is affected by multiple parameters (see :ref:`search behaviour<search_settings>`).
For details on the available options for geometrical, vibrational and occupational displacement vectors see the entry on the :ref:`DISPLACEMENTS file<displacements>`.
There are some caveats to the structure optimization in tensor LEED in general and the implementation in TensErLEED in particular.
See :ref:`the section on structure search in tensor leed<tensor_leed_search>` for details.


Structure Search in ViPErLEED
=============================

The structure search in ViPErLEED is based on the structure search section of TensErLEED, though as with the other sections, ViPErLEED takes care of handling the input and output for the legacy code.
The search is also the only part of ViPErLEED and TensErLEED that makes use of processes communicating via :term:`MPI`, if available (highly recommended).
When a structure search is executed in ViPErLEED the following main steps are performed before the actual calculation starts:

#.  The :ref:`file DISPLACEMENTS<displacements>` is read and interpreted.
#.  The current :ref:`delta files<deltaszip>` are loaded and checked for compatibilitiy.
#.  ViPErLEED generates and writes the TensErLEED inpute files ``rf.info``, ``PARAM`` and ``search.steu`` based on the slab and :ref:`EXPBEAMS file<expbeams>`.
#.  Based on the slab symmetry and the :ref:`symmetry settings<symmetry_settings>`, ViPErLEED determines the symmetry-linked parameters and writes the parameter-space input file control.chem.
#.  ViPErLEED will then, based on :ref:`N_CORES<ncores>` and the presece of ``mpirun`` fetch the corresponding TensErLEED source code files and compile them **at run-time**. **Note** that this will require the pre-compiled object files random_.o or MPIrandom_.o to be present. See the :ref:`getting started section<getting_started>` for details.
#.  The :ref:`search log file<log_files>` ``search-$timestamp`` is created and will be filled with progress information as the search continues.

With the preparation finished, the search is now executed (via ``mpirun`` if available).
Trial surface structures will be sampled using the algorithm described by :cite:t:`kottckeNewApproachAutomated1997`, with a starting configuration as defined by :ref:`SEARCH_START<searchstart>` and :ref:`SEARCH_POPULATION<searchpop>` parallel trial individuals.
ViPErLEED periodically monitors the search progress by reading the :ref:`SDTL` file and will report on the current best R-factor and the amount of sampled structures.
From this information, the files :ref:`search-progress.pdf<searchprogresspdf>` and :ref:`search-report.pdf<searchreportpdf>` will be generated and updated, which provides a graphical overview of the structure search progress and convergence.
If multiple structure optimizatios are run, the file :ref:`search-progress.pdf<searchprogresspdf>` will contain information related to the current run, whereas :ref:`search-report.pdf<searchreportpdf>` summarizes overall progress.

Once all required convergence criteria are met, the search will be cleanly aborted, the resulting files will be processed and :ref:`search-progress.pdf<searchprogresspdf>` and :ref:`search-report.pdf<searchreportpdf>` will be updated one last time with the final values.
After this, the structure search section is finished and ViPErLEED will continue with the next section as defined in the :ref:`RUN parameter<run>` (or stop if there are none).

.. warning:: 
  **Remember** to call the :ref:`bookkeeper utility<bookkeeper>` with the ``-c`` flag after a ViPErLEED run containing a strucutre search, if you want to continue from the found best-fit structure.
  **Otherwise the progress will be discarded** and following runs will start again from the refercence structure.
