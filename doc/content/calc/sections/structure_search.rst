.. include:: /substitutions.rst

.. _sec_search:

================
Structure search
================

The structure search (referred to as ``search`` in the code) is
the third part of a :ref:`tensor-LEED<tensor_leed>` calculation as
implemented in ViPErLEED (:ref:`RUN = 3<run>`).
It must follow a :ref:`reference calculation<ref-calc>` and a
:ref:`delta-amplitudes calculation<sec_deltas>`. It requires stored
:ref:`delta files<deltaszip>` to run. In almost all cases, it is the
most computationally expensive part of a |LEED-IV| calculation. The
complexity depends on the size of the parameter space defined in the
:ref:`DISPLACEMENTS` file.

During the structure search, the optimization algorithm samples surface
structures in the configuration space defined in :ref:`DISPLACEMENTS`.
Diffraction intensities and the corresponding |R factor| are calculated
for these structures based on combinations of the precomputed delta
amplitudes. The optimization tries to find the combination of parameters
yielding the smallest possible |R factor|.

The behavior of the structure-optimization algorithm is affected by multiple
parameters (see :ref:`search behaviour<search_settings>`). See the section
describing the :ref:`DISPLACEMENTS` file for details on the available options
for geometric, vibration, and occupation displacements. There are some
caveats to the structure optimization in tensor LEED in general, and to
the implementation in TensErLEED in particular.
See :ref:`the section on structure search in tensor LEED<tensor_leed_search>`
for details.


.. _optimization_algorithm:

Optimization algorithm
======================

The rough |R-factor| surface :cite:p:`rousTensorLEEDApproximation1992` greatly
limits the pool of applicable optimization algorithms. TensErLEED employs
a modified random-sampling strategy with a down-step criterion, as described
by :cite:t:`kottckeNewApproachAutomated1997`. The optimization is performed
in parallel for a set of individuals (i.e., independent parameter
combinations), defined in ViPErLEED by the :ref:`searchpop` parameter. The
initial configurations for the optimization are defined by :ref:`searchstart`.

For each search step (called "generation", based on the terminology of genetic
algorithms), a new grid point in the parameter space is selected *randomly*,
but based on a probability distribution centered on the current position.
The |R factor| is calculated for the selected parameter combination and the new
parameter set is accepted **only if** the |R factor| for the new configuration
is smaller than for the previous configuration.
The width of the probability distribution is determined by the current
|R factor|, the number of displacements, and the ``gaussian`` flag of
the :ref:`SEARCH_CONVERGENCE` parameter.

ViPErLEED enables more sophisticated control over the search process than
is possible with TensErLEED alone. Different types of convergence criteria
and an automatic scaling of the probability distribution can be set via
:ref:`SEARCH_CONVERGENCE`. Furthermore, whenever
:ref:`partial convergence<search_convergence>` is reached, a portion of the
search population can be dropped and reinitialized to escape from local minima.
By default, the search population is partially reinitialized using a custom
genetic algorithm (see :ref:`SEARCH_CULL` for details).



Structure search in ViPErLEED
=============================

The structure search in ViPErLEED is based on the corresponding section of
TensErLEED. As with the other sections, ViPErLEED takes care of handling the
input and output for the legacy code. The search is also the only part of
ViPErLEED and TensErLEED that makes use of processes communicating via
:term:`MPI`, if available (highly recommended).

When a structure search is executed in ViPErLEED, the following main steps
are performed before the actual calculation starts:

#.  The :ref:`displacements` file is read and interpreted.
#.  The current :ref:`delta files<deltaszip>` are loaded and checked for
    compatibility.
#.  The TensErLEED input files :file:`rf.info`, :file:`PARAM`, and
    :file:`search.steu` are prepared based on the slab and on the
    :ref:`EXPBEAMS file<expbeams>`. (They are available for inspection among
    the :ref:`supplementary files<supp_files>` in the :file:`SUPP` directory.)
#.  Based on the slab symmetry and the
    :ref:`symmetry settings<symmetry_settings>`, symmetry-linked parameters
    are identified. The :file:`control.chem` input file for TensErLEED is
    written.
#.  Based on the :ref:`NCORES` parameter and on the presence of
    :program:`mpirun` on the machine, the correct TensErLEED source files
    are collected. They are then compiled **at run time** (compilation
    information saved to :file:`SUPP/compile_logs/compile-search.log`).
    When using early TensErLEED versions (< v1.7.4) the precompiled object
    files :file:`random_.o` or :file:`MPIrandom_.o` must be available at
    this stage. More information on compiling these files can be found in
    :ref:`this<mpirandom>` part of the :ref:`installation` section.
#.  The :ref:`search log file<log_files>` :file:`search-<timestamp>.log` is
    created. It will be filled with progress information as the search
    continues.

With the preparation finished, the TensErLEED search program
is started (using :program:`mpirun`, if available).

ViPErLEED periodically monitors the progress of the search, reporting the best
|R factor| achieved up to that point and the number of sampled configurations.
Based on the information read from the :ref:`sdtl` file produced by TensErLEED,
the :ref:`searchprogresspdf` and :ref:`searchreportpdf`  files are generated
and periodically updated. They provide a graphical overview of the progress of
the structure search and of its convergence. :ref:`searchprogresspdf` contains
information related exclusively to the current TensErLEED structure
optimization, i.e., one block in the :ref:`displacements` file. ViPErLEED
enables the user to chain or loop multiple TensErLEED structure optimizations
(see the :ref:`displacements` syntax for details). In that case,
:ref:`searchreportpdf` summarizes the overall progress, including
all optimization runs.

Once all the convergence criteria (see :ref:`SEARCH_CONVERGENCE`) are
met, the search is cleanly aborted, the results are processed, and the
:ref:`searchprogresspdf` and :ref:`searchreportpdf` files are updated
one last time with the final values. This concludes the structure-search
section. ViPErLEED proceeds by first executing a :ref:`super_pos`, then
continuing to the next segment as defined in the  :ref:`RUN` parameter
(or it stops if there is none).

.. warning::
  **Remember** to call the :ref:`bookkeeper utility<bookkeeper>` with
  the ``-c`` flag after a ViPErLEED run containing a structure search,
  if you want to continue from the found best-fit structure.
  **Otherwise the progress will be discarded** and following runs will
  start again from the reference structure, unless :ref:`POSCAR` and
  :ref:`VIBROCC` are manually copied from the ``OUT`` directory.

