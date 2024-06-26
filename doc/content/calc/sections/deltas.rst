.. _sec_deltas:

============================
Delta-amplitudes calculation
============================

The delta-amplitudes calculation (called with :ref:`RUN = 2<run>`) is the
second part of a :ref:`tensor-LEED<tensor_leed>` calculation as implemented
in ViPErLEED. It requires the :ref:`tensor files<tensorszip>` produced by a
:ref:`reference calculation<ref-calc>`.

The delta calculation combines the :ref:`tensor files<tensorszip>`, i.e., the
tensor quantities :math:`T^\mathrm{ref}_{i,\mathbf{g};l,m;l',m'}` (see
:ref:`tensor_leed`), which only depend on the reference structure, with
the requested perturbations (see file :ref:`DISPLACEMENTS` for details)
to produce the corresponding changes in the scattered amplitudes of all
the beams :math:`\mathbf{g}`.

Delta-amplitudes calculation in ViPErLEED
-----------------------------------------

When a delta calculation is requested, ViPErLEED loads the
:ref:`tensor files<tensorszip>` and reads the :ref:`DISPLACEMENTS` file. Both 
are checked for consistency before proceeding. The :ref:`PHASESHIFTS` file —
needed to handle chemical substitutions — is also loaded, or generated
anew if required.

ViPErLEED then determines for which atoms delta amplitudes need to be
produced. If it is found that a previously generated delta file (stored
in the :file:`Deltas` directory) already contains some of the requested
amplitudes, their computation is skipped and a message is written in the
:file:`calc-<timestamp>.log` :ref:`log file<log_files>` accordingly.

As for the :ref:`reference calculation<ref-calc>`, ViPErLEED then moves the
needed TensErLEED source files to a temporary subfolder of the work directory
(\ :file:`Delta_Compile_<index>`) and compiles them **at run time**, storing
the :file:`SUPP/compile_logs/Delta_Compile_<index>.log`
:ref:`log file<log_files>`.
The input files (based on the :ref:`DISPLACEMENTS` file and the system 
symmetry) are written in the format expected by TensErLEED. Then
calculations for each atom are executed independently. Again, the 
parameter :ref:`ncores` determines how many processes run in parallel.

Once finished, ViPErLEED packs the resulting :ref:`delta files<deltaszip>`
into a ``.zip`` archive in the ``Deltas`` directory. Messages produced by
the :program:`delta.f` TensErLEED Fortran program are written to the
:file:`SUPP/delta-<timestamp>.log` :ref:`log file<log_files>`.
