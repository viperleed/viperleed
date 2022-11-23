.. _sec_deltas:

=================
Delta Calculation
=================

The delta amplitude calculation (called with :ref:`RUN = 2<run>`) is the 
second part of a :ref:`Tensor LEED<tensor_leed>` calculation as implemented 
in ViPErLEED.
It requires a :ref:`refercence calculation<ref-calc>` to be completed 
beforehand and create the :ref:`tensor files<tensorszip>` which are stored.

In the delta calculation, these :ref:`tensor files<tensorszip>`, i.e. the
tensor quantities :math:`T^{ref}_{i;l,m;l',m'}` (see :ref:`tensor leed<tensor_leed>`),
which only depend on the refercence structure, are combined with the requested 
pertubations (see file :ref:`DISPLACEMENTS<displacements>` for details)
to calculate amplitude changes.

Delta Calculation in ViPErLEED
------------------------------

When a delta calculation is requested, ViPErLEED will load the :ref:`tensor files<tensorszip>` and read the :ref:`DISPLACEMENTS file<displacements>`.
Both are checked for consistency before proceeding.
The :ref:`PHASESHIFTS file<phaseshifts>` (which needs to include phaseshifts for chemical substitutions) will be loaded as well or generated if required.

ViPErLEED will then determine for which atoms delta-amplitudes need to be
generated. If it is found that a previously generated delta file (stored
in the ``Deltas`` directory) already contains some of the requested amplitudes, their
computation will be skipped and a message will be put in the log accordingly.

As for the :ref:`refercence Calculation<ref-calc>`, ViPErLEED will then move the needed TensErLEED source files to a temporary directory and compile them **at run-time**.
The input files (based on the :ref:`DISPLACEMENTS fil<displacements>` and the system symmetry) will be written in the format expected by TensErLEED and calculations for each atom will be executed independently.
Again, the parameter :ref:`N_CORES<ncores>` determines how many processes are run concurrently.

Once finished, ViPErLEED will pack the resulting :ref:`delta files<deltaszip>` into a ``.zip`` archive in the ``Deltas`` directory.