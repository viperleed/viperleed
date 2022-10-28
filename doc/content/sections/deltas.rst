.. _sec_deltas:

=================
Delta Calculation
=================

The delta amplitude calculation (called with :ref:`RUN = 2<run>`) is the 
second part of a :ref:`Tensor LEED<tensor_leed>` calculation as implemented 
in ViPErLEED.
It requires a :ref:`refercence Calculation<ref-calc>` to be completed 
beforehand and :ref:`tensor files<tensorszip>` to be stored.

In the delta calculation, these :ref:`tensor files<tensorszip>`, i.e. the
tensor quantities :math:`T^{ref}_{i;l,m;l',m'}`, that only depend on the 
refercence structure are combined with the requested pertubations 
(see file :ref:`DISPLACEMENTS<displacements>` for details)
to calculate amplitude changes.

The individual pertubations to the reference structure may be (arbitrary combinations of) 
geometrical  displacements, changes in the vibrational amplitudes or 
chemical substitutions. 
All of these pertubations are considered on an atom-by-atom basis. For each
atom :math:`i` and for each requested pertubation :math:`n` to that atom,
first the perturbed :math:`t`-matrix :math:`\tilde{t_i} = t_i + \delta \tilde{t_i}` and then the 
expression

.. math:: 

    \delta \tilde{A}_{i,n}^{per} = \sum_{l,m;l',m'} T^{ref}_{i;l,m;l',m'} \braket{\vec{r_i},l,m| \delta t_{i,n} |\vec{r_i},l',m'}

are evaluated to calculate linearly independent amplitude changes 
:math:`\delta \tilde{A}_{i,n}^{per}`.

The resulting delta-amplitudes are stored in the :ref:`delta files<deltaszip>`
and will be used in the :ref:`structure search<sec_search>` to calculate
the intensities and subsequently the :ref:`R-factor<r-factor_calculation>` 
for each structure candidate. :cite:p:`blumFastLEEDIntensity2001a`

.. note:: 
    Depending on the size of the unit cell and the requested pertubations,
    the parameter space (and the :ref:`delta files<deltaszip>`) may become
    very big.

Delta Calculation in ViPErLEED
------------------------------

When a delta calculation is requested, ViPErLEED will load the 
:ref:`tensor files<tensorszip>` and read the :ref:`DISPLACEMENTS file<displacements>`.
Both are checked for consistency before proceeding.
The :ref:`PHASESHIFTS file<phaseshifts>` (which needs to include phaseshifts for 
possible chemical substitutions) will be loaded as well or generated if 
required.

ViPErLEED will then determine for which atoms delta-amplitudes need to be
generated. If it is found, that a previously generated delta file (stored
in the ``Deltas`` directory) already contains some of the requested amplitudes, their
computation will be skipped and a message will be put in the log accordingly.

As for the :ref:`refercence Calculation<ref-calc>`, ViPErLEED will then 
move the needed TensErLEED source files to a temporary directory and 
compile them **at run-time**.
The input files will be written in the format expected by TensErLEED and 
calculations for each atom will be executed independently.
Again, the parameter :ref:`N_CORES<ncores>` determines how many
processes are run concurrently.

Once finished, ViPErLEED will pack the resulting :ref:`delta files<deltaszip>`
into a ``.zip`` archive in the ``Deltas`` directory.