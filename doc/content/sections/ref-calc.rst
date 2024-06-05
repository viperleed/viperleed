.. _ref-calc:

=====================
Reference calculation
=====================

ViPErLEED, which is based on the TensErLEED package, follows
the same scheme as TensErLEED for the calculation of diffraction intensities.
Taking the reference structure (:ref:`POSCAR file<POSCAR>`) and pre-computed
scattering phase shifts (\ :ref:`PHASESHIFTS file<phaseshifts>`) as input,
a so-called reference calculation is performed.
In the code, the reference calculation is often shortened to ``ref-calc``.

Following the reference calculation, a local structure optimization can be
performed using the :ref:`tensor-LEED approach<tensor_leed>`.
This is accomplished by running a :ref:`delta amplitude calculation<sec_deltas>`
and a :ref:`structure search<sec_search>`.

Reference calculation in ViPErLEED
----------------------------------

The reference calculation in ViPErLEED is implemented as a smart wrapper
for the TensErLEED reference calculation and called by setting the :ref:`RUN`
parameter = 1. Calculations for each energy step (defined by
:ref:`THEO_ENERGIES<theo_energies>`) are  performed **independently** and
thus the computation time scales roughly linearly with the number of total
energy steps. Additionally, :ref:`N_CORES<ncores>` reference calculations
are executed simultaneously.

.. note::
    The reference calculation can be memory intensive for large unit cells.
    If using a large :ref:`N_CORES<ncores>`, make sure you are not running
    into memory limitations.

For each energy step, ViPErLEED will determine at which angular momentum
quantum number :math:`l_{max}` the sums can be truncated. See parameters
:ref:`LMAX` and :ref:`PHASESHIFT_EPS<PHASESHIFTMIN>` for details.
ViPErLEED will then, in a first step, compile the required TensErLEED source
files for all required values of :math:`l_{max}` **at run-time** in temporary
directories called ``refcalc-compile_LMAX=n``.

ViPErLEED will then go through all required energy steps from highest
energy to lowest, performing  calculations in temporary directories called
``refcalc-part_xxeV``. The output amplitudes, intensities and tensors are
stored in the work directory before further processing. A log will be written
to ``refcalc-$timestamp.log``.

Once the TensErLEED reference calculation has concluded, ViPErLEED will
collect all files, remove temporary directories and combine the results
into a :ref:`THEOBEAMS.csv<THEOBEAMS>` file. By default, the theoretical
beams will also be plotted for inspection in :ref:`THEOBEAMS.pdf<theobeams>`.

Finally, unless the tensor output was disabled with the
:ref:`TENSOR_OUTPUT<toutput>` parameter, ViPErLEED will
collect the created :ref:`tensor files<tensorszip>` in
compressed form in the ``Tensors`` directory.
