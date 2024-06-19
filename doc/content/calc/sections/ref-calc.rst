.. include:: /substitutions.rst

.. _ref-calc:

=====================
Reference calculation
=====================

ViPErLEED, which is based on the :term:`TensErLEED` package, follows the same
scheme as TensErLEED for the calculation of diffraction intensities. Taking the
reference structure (:ref:`POSCAR` file) and precomputed scattering
phase shifts (\ :ref:`PHASESHIFTS` file) as input, a so-called
reference calculation is performed. In the code, the reference calculation
is often shortened to ``refcalc``.

Following the reference calculation, a local structure optimization can
be performed using the :ref:`tensor-LEED approach<tensor_leed>`. This is
accomplished by running a :ref:`delta-amplitude calculation<sec_deltas>`
and a :ref:`structure search<sec_search>`.

Reference calculation in ViPErLEED
----------------------------------

The reference calculation in ViPErLEED is implemented as a smart wrapper
for the TensErLEED reference calculation and is called by setting ``RUN = 1``
(see also the :ref:`RUN` parameter). Calculations for each energy step
(defined by :ref:`THEO_ENERGIES`) are (essentially) independent from one
another. Thus, the computation time scales roughly linearly with the total
number of energy steps.\ [1]_ Additionally, :ref:`ncores` energy steps for the
reference calculation are executed simultaneously.

.. note::
    The reference calculation can be memory intensive for large unit cells.
    If using a large :ref:`ncores` value, make sure not to run into memory
    limitations.

For each energy step, ViPErLEED determines the angular-momentum quantum number
|lmax| at which the spherical-harmonics expansion can be truncated. See
parameters :ref:`LMAX` and :ref:`PHASESHIFTMIN` for details.
ViPErLEED will then compile the required TensErLEED source files for all
required values of |lmax| **at run time** in temporary directories named
:file:`refcalc-compile_LMAX<lmax>`. A :file:`refcalc_compile_LMAX<lmax>.log`
:ref:`log file<log_files>` is stored in the :file:`SUPP/compile_logs`
folder for each |lmax| value.

ViPErLEED will then go through all required energies from highest
to lowest, performing calculations in temporary directories named
:file:`refcalc-part_<energy>eV`. The output amplitudes, intensities
and tensors are stored in the work directory before further processing.
A :ref:`log file<log_files>` containing messages emitted by the
:program:`ref-calc.f` TensErLEED Fortran program can be found under
:file:`SUPP/refcalc-<timestamp>.log`.

Once the TensErLEED reference calculation has concluded, ViPErLEED collects
all files, removes temporary directories, and combines the results into a
:ref:`THEOBEAMS.csv<THEOBEAMS>` file. The calculated beams are also plotted
for inspection in :ref:`THEOBEAMS.pdf<theobeams>`.

Finally, unless the tensor output was disabled via :ref:`TENSOR_OUTPUT`,
ViPErLEED collects the :ref:`tensor files<tensorszip>` created by TensErLEED.
They are stored, in compressed form, in the :file:`Tensors` directory.


.. todo::

    Check if note [1] is indeed correct by experimenting. I think it is as the
    number of beams scales with the cross section of the Ewald sphere, whose
    diameter scales with sqrt(E). The only real question is how does the number
    of non-evanescent beams scales with Emax.


.. [1] Notice, however, that larger energies involve a larger number
       of scattered beams. They are therefore slower: the size of
       :math:`d\times{}d` matrices containing scattered beams scales as
       :math:`d \sim{}E_\mathrm{max}` at fixed |lmax|. The computation
       time thus scales as :math:`d^{n+1} = E_\mathrm{max}^{n+1}`,
       with :math:`n=2.38`–\ :math:`3` due to matrix inversion. This means
       that the computation time scales roughly with the fourth power of
       the maximum value of the energy in the calculation.
