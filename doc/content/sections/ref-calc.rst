.. _ref-calc:

=====================
Reference calculation
=====================

The goal of any :term:`LEED-I(V)` calculation is the calculation of 
energy dependent electron scattering amplitudes and intensities. These 
intensity curves (often referred to as :math:`I(V)` curves or spectra) 
are very sensitive to the precise position and vibrational amplitudes in
the surface unit cell.
For more details consult works that cover the basics of :term:`LEED` and 
:term:`LEED-I(V)`, e.g.
Chapter 4 in :cite:t:`fausterSurfacePhysicsFundamentals2020` or 
the overview by:cite:t:`heinzElectronBasedMethods2013`.
For computational details please have a look at the ViPErLEED paper 
(**TODO**) and the original work describing TensErLEED by Blum and Heinz 
:cite:p:`blumFastLEEDIntensity2001a`.
In the code, the refercence calculation is often shorted to ``ref-calc``.


ViPErLEED, which is built upon the TensErLEED package as a backend, follows
the same scheme as TensErLEED for the calculation of diffraction intensities.
Taking the reference structure (:ref:`POSCAR file<POSCAR>`) and pre-computed 
scattering phaseshifts (:ref:`PHASESHIFTS file<phaseshifts>`)as input, 
a so called refercence calculation is performed. 

Following the refercence calculation, a local strucutre optimization can be 
performed using the Tensor LEED approach summarized :ref:`below<tensor_leed>`.
This is accomplished by running a :ref:`delta amplitude calculation<sec_deltas>`
and a :ref:`structure search<sec_search>`.

.. _tensor_leed:

The Tensor LEED approach
------------------------

In the reference calculation, the full dynamic, multi-scattering of the 
incident electron wave (with complex amplitude :math:`A_{in}`) at the reference 
structure is calculated. 
In principle, this calculation only yields the scattering amplitudes
:math:`A_{out}^{ref}` (and intensities) for the requested reference structure.
However, as shown by Rous and Pendry :cite:p:`rousTensorLEEDTechnique1986`, 
using a first-order pertubation theory approach, it is possible to obtain
accurate diffraction amplitudes for small deviations from this refercence structure.

These deviations may be geometrical (altered atom positions), changes to 
the vibrational amplitude or chemical substitutions.
In the reference calculation, each atom :math:`i` is assigned
an atomic :math:`t`-matrix, :math:`t_i` based on phaseshifts and positions within the unit cell.
The perturbed structure is consequently characterized by altered atomic 
:math:`t`-matricies :math:`\tilde{t_i} = t_i + \delta \tilde{t_i}`.

In this case, the diffraction amplitudes for a perturbed structure can be written 
as the refercence amplitudes plus a sum of delta-amplitudes for the 
altered atoms (index :math:`i`):

.. math:: 

    \tilde{A}^{per} = A^{ref} + \sum_{i} \delta \tilde{A}_{i}^{per}

These delta-amplitudes can be expressed as 

.. math:: 

    \delta \tilde{A}_{i}^{per} = \sum_{l,m;l',m'} T^{ref}_{i;l,m;l',m'} \braket{\vec{r_i},l,m| \delta t_i |\vec{r_i},l',m'}

using the perturbed atomic :math:`t`-matricies :math:`\delta t_i` and the
tensor quantities :math:`T^{ref}_{i;l,m;l',m'}`. The sum runs over angular 
momentum and magnetic quantum numbers :math:`l` and :math:`m`.
For a more rigourous derivation, refer to the original work by Rous and Pendry 
:cite:p:`rousTensorLEEDTechnique1986` and the TensErLEED paper by Bulm and 
Heinz :cite:p:`blumFastLEEDIntensity2001a`.

The quantities :math:`T^{ref}_{i;l,m;l',m'}` only depend on the reference structure
and are commonly just referred to as tensors.
Importantly, the tensors can be calculated in the refercence calculation, 
and saved in the :ref:`tensor files<tensorszip>`. 
They are the starting point for the subsequent :ref:`delta amplitude calculation<sec_deltas>`
and :ref:`structure search<sec_search>`.


Reference Calculation in ViPErLEED
----------------------------------

The refercence calculation in ViPErLEED is implemented as a smart wrapper 
for the TensErLEED refercence calculation and called by setting the :ref:`RUN` parameter = 1.
Calculations for each energy step (defined by :ref:`THEO_ENERGIES<REFENERGIES>`) are 
performed **independently** and thus the computation time scales roughly 
linearly with the number of total energy steps.
Additionally, :ref:`N_CORES<ncores>` refercence calculations are executed 
simultaneously.

.. note:: 
    The refercence calculation can be memory intensitve for large unit cells.
    If using a large :ref:`N_CORES<ncores>` make sure you are not running into 
    memory limitations.

For each energy step, ViPErLEED will determine at which angular momentum 
quantum number :math:`l_{max}` the sums can be truncated. See parameters
:ref:`LMAX` and :ref:`PHASESHIFT_EPS<PHASESHIFTMIN>` for details.
ViPErLEED will then, in a first step, compile the required TensErLEED source 
files for all required values of :math:`l_{max}` **at run-time** in temporary
directories called ``refcalc-compile_LMAX=n``.

ViPErLEED will then go through all required energy steps from highest 
energy to lowest, performing  calculations in temporary directories called 
``refalc-part_xxeV``. The raw TensErLEED input and output files 
:ref:`refcalc-fin` and :ref:`fd.out<fd_out>` are found in these directories
together with the executables.
A log will be written to ``refcalc-$timestamp$.log``.

Once the TensErLEED refercence calculation has concluded, ViPErLEED will
collect all files, remove temporary directories and combine the results 
into a :ref:`THEOBEAMS.csv<THEOBEAMS>` file. By default, the theoretical 
beams will also be plotted for inspection in :ref:`THEOBEAMS.pdf<theobeams>`.

Finally, unless the tensor output was disabled with the :ref:`TENSOR_OUTPUT<toutput>`
parameter, ViPErLEED will collect the created :ref:`tensor files<tensorszip>`
in the ``Tensors`` directory.

