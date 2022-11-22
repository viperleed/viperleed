.. _tensor_leed:

========================
The Tensor LEED approach
========================

Here, we provide a very rudimentary introduction to the tensor LEED
approach employed by TensErLEED and consequently also ViPErLEED.
Please note, that this is not, and does not aim to be a comprehensive or rigourous 
introduction to the topic.
The descriptions below are only intended to provide a quick overview of 
the method and serve as explanation and motivation for the various sections of 
a LEED :math:`I(V)` calculation in ViPErLEED.

For an in-depth description of all parts of the tensor LEED approach, we refer to the 
papers by Rous and Pendry (see Refs. :cite:p:`rousTensorLEEDTechnique1986,rousTheoryTensorLEED1989`)
and the explaination of how to approach is implemented in TensErLEED in 
the software's original publication (see Ref. :cite:p:`blumFastLEEDIntensity2001a`).

Reference Calculation
=====================

In the :ref:`reference calculation<ref-calc>`,
the full dynamic, multi-scattering of the 
incident electron wave (with complex amplitude :math:`A_{in}`) at the reference 
structure is calculated. 
In principle, this calculation only yields the scattering amplitudes
:math:`A_{out}^{ref}` (and intensities) for the requested reference structure.
However, as shown by Rous and Pendry :cite:p:`rousTensorLEEDTechnique1986`, 
using a first-order pertubation theory approach, it is possible to obtain
accurate diffraction amplitudes for small deviations from this refercence structure.
These deviations may be geometrical (altered atom positions), changes to 
the vibrational amplitude or chemical substitutions.

Each atom :math:`i` is assigned
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
:cite:p:`rousTensorLEEDTechnique1986` and the TensErLEED paper by Blum and 
Heinz :cite:p:`blumFastLEEDIntensity2001a`.

The quantities :math:`T^{ref}_{i;l,m;l',m'}` only depend on the reference structure
and are commonly just referred to as tensors.
Importantly, the tensors can be calculated in the refercence calculation, 
and saved in the :ref:`tensor files<tensorszip>`. 
They are the starting point for the subsequent :ref:`delta amplitude calculation<sec_deltas>`
and :ref:`structure search<sec_search>`.


Delta-Amplitude Calculation
===========================

The individual pertubations to the reference structure may be (arbitrary combinations of) 
geometrical  displacements, changes in the vibrational amplitudes or 
chemical substitutions.
As tensor LEED is based on first-order pertubation theory approach,
we can treat these pertubations on an atom-by-atom basis with the resulting 
amplitude changes being considered linearly independent.

For each
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

.. _tensor_leed_search:

Structure Search
================

Once the amplitude changes for all required pertubations have been obtained,
the final diffraction amplitudes can be calculated using a simple combination.
Essentially, for any perturbed structure, we compute the amplitudes by 
simply summing up amplitude changes (deltas) for all affected atoms.

Consequently, using these resulting amplitudes (and intensities via :math:`I = |A|^2`), 
an R-factor vs. the experimental intensities can now be obtained for any
structure in the configuration-space. 
All that is left then, to finding the best-fit structure is an optimization
in the configuration space over the R-factor.

While conceptually simple, this optimization can be practically and computationally 
very challenging, and generally constitutes the computationally most expensive
part of a LEED :math:`I(V)` calculation. Still, using the TensErLEED approach,
the problem generally remains tractable, even for relatively large unit cells.
Running a full-dynamic calculation for every configuration is usually orders
of magnitude more expensive. :cite:p:`rousTensorLEEDTechnique1986`

That being said, there remain some fundamental caveats to the structure optimization 
in the tensor LEED approximation and also LEED :math:`I(V)` in general:

-   Since the tensor LEED method is a perturbative approach, it only works reliably for
    *small* pertubations.
    What constitutes a *small* pertubation is naturally system-dependent, but 
    generally, the limit lies in the range of 200 to 300 pm at best :cite:`rousTensorLEEDTechnique1986`.

    To extend the range of the structural search, as work-around, it is possible to run a new refercence calculation and delta-amplitudes calculation when the structure optimization trajectory approaches this limit.
    You can use the the :ref:`RUN parameter<run>`
    to execute multiple reference calculations, delta-amplitude calculations,
    and structure searches in series.

-   The parameter space grows quickly for larger unit cells.
    Luckily, many symmetries inherent to the surface structure can be exploited to eliminate redundent parameters.
    For example, geometric displacements of symmetry-linked atoms must always happen in a concerted fashion.
    If that were not the case, the symmetry would be broken and usually\ [.#] a different LEED pattern would result.

    To make use for these symmetries and the resulting reduced parameter space, it is necessary to know and enforce the surface slab symmetry.
    While manually finding out the surface slab symmetry is generally an easy task, enforcement is not.
    This would require manually going over every symmetry-linked atom and defining matching displacement vectors.

    *Fortunately for the user*, automatic symmetry-detection and enforcement is one of the **main features** of ViPErLEED.
    See the ViPErLEED paper for details (**TODO**).

-   The R-factor hyper-surfaces tend to be inherently non-smooth. This is 
    a consequence of how the various R-factors are designed.

-   As described above, the tensor LEED implementation in TensErLEED separates the calculation of delta-amplitudes and the structure optimization into two mostly-indepenent stages.
    As a direct consequence, the optimization can **only** be performed on a pre-defined grid of pertubation vectors (as given by the :ref:`DISPLACEMENTS file<displacements>`).
    Further, to achieve the best possible fit, this grid
    generally makes it necessary to run multiple sets of delta-amplitude
    calculations and structure optimizations with varying grid-densities.

-   The structure search implemented in TensErLEED has the additional limitation that geometrical displacements are limited to one dimension per atom.
    Per search run, atoms can only be displaced along a pre-defined parametrized curve, rather than freely in 3D space.
    To optimize the position of atoms in 3 dimensions, multiple sequential search runs are required.
    See the entry on the :ref:`DISPLACEMENTS file<displacements>` for details and work-arounds (such as looping searches).
    
-    The rough R-factor surface, together with its grid-based nature greatly limits the pool of applicable optimzation algorithms.
    TensErLEED employes the strategy as descirbed by Kottcke and Heinz :cite:p:`kottckeNewApproachAutomated1997`.
    In each search step (i.e. "generation")...



.. [1] There are exceptions, in which the same LEED pattern can result. For example, on an fcc(111) surface, a (2:math:`\times`2) reconstruction and a (1:math:`\times`2) with domains would give the same qualitative pattern.
.. [2] For details on the algorithm used in TensErLEED, see ref. :cite:p:`kottckeNewApproachAutomated1997`.
