.. include:: /substitutions.rst

.. _tensor_leed:

========================
The tensor-LEED approach
========================

Here, we provide a very rudimentary introduction to the tensor LEED
approach employed by TensErLEED and consequently also ViPErLEED. Please
note that this neither is, nor aims to be a comprehensive or rigorous
introduction to the topic. The descriptions below are only intended to
provide a quick overview of the method and serve as explanation and
motivation for the various sections of a |LEED-IV| calculation in ViPErLEED.

For an in-depth description of all parts of the tensor LEED approach,
we refer to the papers by Rous and Pendry (see
Refs. :cite:p:`rousTensorLEEDTechnique1986,rousTheoryTensorLEED1989`)
and the explanation of how the approach is implemented in TensErLEED
in the original publication of the TensErLEED software
:cite:p:`blumFastLEEDIntensity2001a`.

Reference calculation
=====================

In the :ref:`reference calculation<ref-calc>`, the full dynamic scattering
(i.e., including all multiple-scattering contributions) of the incident
electron wave (with complex amplitude :math:`A_{in}`) at the reference
structure is calculated. In principle, this calculation only yields the
scattering amplitudes :math:`A_{\mathrm{out}}^{\mathrm{ref}}` (and intensities
via :math:`I = |A|^2`) of all diffracted beams of interest for the requested
reference structure.
However, as shown by Rous and Pendry :cite:p:`rousTensorLEEDTechnique1986`,
using a first-order perturbation theory approach, it is possible to obtain
rather accurate diffraction amplitudes for small deviations from this reference
structure. These deviations may be geometrical (altered atom positions),
changes to the vibrational amplitude or chemical substitutions
:cite:p:`kottckeNewApproachAutomated1997`.

Each atom :math:`i` is assigned an atomic :math:`t`-matrix, :math:`t_i`, based
on phase shifts and positions within the unit cell. The perturbed structure is
consequently characterized by altered atomic
:math:`t`-matrices :math:`\tilde{t_i} = t_i + \delta \tilde{t_i}`.

In this case, the diffraction amplitudes for beam :math:`n` a perturbed
structure can be written  as the reference amplitudes plus a sum of delta
amplitudes for the altered atoms (index :math:`i`):

.. math::

    \tilde{A}^{\mathrm{per}}_{n} = A^{\mathrm{ref}}_{n} + \sum_{i} \delta \tilde{A}_{i,n}^{\mathrm{per}}

These delta-amplitudes can be expressed as

.. math::

    \delta \tilde{A}_{i,n}^{\mathrm{per}} = \sum_{l,m;l',m'} T^{\mathrm{ref}}_{i,n;l,m;l',m'} \braket{\vec{r_i},l,m| \delta t_{i,n} |\vec{r_i},l',m'}

using the perturbed atomic :math:`t`-matrices :math:`\delta t_i` and the
tensor quantities :math:`T^{\mathrm{ref}}_{i,n;l,m;l',m'}`. The sum runs over
angular momentum and magnetic quantum numbers :math:`l` and :math:`m`.
For a more rigorous derivation, refer to the original work by Rous and Pendry
:cite:p:`rousTensorLEEDTechnique1986` and the TensErLEED paper by Blum and
Heinz :cite:p:`blumFastLEEDIntensity2001a`.

The quantities :math:`T^{\mathrm{ref}}_{i;l,m;l',m'}` only depend on the
reference structure and are commonly just referred to as tensors.
Importantly, the tensors can be calculated in the reference calculation,
and saved in the :ref:`tensor files<tensorszip>`. They are the starting
point for the subsequent :ref:`delta amplitude calculation<sec_deltas>`
and :ref:`structure search<sec_search>`.


Delta-amplitude calculation
===========================

The individual perturbations to the reference structure may be (arbitrary
combinations of) geometrical  displacements, changes in the vibrational
amplitudes or chemical substitutions. As tensor LEED is based on first-order
perturbation theory, we can treat these perturbations and the resulting
amplitude changes on an atom-by-atom basis.

For each
atom :math:`i` and for each requested perturbation :math:`p` to that atom,
first the perturbed :math:`t`-matrix
:math:`\tilde{t}_{i,p} = t_i + \delta \tilde{t}_{i,p}`
and then the expression

.. math::

    \delta \tilde{A}_{i,n,p}^{\mathrm{per}} = \sum_{l,m;l',m'} T^{\mathrm{ref}}_{i,n;l,m;l',m'} \braket{\vec{r_i},l,m| \delta \tilde{t}_{i,n,p} |\vec{r_i},l',m'}

are evaluated to calculate amplitude changes
:math:`\delta \tilde{A}_{i,n,p}^{\mathrm{per}}`.

The resulting delta-amplitudes are stored in the :ref:`delta files<deltaszip>`
and will be used in the :ref:`structure search<sec_search>` to calculate the
intensities and subsequently the :math:`R` :ref:`factor<r-factor_calculation>`
for each structure candidate. :cite:p:`blumFastLEEDIntensity2001a`

.. note::
    Depending on the size of the unit cell and the requested perturbations,
    the parameter space (and the :ref:`delta files<deltaszip>`) may become
    very big.

.. _tensor_leed_search:

Structure search
================

Once the amplitude changes for all required perturbations have been
obtained, the final diffraction amplitudes can be calculated using a
simple superposition. Essentially, for any perturbed structure, we
compute the amplitudes by simply summing up amplitude changes (deltas)
for all affected atoms.

Consequently, using these resulting amplitudes and intensities,
an |R factor| vs. the experimental intensities can now be obtained
for any structure in the configuration space. Then, the best-fit
structure must be found by an optimization (minimization of the
R factor) in the configuration space.


While conceptually simple, this optimization can be practically and
computationally very challenging, and usually constitutes the
computationally most expensive part of a |LEED-IV| calculation.
Still, using the TensErLEED approach, the problem generally remains
tractable, even for relatively large unit cells. Running a full-dynamic
calculation for every configuration is usually orders of magnitude more
expensive. :cite:p:`rousTensorLEEDTechnique1986`

That being said, there remain some fundamental caveats to the structure
optimization in the tensor LEED approximation and also |LEED IV| in general:

-   Since the tensor LEED method is a perturbative approach, it only works
    reliably for *small* perturbations. What constitutes a *small* perturbation
    is naturally system-dependent, but generally, the limit lies in the range
    of 0.2 Å to 0.3 Å at best :cite:`rousTensorLEEDTechnique1986`. For larger
    displacements the search might still give the right trends but equally
    might be misleading.

    To extend the range of the structural search, it is possible to run a new
    reference calculation and delta-amplitudes calculation when the structure
    optimization trajectory approaches this limit. You can use the the
    :ref:`RUN parameter<run>` to chain multiple reference calculations,
    delta-amplitude calculations and structure searches together.

-   The parameter space grows quickly for larger unit cells.
    Luckily, many symmetries inherent to the surface structure can be exploited
    to eliminate redundant parameters. For example, geometric displacements of
    symmetry-linked atoms must always happen in a concerted fashion. If that
    were not the case, the symmetry would be broken and a different LEED
    pattern would result.

    To make use for these symmetries and the resulting reduced parameter space,
    it is necessary to know and enforce the surface slab symmetry. While
    manually finding out the surface slab symmetry is generally an easy task,
    enforcement is not. This would require manually going over every
    symmetry-linked atom and defining matching displacement vectors.

    *Fortunately for the user*, automatic symmetry-detection and enforcement
    is one of the **main features** of ViPErLEED.


-   When using Pendry's |R factor|, the |R factor| hyper-surfaces tend to be
    inherently non-smooth :cite:p:`rousTensorLEEDApproximation1992`. Users
    should be aware that local minima are possible and that the optimization
    algorithm might get stuck in these minima if the parameter space is not
    opened up sufficiently.

-   As described above, the tensor LEED implementation in TensErLEED separates
    the calculation of delta-amplitudes and the structure optimization into two
    mostly independent stages. As a direct consequence, the optimization can
    **only** be performed on a pre-defined grid of perturbation vectors (as
    given by the :ref:`DISPLACEMENTS file<displacements>`). Further, to achieve
    the best possible fit, the grid based nature makes it necessary to run
    multiple sets of delta-amplitude calculations and structure optimizations
    with increasingly finer grids.

    .. note::
        Starting with a fine grid over a large variation range is not
        recommended since too many grid points per parameter will significantly
        slow down the fit.

-   The structure search implemented in TensErLEED has the additional
    limitation that geometrical displacements are limited to one dimension
    per atom. Per search run, atoms can only be displaced along a pre-defined
    parametrized curve, rather than freely in 3D space. To optimize the
    position of atoms in 3 dimensions, multiple sequential search runs are
    required. See the entry on the :ref:`DISPLACEMENTS file<displacements>`
    for details and workarounds (such as looping searches).

Optimization algorithm
======================

.. _optimization_algorithm:

The rough |R-factor| surface, together with its grid-based nature greatly
limits the pool of applicable optimization algorithms. TensErLEED employs
a modified random sampling strategy with a down-step criterion as described
by Kottcke and Heinz :cite:p:`kottckeNewApproachAutomated1997`.
The optimization is performed in parallel for a set of individuals
(= independent parameter combinations), as defined by the parameter
:ref:`SEARCH_POPULATION<searchpop>`. The starting points for the
optimization is defined by :ref:`SEARCH_START<searchstart>`.

For each search step (called "generation" based on the terminology of genetic
algorithms), a new grid point in the parameter space is selected *randomly*,
but based on a probability distribution centered on the current position.
The |R factor| is calculated for the selected parameter combination and the new
parameter set is accepted **only if** the |R factor| for the new configuration
is lower then for the previous configuration. The width of the probability
distribution is determined by the current |R factor| and the parameter
:ref:`SEARCH_CONVERGENCE<search_convergence>` (in particular, the ``gaussian``
flag).

ViPErLEED enables more sophisticated control over the search process than
is possible with TensErLEED alone. Different types of convergence criteria
and an automatic scaling of the probability distribution can be set using
:ref:`SEARCH_CONVERGENCE<search_convergence>`. Furthermore, as defined by
the parameter :ref:`SEARCH_CULL<search_cull>`, whenever
:ref:`partial convergence<search_convergence>` is reached, a portion of the
search population can be dropped and re-initialized to get out of local minima.
By default, the search population is partially re-initialized using a custom
genetic algorithm (see :ref:`SEARCH_CULL<search_cull>` for details).
