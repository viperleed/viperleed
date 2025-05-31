.. include:: /substitutions.rst

.. _tensor_leed:

========================
The tensor-LEED approach
========================

This section provides a rudimentary introduction to the tensor-LEED approach
employed by TensErLEED and, consequently, also ViPErLEED. Note that this
is neither, nor aims to be a comprehensive or rigorous introduction to the
topic. The descriptions below are only intended to provide a quick overview
of the method, and serve as explanation and motivation for the various sections
of a |LEED-IV| calculation in ViPErLEED.

An in-depth description of all parts of the tensor-LEED approach is presented
in Refs. :cite:alp:`rousTensorLEEDTechnique1986,rousTheoryTensorLEED1989`.
The original publication of the TensErLEED software explains 
how the tensor-LEED approach is implemented in TensErLEED
:cite:p:`blumFastLEEDIntensity2001a`.

Reference calculation
=====================

The :ref:`reference calculation<ref-calc>` determines the full-dynamic
scattering (i.e., including all multiple-scattering contributions) of an
electron wave incident at a "reference" structure. This calculation yields
the complex scattering amplitudes :math:`A^{\mathrm{ref}}_{\mathbf{g}}` (and
intensities, via :math:`I = |A|^2`) of all diffracted beams :math:`\mathbf{g}`
that are of interest for the requested structure.

As discussed in more detailed elsewhere
:cite:p:`heinzElectronBasedMethods2013,kraushoferViPErLEEDPackageCalculation2025`,
the full-dynamic calculation is computationally demanding. However, it
is possible to obtain accurate diffraction amplitudes for small deviations
from the reference structure by using a first-order-perturbation approach
:cite:p:`rousTensorLEEDTechnique1986`. The deviations from the
reference structure may be geometric (i.e., altered atom positions),
changes to atomic vibration amplitudes, or chemical substitutions
:cite:p:`kottckeNewApproachAutomated1997`.

Each atom :math:`i` is assigned a "\ :math:`t` matrix", :math:`t_i`, based on
electron-scattering phase shifts and positions within the unit cell.
The perturbed structure is consequently characterized by altered atomic
:math:`t` matrices, :math:`\tilde{t_i} = t_i + \delta t_i`.

The amplitude for a beam :math:`\mathbf{g}`  diffracted at the perturbed
structure can be written as

.. math::

    A^{\mathrm{per}}_{\mathbf{g}} = A^{\mathrm{ref}}_{\mathbf{g}} + \sum_{i} \delta{}A_{i,\mathbf{g}} ,

that is, as the reference amplitudes plus a sum of delta amplitudes for
the altered atoms. These delta amplitudes can be expressed as

.. math::

    \delta{}A_{i,\mathbf{g}} = \sum_{l,m;l',m'} T^{\mathrm{ref}}_{i,\mathbf{g};l,m;l',m'} \braket{\mathbf{r}_i,l,m| \delta{}t_{i} |\mathbf{r}_i,l',m'}

using the perturbations of the atomic :math:`t` matrices, :math:`\delta t_i`,
the tensorial quantities :math:`T^{\mathrm{ref}}_{i,\mathbf{g};l,m;l',m'}`, and
the unperturbed positions of the atoms :math:`\mathbf{r}_i`. The sum runs over
two sets of angular-momentum (\ :math:`l`, :math:`l'`) and magnetic
(\ :math:`m`, :math:`m'`) quantum numbers. For a rigorous derivation, refer to
the original work by :cite:t:`rousTensorLEEDTechnique1986` and to the
TensErLEED paper by :cite:t:`blumFastLEEDIntensity2001a`.

The quantities :math:`T^{\mathrm{ref}}_{i,\mathbf{g};l,m;l',m'}` only depend
on the reference structure and are commonly just referred to as "tensors".
Importantly, the tensors can be calculated during the reference calculation.
They are the starting point for the subsequent
:ref:`delta amplitude calculation<sec_deltas>`
and :ref:`structure search<sec_search>`.


Delta-amplitude calculation
===========================

The individual perturbations to the reference structure may be combinations of
geometric displacements, changes in the vibration amplitudes, or chemical
substitutions. As tensor LEED is based on first-order perturbation theory,
these perturbations — and the resulting amplitude changes — can be treated
on an atom-by-atom basis.

For each atom :math:`i` and for each requested perturbation :math:`p` to that
atom, the delta-amplitude calculation evaluates the perturbed :math:`t` matrix
:math:`\tilde{t}_{i,p} = t_i + \delta t_{i,p}`
and the corresponding amplitude changes

.. math::

    \delta{}A_{i,\mathbf{g},p} = \sum_{l,m;l',m'} T^{\mathrm{ref}}_{i,\mathbf{g};l,m;l',m'} \braket{\mathbf{r}_i,l,m| \delta t_{i,p} |\mathbf{r}_i,l',m'} .

The resulting delta-amplitudes are used in the
:ref:`structure search<sec_search>` to calculate the
perturbed intensities for each structure candidate
:cite:p:`blumFastLEEDIntensity2001a`.

.. note::
    Depending on the size of the unit cell and the requested
    perturbations, the parameter space may become very big.

.. _tensor_leed_search:

Structure search
================

.. todo:
    Add info on parameter redundancy and how many parameters should/can be
    optimized and fitted.
    Mention rule of thumb of ~100eV energy range per parameter (because of
    parameter correlation).
    Mention that everything optimized counts as a parameter, i.e. geometry,
    vibrations and chemical perturbations but also V0r shift, incidence angles,
    and all FD parameters.

Once the amplitude changes for all required perturbations have been
obtained, the final diffraction amplitudes can be calculated using
a simple superposition: the overall amplitude change is the sum

.. math::

    \delta{}A^\mathrm{tot}_\mathbf{g} = \sum_{i,p} \delta{}A_{i,\mathbf{g},p}

of the amplitude changes for all the affected atoms. The total diffracted
amplitudes then result by adding the amplitude changes to the reference
amplitudes.

This way, scattered intensities can be obtained for any structure in the
configuration space. Different candidate structures in this configuration
space are compared to experimental data via the
:math:`R` :ref:`factor<r-factor_calculation>`. The best-fit structure
results from a minimization of the |R factor| in the configuration space.

While conceptually simple, this optimization is practically very challenging,
and usually constitutes the computationally most expensive part of a |LEED-IV|
calculation. By using the tensor-LEED approach, the problem is tractable, even
for systems with relatively large unit cells. Running a full-dynamic
calculation for every configuration is orders of magnitude more
computationally expensive :cite:p:`rousTensorLEEDTechnique1986`.

There are some fundamental caveats concerning the structure optimization in
the tensor-LEED approximation:

-   Since the tensor-LEED method is perturbative, it only works reliably for
    *small* perturbations. What exactly constitutes a *small* perturbation
    depends on the system. The applicability of tensor LEED is  normally
    limited to displacements of at most 0.2 Å (only in very simple cases
    up to 0.3 Å) :cite:`rousTensorLEEDTechnique1986`.
    For larger displacements, the optimization might give
    the right trends but the results should be taken with
    caution (see also :ref:`tensor_leed_errors`).

    When the trajectory of the structure optimization approaches the limit
    of applicability of tensor LEED, the range of the structural search can
    be extended by running new reference and delta-amplitudes calculations.

-   The parameter space grows quickly with the number of atoms in the unit
    cell. Luckily, many symmetries inherent to the surface structure can be
    exploited to eliminate redundant parameters. Specifically, displacements
    of symmetry-linked atoms must always happen in a concerted fashion. If
    that were not the case, the symmetry would be broken and a different
    LEED pattern would result.

    To make use of these symmetries and the resulting reduction of the
    parameter space, it is necessary to know and constrain the surface
    symmetry of the slab. While manually finding out the surface symmetry
    is normally an easy task, maintaining it is not. This is especially true
    for geometric displacements: It would require manually defining matching
    displacement vectors for all symmetry-linked atoms.

    *Fortunately for the user*, automatic symmetry detection and constraint
    is one of the **main features** of ViPErLEED.

-   When using Pendry's |R factor|, the |R factor| hypersurfaces are inherently
    rough :cite:p:`rousTensorLEEDApproximation1992`. Users should be aware that
    local minima are possible and that the optimization algorithm might
    get trapped in these minima if the parameter space is not opened up
    sufficiently. Simultaneous optimization over too many parameters is
    also a common cause of trapping in local minima.

-   As described above, the tensor-LEED implementation in TensErLEED separates
    the calculation of delta amplitudes and the structure optimization into two
    independent stages. As a direct consequence, the optimization can **only**
    be performed on a predefined grid of perturbations. Further, to achieve the
    best possible fit, the grid-based nature makes it necessary to run multiple
    sets of delta-amplitude calculations and structure optimizations with
    increasingly finer pitch.

    .. note::
        Starting with a fine grid over a large variation range
        is not recommended: Too many grid points per parameter
        significantly slow down the convergence.

-   The structure search implemented in TensErLEED has the additional
    limitation that geometric displacements are limited to one dimension
    per atom. During each search run, atoms can only be displaced along a
    predefined curve rather than freely in 3D space. To optimize the position
    of atoms in three dimensions, multiple sequential search runs are needed.
    See the entry on the :ref:`DISPLACEMENTS` file for details and workarounds
    (such as looping searches).


.. _tensor_leed_errors:

Tensor-LEED errors
==================

Since the tensor-LEED approach is based on first-order perturbation theory, it
is inherently limited to small perturbations.
The larger the perturbation, the larger the error incurred by the approximation
and the less reliable the result.

This should be kept in mind when interpreting the results of any ViPErLEED
segment that uses the tensor-LEED approach (i.e., the
:ref:`structure search<sec_search>` and the
:ref:`error calculation<error_calculation>`).
In particular, it is **strongly** recommended to run a new reference calculation
after the structure optimization has converged to get rid of any accumulated
errors. It is also usually necessary to iterate between structure search and
reference calculation to obtain the best possible fit.

There is one case, however, in which a full-dynamic calculation can yield more
erroneous results than tensor LEED. The full-dynamic reference calculation
cannot provide exact results when an atom has mixed chemical composition
and the elements have different optimized positions. This is because only
one position can be specified for each atom. In this case, the tensor-LEED
approximation is the only viable alternative. It should anyway be used with
care. In particular, the position deviations of the different chemical species
from the "mean" position used for the full-dynamic calculation should be small.
