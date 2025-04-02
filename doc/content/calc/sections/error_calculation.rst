.. include:: /substitutions.rst


.. |varR|    replace:: :math:`\mathrm{var}(R_\mathrm{P})`
.. |R+varR|  replace:: :math:`R_{\mathrm{P,min}} + \mathrm{var}(R_\mathrm{P,min})`

.. _error_calculation:

Error calculation
=================

Once a best-fit structure has been determined, it is useful to see how strongly
small changes to specific parameters affect the |R factor|.
This can be obtained via a so-called "error calculation". In the error
calculation, displacements are applied to one parameter at a time, and the
|R factor| is determined for each step of the variation range, resulting
in an "error curve".

When using the :ref:`Pendry<pendry_r>` |R factor|, its standard error,
historically called |varR|,\ [1]_ provides an estimate of the uncertainty of
a given parameter: The point of intersection between the error curve for
one parameter and |R+varR| can be taken as a measure of the statistical
error of the parameter itself :cite:p:`heinzElectronBasedMethods2013`.
This gives the "error calculation" its name. |varR| is defined as

.. math::
    \mathrm{var}(R_\mathrm{P}) = R_\mathrm{P} \sqrt{ \frac{8 |V_{0\mathrm{i}}| }{ \delta E} },

where |RP| is the value of the Pendry |R factor|, |V0i| the imaginary part of
the inner potential, and :math:`\delta E` the total energy range of the beam
set (i.e., the sum of the energy ranges for all the symmetry-inequivalent
beams). For more information, see
Ref. :cite:alp:`pendryReliabilityFactorsLEED1980`.

.. note::
    Note that, when multiple parameters are varied separately, they are
    assumed to be statistically independent from one another. This is
    not always the case. For example, if a site can be occupied by
    different chemical elements, the site occupation (i.e., the element
    concentrations) and vibration amplitude can be strongly
    correlated :cite:p:`blumSegregationOrderingFe12001`.
    In such a case, the increase of the |R factor| when changing only one
    of these parameters does not provide a good estimate for the accuracy
    of that parameter.

Aside from the estimate of the accuracy of a displacement parameter, error
curves can be used to judge the impact of a certain displacement on the
overall |R factor|. While estimates on the parameter accuracy only apply
to the Pendry |R factor|, the impact of a parameter on the |R factor| is
also meaningful for other |R-factor| :ref:`definitions<r-factor_calculation>`.

-  You may find the |R factor| to be very insensitive to the displacement of
   some atoms (i.e., much less sensitive than for other atoms with a similar
   depth and similar scattering properties). This is an indication that the
   atom is absent, or that its position is **far away** from the correct one.
   In this case, consider increasing the displacement range for this atom. Note
   also that what is "far away" depends on how strongly the atom scatters
   (i.e., on the chemical species and on the depth within the solid), and
   in some cases may be as small as 0.1 Å, e.g., for displacements
   along the surface-normal, |z| direction of a surface atom.
-  Hydrogen is a very weak electron scatterer: the |R factor| depends only
   weakly on its position.

.. attention::
  The |R factor| values obtained in the error calculation contain
  errors from the :ref:`tensor-LEED approximation<tensor_leed_errors>`.


Error calculation in ViPErLEED
------------------------------

.. todo::
    This part is a bit unclear as to what one can actually do. E.g.,
    can one CONSTRAIN to vary stuff simultaneously? A few example
    inputs for DISPLACEMENTS in non-trivial cases would help.

To run the error calculation, set :ref:`RUN = 5<RUN>` in the :ref:`PARAMETERS`
file. It is recommended to first run a reference calculation (e.g., via
``RUN = 1 5``), as the error curves may not be centered otherwise.

Aside from the :ref:`PARAMETERS` file, the required input is the same as for
running a :ref:`sec_deltas` and a :ref:`sec_search`: structural input files
(\ :ref:`POSCAR`, and, optionally, :ref:`VIBROCC`),
:ref:`experimental beams<EXPBEAMS>`, a set of :ref:`Tensors<Tensorszip>`
from a :ref:`reference calculation<ref-calc>`, and a :ref:`DISPLACEMENTS`
file defining which parameters should be varied. If multiple parameters
are linked (e.g., by symmetry), they are treated as one and varied together.

Note that defining multiple sections
in the :ref:`DISPLACEMENTS` file, as is possible for the search, is not
allowed here. Only one section of the :ref:`DISPLACEMENTS` file is read: the
last one if the error calculation is run following a :ref:`sec_search`, the
first one otherwise. Defining geometric, vibration, and occupation variations
all in the same :ref:`DISPLACEMENTS` file is allowed, but the different
variations are split up, so the result is the same as executing multiple
error calculations. This means that you cannot have simultaneous error
calculations for multiple geometric-displacement directions (e.g.,
|x| *and* |z|) for the same atom, since this would require multiple
consecutive blocks in the DISPLACEMENTS file.

.. todo::
    The next note used to say "Simultaneous geometric AND vibration..."
    but was in contraction with the "we read one block", or was unclear
    how the input would be. This is a typical case in which an explicit
    example would be useful.

.. tip::
    Simultaneous geometric or vibration variations of multiple chemical
    elements occupying the same site is possible, and the displacement values
    for the different elements may differ. However, all displacement ranges
    must have the same number of steps.

.. (unless some of the displacements are :ref:`explicitly linked<searchconstraints>`).  IS THIS TRUE? If yes, it should go between "calculations." and "This means that you cannot"

The error calculation does *not* require a set of
:ref:`Delta files<Deltaszip>`, since the normal delta-calculation routines
mix geometric and vibration displacements. Instead, the error calculation
runs the required delta calculations automatically, splitting the
geometric and vibration variations into separate delta files to
reduce computational cost.

The results of the error calculation consists of the :ref:`errorspdf_header`
files, stored in the :file:`OUT` directory.

If the Pendry |R factor| is used, the value of |R+varR| is calculated
for each error type and drawn as a horizontal line in :ref:`errorspdf`.
If  ViPErLEED finds an intersection between the error curve for any
parameter :math:`p` and the line |R+varR|, the intersection points are
used to calculate the (estimate of the) statistical uncertainty for
:math:`p`. This is written to :ref:`errorscsv` and plotted in :ref:`errorspdf`.


.. [1] It is a common mistake to refer to |varR| as the "variance" of the
       |R factor|. This is **wrong**, as, from a dimensional point of view,
       it does not correspond to a variance (it is proportional to |RP|,
       not to its square). If it were a variance, it would not make sense
       to calculate |R+varR|.
       The original "\ :math:`\mathrm{var}`" abbreviation was
       probably intended as a contraction for "variation"
       :cite:p:`pendryReliabilityFactorsLEED1980`.
