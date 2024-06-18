.. include:: /substitutions.rst

.. |t-LEED| replace:: :ref:`the tensor-LEED approximation<tensor_leed>`

.. _work-segments:

ViPErLEED work segments
-----------------------

ViPErLEED operates using a set of self-contained work segments (see also the
:ref:`RUN` parameter). The three main segments, following the logic of
calculations using |t-LEED|, are:

#. :ref:`ref-calc`: Full-dynamic LEED calculation, which outputs a set
   of :ref:`theoretical beams<THEOBEAMS>` for a given structure and
   :ref:`the "Tensors"<Tensorszip>`.
#. :ref:`sec_deltas`:
   The delta amplitudes specify how parameter changes affect the scattering
   amplitudes within |t-LEED|. This calculation is based on the
   :ref:`tensors<Tensorszip>` and a set of parameter variations
   :ref:`specified by the user<DISPLACEMENTS>`. The output of a
   delta-amplitudes calculation are the :ref:`"Delta files"<Deltaszip>`.
#. :ref:`sec_search`: Using the :ref:`Delta files<Deltaszip>` to vary
   the theoretical beams, looks for a set of parameters such that the
   :math:`R` :ref:`factor<r-factor_calculation>` between the theoretical
   beams and :ref:`a given set of experimental beams<EXPBEAMS>` is minimized.

Which of these segments should be executed must be specified using the
:ref:`RUN` parameter, using the segment numbers in the list above or a
contraction of their names. More information on the allowed contractions
are found in the documentation for :ref:`RUN`.

Besides these three main segments, there are also the following minor segments,
which are inserted automatically during normal ViPErLEED execution when
appropriate (but can also be explicitly selected via :ref:`RUN`):

-  :ref:`initialization`: Always runs at the beginning. Reads and checks input
   files, runs symmetry search, generates derived input files if appropriate.
-  :ref:`super_pos`: Automatically runs after the :ref:`search<sec_search>`.
   Generates a set of theoretical beams for the actual best-fit configuration
   based on |t-LEED|.
-  |R-factor| :ref:`calculation<r-factor_calculation>`:
   Automatically runs after the :ref:`reference-calculation segment<ref-calc>`
   if an :ref:`experimental-beams file<EXPBEAMS>` is present, and after the
   :ref:`superpos section<super_pos>`. Calculates the |R factor| per beam and
   for the entire set of beams, and outputs an
   :file:`Rfactor_plots_<section>.pdf` :ref:`file<Rfactorplots>`.

Further specialized segments include:

-  :ref:`error_calculation`: Based on a given reference structure (i.e., after
   a :ref:`reference calculation<ref-calc>`), calculates one-dimensional error
   curves for variation of a single parameter. Effectively, this produces
   :ref:`delta amplitudes<sec_deltas>` for variations of a single parameter,
   and outputs the |R factor| for every single configuration along that axis.
-  :ref:`fdoptimization`: Optimizes parameters that are not accessible to
   |t-LEED|, like :ref:`BEAMINCIDENCE`, :ref:`V0_IMAG`, or unit-cell scaling. 
   This is achieved by performing multiple full-dynamic (i.e., "reference") 
   calculations (but without producing Tensor files). The behavior is 
   controlled by the :ref:`OPTIMIZE` parameter.

The pages listed above cover normal operation, in which the theoretical beams
correspond to only one surface structure. If multiple structures coexist on
the sample, the same segments need to be executed, but their behavior is
somewhat different, as described in

-  :ref:`Domain calculations<domain_calculation>`:
   :ref:`Reference calculations<ref-calc>` are run separately for the different
   domains (if necessary) and :ref:`delta amplitudes<sec_deltas>` are generated
   independently. The :ref:`search<sec_search>` then combines the optimization
   of the different structures — weighted by their area fraction — for the best
   overall |R factor| with respect to the
   :ref:`experimental beam set<EXPBEAMS>`.


.. toctree::
   :maxdepth: 1
   :caption: Sections

   sections/initialization
   sections/ref-calc
   Delta amplitudes<sections/deltas>
   Search<sections/structure_search>
   Superpos<sections/superpos>
   sections/domain_calculation
   sections/error_calculation
   sections/fdoptimization
   R factors<sections/r-factor_calculation>
   sections/bookkeeper
