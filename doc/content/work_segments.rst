.. include:: /substitutions.rst

.. _work-segments:

ViPErLEED work segments
-----------------------

ViPErLEED operates using a set of self-contained work segments (see also the
:ref:`RUN<RUN>`  parameter). The three main segments, following the logic of
calculations using the tensor-LEED approximation, are:

#. :ref:`Reference calculation<ref-calc>`: Full-dynamic LEED calculation,
   which outputs a set of :ref:`theoretical beams<THEOBEAMS>` for a
   given structure and :ref:`the "Tensors"<Tensorszip>`
#. :ref:`Delta-amplitudes generation<sec_deltas>`:
   The delta-amplitudes specify how parameter changes affect the scattering
   amplitudes within the :ref:`tensor-LEED approximation<tensor_leed>`. This
   calculation is based on the tensors and a
   :ref:`set of parameter variations specified by the user<DISPLACEMENTS>`,
   produces :ref:`"Delta files"<Deltaszip>`.
#. :ref:`Search<sec_search>`: Using the :ref:`Delta files<Deltaszip>` to vary
   the theoretical beams, looks for a set of parameters such that the :math:`R`
   :ref:`factor<r-factor_calculation>` between the theoretical beams
   and :ref:`a given set of experimental beams<EXPBEAMS>` is minimized.

Which of these segments should be executed must be specified using the
:ref:`RUN<RUN>`  parameter, using the segment numbers in the list above.
Besides these main three segments, there are also the following minor
segments, which will be inserted automatically during normal ViPErLEED
execution when appropriate:

-  :ref:`Initialization<initialization>`: Always runs at the beginning; reads
   and checks input files, runs symmetry search, generates derived input files
   if appropriate.
-  :ref:`Superpos calculation<super_pos>`: Automatically runs after the search.
   Generates a set of theoretical beams for the actual best fit configuration
   based on the tensor-LEED approximation,
-  |R-factor| :ref:`calculation<r-factor_calculation>`: Automatically runs
   after the :ref:`reference calculation<ref-calc>` and superpos segments,
   if an :ref:`experimental beams file<EXPBEAMS>` is present. Calculates the
   |R factor| per beam and for the entire set of beams, and outputs an
   :ref:`Rfactor_plots pdf file<Rfactorplots>`.

Further specialized segments include:

-  :ref:`Error calculations<error_calculation>`: Based on a given reference
   structure (i.e., after a reference calculation has been run), calculate
   one-dimensional error curves for variation of a single parameter.
   Effectively, this calculates delta amplitudes for variations of a single
   parameter, and outputs the |R factor| for every single configuration
   along that axis.
-  :ref:`Full-dynamic optimization<fdoptimization>`: Optimize parameters that
   cannot be varied during the search, like
   :ref:`BEAM_INCIDENCE<BEAMINCIDENCE>`, :ref:`V0_IMAG<v0_imag>`
   or unit-cell scaling. This is achieved by performing multiple
   full-dynamic (i.e., "reference") calculations (without Tensor output).
   The behavior is controlled by the :ref:`OPTIMIZE<OPTIMIZE>` parameter.

The pages listed above cover normal operation, in which the theoretical beams
correspond to only one surface structure. If multiple structures coexist on
the sample, the same segments need to be executed, but their behavior is
somewhat different, as described in

-  :ref:`Domain calculations<domain_calculation>`: Reference calculations
   are run separately for the different domains (if necessary) and
   Delta-amplitudes are generated independently.

The search then combines the optimization of the different structures –
weighted by their area fraction – for the best overall |R factor|, compared
to the experimental beam set.


.. If we want to, we could hide the below table in the HTML version by
   creating ..only html & ..only latex versions.

.. toctree::
   :maxdepth: 1
   :caption: Sections

   Initialization<sections/initialization>
   Reference calculation<sections/ref-calc>
   Delta-amplitudes<sections/deltas>
   Search<sections/structure_search>
   Superpos<sections/superpos>
   Domain calculations<sections/domain_calculation>
   Error calculations<sections/error_calculation>
   Full-dynamic optimization<sections/fdoptimization>
   R factors<sections/r-factor_calculation>
   Bookkeeper<sections/bookkeeper>
