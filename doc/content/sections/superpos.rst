.. include:: /substitutions.rst

.. _super_pos:

=========================
Superposition calculation
=========================

The superposition calculation (``superpos`` in the code) is the last section
of "classic" TensErLEED and is automatically executed by |calc| after the
:ref:`structure search<sec_search>`.

It takes the :ref:`delta files<deltaszip>` and the result of the structure
search – which at that point is only a set of optimized parameters – and
performs the name-giving superposition, i.e., summing up of delta-amplitudes
to a final set of beam-amplitudes and intensities. Then, it generates an
output of theoretical beams in the same format as  yielded by the
:ref:`refercence Calculation<ref-calc>` (see :ref:`file FITBEAMS<fitbeams>`).
Finally, an |R factor| is calculated based on this file and the files
:ref:`Rfactor_plots_superpos.pdf<rfactorplots>` and
:ref:`Rfactor_analysis_superpos.pdf<rfactoranalysis>` are produced
(see the corresponding sections for details/examples).

.. note::
    It is recommended to run another reference calculation for the final
    structure as the :ref:`FITBEAMS file<fitbeams>` will contain errors
    due to the :ref:`tensor LEED approximation<tensor_leed>`.
