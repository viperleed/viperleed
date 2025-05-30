.. include:: /substitutions.rst

.. _super_pos:

=========================
Superposition calculation
=========================

The superposition calculation (``superpos`` in the code) is the last section
of "classic" TensErLEED and is automatically executed by |calc| after the
:ref:`structure search<sec_search>`.

It takes the :ref:`delta files<deltaszip>` and the result of the structure
search — which at that point is only a set of optimized parameters — and
performs the name-giving superposition: it adds the sum of delta amplitudes
corresponding to the optimal parameters to the reference amplitudes, producing
a final set of total beam amplitudes and intensities in the same format as
yielded by the :ref:`ref-calc`.

ViPErLEED processes the intensities to the :ref:`FITBEAMS` file. It then
calculates the |R factor| between these beams and those in the :ref:`EXPBEAMS`
file, producing the :ref:`Rfactor_plots_superpos.pdf<rfactorplots>`
and :ref:`Rfactor_analysis_superpos.pdf<rfactoranalysis>` files
(see the corresponding sections for details and examples).

.. note::
    It is recommended to run another reference calculation for the 
    final structure, as the :ref:`FITBEAMS` contain errors due to the 
    :ref:`tensor-LEED approximation<tensor_leed_errors>`.
