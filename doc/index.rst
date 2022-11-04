.. _index:

ViPErLEED documentation
=======================

Welcome to the manual for :term:`ViPErLEED` (Vienna Package for Erlangen LEED) 
and the :term:`Python` package :term:`tleedm`. See the :ref:`Getting Started page<getting_started>`.

The ViPErLEED (Vienna Package for Erlangen LEED) project is a set of open-source tools that aims at drastically reducing the effort for LEED :math:`I(V)` studies, both on the computational and on the experimental side.

The package consists of 
    i.  software for calculation of :math:`I(V)` curves for a given structure 
        and structure optimization, by minimizing the difference between 
        the calculated and experimental :math:`I(V)` data, and, on the 
        experimental side, 
    #.  hardware and software for data acquisition, as well as 
    #.  software for extracting :math:`I(V)` curves from the experimental data. 

For details please see the ViPErLEED publication series (**TODO**). This manual primarily deals with part **i.**.

The goal of any :term:`LEED-I(V)` calculation is the calculation of 
energy dependent electron scattering amplitudes and intensities. These 
intensity curves (often referred to as :math:`I(V)` curves or spectra) 
are very sensitive to the precise position and vibrational amplitudes in
the surface unit cell.
For more details consult works that cover the basics of :term:`LEED` and 
:term:`LEED-I(V)`, e.g.
Chapter 4 in :cite:t:`fausterSurfacePhysicsFundamentals2020` or 
the overview by:cite:t:`heinzElectronBasedMethods2013`.
In ViPErLEED, these calculations are performed by the TensErLEED manager :term:`tleedm`,
a Python package that is built based on, and as a comprehensive feature extension to
:term:`TensErLEED`.

TensErLEED is used for the calculation of diffraction intensities of 
surface slabs (see also :ref:`reference calculation<ref-calc>`) 
and structure optimization using the :ref:`Tensor LEED approach<tensor_leed>`.
For computational details please have a look at the ViPErLEED paper 
(**TODO**) and the original work describing TensErLEED by Blum and Heinz 
:cite:p:`blumFastLEEDIntensity2001a`.


.. Table of contents in LaTeX pdf called Contents
.. only:: latex

   .. toctree::
      :caption: Contents

.. toctree:: 

   Getting started<content/getting_started>
   ViPErLEED segments<content/work_segments>
   ASE Interface<content/aseapi>

.. toctree:: 
   :maxdepth: 2
   :caption: Files

   content/input_files
   content/output_files
   content/supp_files


.. only:: html

   .. toctree:: 
      :maxdepth: 1
      :caption: Parameters

      Overview by Name<content/files/input/params_by/param_name>
      Overview by Section<content/files/input/params_by/param_section>
      Overview by Topic<content/files/input/params_by/param_topics>


.. toctree:: 
   :hidden:

   content/param_toc


.. toctree:: 
   :maxdepth: 1
   :caption: Utilities

   Utilities<content/utilities>



.. toctree:: 
   :maxdepth: 1
   :caption: References

   References<references>
   Glossary<content/glossary>