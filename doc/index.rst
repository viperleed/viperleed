.. _index:

.. raw:: latex

   \def\Hy@colorlink#1{\begingroup\fontshape{\default}\selectfont}%
   \begingroup
   \sphinxsetup{%
      InnerLinkColor={rgb}{0,0.120,0.204},
      OuterLinkColor={rgb}{0,0.120,0.204}

   }

   % use heavy boxes for note, hint, important & tip
   \renewenvironment{sphinxnote}[1]
      {\begin{sphinxheavybox}\sphinxstrong{#1} }{\end{sphinxheavybox}}
   \renewenvironment{sphinxhint}[1]
      {\begin{sphinxheavybox}\sphinxstrong{#1} }{\end{sphinxheavybox}}
   \renewenvironment{sphinximportant}[1]
      {\begin{sphinxheavybox}\sphinxstrong{#1} }{\end{sphinxheavybox}}
   \renewenvironment{sphinxtip}[1]
      {\begin{sphinxheavybox}\sphinxstrong{#1} }{\end{sphinxheavybox}}

   \AddEverypageHook{
        \settowidth{\chapterNameLength}{\leftmark}%
        \begin{textblock}{1}(0,0)%first argument {1} is number of blocks horiz
        \,\ \hyperlink{link_content}{$\rightarrow$Contents}%
        \,\ \ \ \Acrobatmenu{GoBack}{$\leftarrow$Back}%
        \,\ \Acrobatmenu{GoForward}{Forward$\rightarrow$}%
        \end{textblock}%
    }%end AddEverypageHook

ViPErLEED documentation
=======================



Welcome to the manual for :term:`ViPErLEED` and the :term:`Python` package :term:`tleedm`.
See the :ref:`Getting Started page<getting_started>`.

The ViPErLEED project is a set of :ref:`open-source<license>` tools that aims at drastically reducing the effort for LEED-:math:`I(V)` studies, both on the computational and on the experimental side.

The package consists of:


    i.  Hardware and software for data acquisition.
    #.  Software for extracting :math:`I(V)` curves from the experimental data ("movies").
    #.  Software for calculation of :math:`I(V)` curves for a given structure and structure optimization, by minimizing the difference between the calculated and experimental :math:`I(V)` data.

For details please see the ViPErLEED publication series (**TODO**).
This manual primarily deals with part **iii**.

The goal of any :term:`LEED-I(V)` calculation is the calculation of energy-dependent electron-scattering amplitudes and intensities.
These intensity curves [often referred to as :math:`I(V)` curves or spectra]  are very sensitive to the precise position and vibrational amplitudes of each atom in the surface unit cell.
For more details consult works that cover the basics of :term:`LEED` and :term:`LEED-I(V)`, e.g.
Chapter 4 in :cite:t:`fausterSurfacePhysicsFundamentals2020,fausterOberflachenphysikGrundlagenUnd2013` or the overview by :cite:t:`heinzElectronBasedMethods2013`.
In ViPErLEED, these calculations are performed by the TensErLEED manager :term:`tleedm`,
a Python package that is based on, and as a comprehensive feature extension to :term:`TensErLEED`.

TensErLEED is used for the calculation of diffraction intensities of surface slabs (see also :ref:`reference calculation<ref-calc>`) and structure optimization using the :ref:`tensor LEED approach<tensor_leed>`.
For computational details please have a look at the ViPErLEED paper (**TODO**) and the original work describing TensErLEED by Blum and Heinz :cite:p:`blumFastLEEDIntensity2001a`.


.. Table of contents in LaTeX pdf called Contents

.. only:: latex

   .. toctree::
      :caption: Contents

.. raw:: latex

   \hypertarget{link_content}{}

.. toctree:: 

   Getting started<content/getting_started>
   ViPErLEED segments<content/work_segments>
   ASE Interface<content/aseapi>

.. toctree:: 
   :maxdepth: 2
   :caption: Files

   Input files<content/input_files>
   Output files<content/output_files>
   Supplementary files<content/supp_files>


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
   License<content/license>

.. raw:: latex

   \endgroup
