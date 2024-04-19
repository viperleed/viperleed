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

   % renew table of contents to include link
   \renewcommand{\sphinxtableofcontents}{%
   \pagenumbering{roman}%
   \begingroup
      \parskip \z@skip
      \sphinxtableofcontentshook
      \tableofcontents
      \hypertarget{link_content}{}
   \endgroup
   % before resetting page counter, let's do the right thing.
   \if@openright\cleardoublepage\else\clearpage\fi
   \pagenumbering{arabic}%
   }

=======================
ViPErLEED documentation
=======================

Welcome to the documentation for :term:`ViPErLEED` and the :term:`Python`
package :term:`tleedm`.
The ViPErLEED project is a set of :ref:`open-source<license>` tools that aims at
drastically reducing the effort for a intensity [LEED-:math:`I(V)`] analysis,
both on the computational and on the experimental side.

.. todo:: Add: See also the ViPErLEED publication series (**TODO: add link **).

The ViPErLEED package consists of:

- :ref:`The viperleed calc package<viperleed_calc>`
      A Python package for the calculation of :math:`I(V)` curves,
      quantitative analysis of :term:`LEED` data, and surface structure 
      optimization.
      See the :ref:`Getting Started page<getting_started>`.
- :ref:`The ViPErLEED Spot-Tracker<spot_tracker>`
      Software for extracting :math:`I(V)` curves from the experimental data
      ("movies").
-  The ViPErLEED hardware
      A set of hardware, firmware and control software for the easy
      acquisition of LEED-:math:`I(V)` data with pre-existing LEED systems.

.. _toc_figure:
.. figure:: /_static/paper_figures/ViPErLEED-overview.svg
   :width: 70%
   :align: center

   Overview of the parts of the ViPErLEED project.


.. Table of contents in LaTeX pdf called Contents

.. only:: latex

   .. toctree::
      :caption: Contents

.. toctree:: 
   :caption: viperleed calc

   viperleed calc<content/viperleed_calc>
   Getting started<content/getting_started>
   Examples<content/examples>
   ViPErLEED segments<content/work_segments>
   Files<content/files>

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

   Utilities<content/utilities>
   API<content/api>

.. toctree:: 
   :maxdepth: 1
   :caption: Spot-Tracker

   Spot-Tracker<content/spot_tracker>

.. toctree:: 
   :maxdepth: 1
   :caption: References

   Citing ViPErLEED<content/citing>
   References<references>
   Glossary<content/glossary>
   License<content/license>

.. raw:: latex

   \endgroup
