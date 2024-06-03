.. _index:

.. raw:: latex

   % Modify color of links to conform with those in html
   \def\Hy@colorlink#1{\begingroup\fontshape{\default}\selectfont}%
   \begingroup
   \sphinxsetup{%
      InnerLinkColor={RGB}{41,128,185},
      OuterLinkColor={RGB}{41,128,185}
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

.. include:: /substitutions.rst

=======================
ViPErLEED documentation
=======================

Welcome to the documentation for :term:`ViPErLEED` and for the ``viperleed``
:term:`Python` package.
The ViPErLEED project is a set of :ref:`open-source<license>` tools that aims
at drastically reducing the effort for an intensity [i.e., |LEED-IV|] analysis,
both on the computational and on the experimental side.

.. todo::
    See also the ViPErLEED publication series (LINKS).

:numref:`toc_figure` shows an overview of the contents of the ViPErLEED package.
It consists of:

-  :ref:`The ViPErLEED hardware<hardware>`
      A set of hardware, firmware and control software for the easy
      acquisition of |LEED-IV| data with pre-existing LEED systems.
- :ref:`The ViPErLEED Spot Tracker<spot_tracker>`
      Software for extracting :math:`I(V)` curves from the experimental data
      ("movies").
- :ref:`The viperleed.calc package<viperleed_calc>`
      A Python package for the calculation of :math:`I(V)` curves,
      quantitative analysis of :term:`LEED` data, and surface structure
      optimization.
      See the :ref:`Getting Started page<getting_started>`.

.. _toc_figure:
.. figure:: /_static/paper_figures/ViPErLEED-overview_embedded.svg
   :width: 70%
   :align: center

   Overview of the parts of the ViPErLEED project.


.. Table of contents in LaTeX pdf called Contents

.. only:: latex

   .. toctree::
      :caption: Contents

.. toctree::
   :caption: viperleed.calc

   viperleed.calc<content/viperleed_calc>
   Examples<content/examples>
   Getting started<content/getting_started>
   Segments<content/work_segments>
   Files<content/files>
   Parameters<content/param_topics>
   Utilities<content/utilities>
   API<content/api>

.. toctree::
   :maxdepth: 1
   :caption: Spot tracker

   Spot-Tracker<content/spot_tracker>

.. toctree::
   :maxdepth: 1
   :caption: Hardware and measurements

   Hardware<content/hardware>

.. toctree::
   :maxdepth: 1
   :caption: References

   Citing ViPErLEED<content/citing>
   References<references>
   Glossary<content/glossary>
   License<content/license>

.. raw:: latex

   \endgroup
