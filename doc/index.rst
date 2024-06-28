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

Welcome to the documentation of :term:`ViPErLEED` and of the ``viperleed``
:term:`Python` package.
The ViPErLEED project comprises a set of :ref:`open-source<license>` tools that
aims at drastically reducing the effort for quantitative low-energy electron
diffraction [i.e., |LEED-IV|] on both experimental and computational fronts.


The ViPErLEED source code is hosted on
GitHub at `<https://github.com/viperleed>`__.

.. todo::
    See also the ViPErLEED publication series (LINKS).

:numref:`toc_figure` shows an overview of the tools provided by the ViPErLEED
project:

-  The :ref:`ViPErLEED hardware<hardware>`
      A set of hardware, firmware and control software for the easy
      acquisition of |LEED-IV| data with pre-existing LEED systems.
-  The :ref:`imagej_plugins`
      Software for extracting |IV| curves from the experimental data
      ("movies").
-  The :ref:`viperleed.calc package<viperleed_calc>`
      A Python package for the calculation of |IV| curves, quantitative
      analysis of :term:`LEED` data, and surface-structure optimization.

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
    :maxdepth: 2
    :caption: Getting started

    content/installation
    content/command_line_tools
    content/background

.. toctree::
    :maxdepth: 2
    :caption: viperleed.calc

    viperleed calc<content/viperleed_calc>

.. toctree::
    :maxdepth: 2
    :caption: ViPErLEED ImageJ plugins

    content/imagej_plugins

.. toctree::
    :maxdepth: 2
    :caption: Hardware and measurements

    Hardware<content/hardware>

.. toctree::
    :maxdepth: 1
    :caption: ViPErLEED Python API

    content/api


..
    HANDLING OF APPENDICES
    The next part is a bit ugly. Its purpose is to make only one \part
    in LaTeX named "Appendix" (see content/appendix.rst) and have its
    contents be "chapters". Using for latex the same definition as
    for the non-latex case would make each chapter of "Appendix" a
    part by itself, which is ugly. On the other hand, we really want
    the contents of the appendix to appear at the top level in the
    navigation bar in HTML.


.. only:: latex

    .. toctree::

        content/appendix

.. only:: not latex

    .. toctree::
        :maxdepth: 1
        :caption: Appendix

        content/appendix/citing
        content/appendix/references
        content/appendix/glossary
        License<content/appendix/license>
        content/appendix/notes_for_developers


.. raw:: latex

   \endgroup
