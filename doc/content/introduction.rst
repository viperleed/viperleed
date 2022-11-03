.. _introduction:

Introduction
============

Welcome to the manual for :term:`ViPErLEED`` (Vienna Package for Erlangen LEED) 
and the :term:`Python` package :term:`tleedm`.

To get started see:

- Installation
- Getting started
- Basics and conventions
- ViPErLEED sections
- Examples
- ASE interface


The ViPErLEED (Vienna Package for Erlangen LEED) project is a set of 
open-source tools that aims at drastically reducing the effort for LEED $I(V)$ studies, both on the computational and on the experimental side.

The package consists of 
    i.  software for calculation of $I(V)$ curves for a given structure 
        and structure optimization, by minimizing the difference between 
        the calculated and experimental $I(V)$ data, and, on the 
        experimental side, 
    #.  hardware and software for data acquisition, as well as 
    #.  software for extracting $I(V)$ curves from the experimental data. 

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


