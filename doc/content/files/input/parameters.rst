.. _parameters:

PARAMETERS
==========

The PARAMETERS file lists parameters defining what to run, how to 
interpret the other input files, and what to pass on to the TensErLEED 
scripts.
Generally, elements are separated by whitespace, and ``!``, ``#`` or 
``%`` mark the beginning of a comment.
Parameter names should be given at the beginning of a line, e.g.:

..  code-block:: none


   LAYER_CUTS = 0.09 0.19 0.29 0.39 0.49 
   N_BULK_LAYERS = 2

   ! a comment

   SITE_DEF Fe = siteA 38 41, siteB 50-55, topSite top(2)   !another comment
   SITE_DEF O = topSite top(2)

The order of parameters is not important, but defining a parameter 
twice will generally overwrite it.
An exception are parameters like ``SITE_DEF``, which can be defined for 
different elements.

.. only:: latex

.. toctree:: 
   
   params_by/param_name
   params_by/param_topics
   params_by/param_section