.. _index:

ViPErLEED documentation
=======================

The Vienna Package for TensErLEED (ViPErLEED) is a collection of tools for modern LEED-I(V).

*TODO*: add nice introduction

.. Table of contents in LaTeX pdf called Contents
.. only:: latex

   .. toctree::
      :caption: Contents

.. toctree:: 
   :maxdepth: 2
   :caption: Segments

   content/work_segments


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

   Atomic Simulation Environment API<content/utilities/aseapi>
   Bookkeeper<content/utilities/bookkeeper>
   Further Utilities<content/utilities>

.. toctree:: 
   :glob:
   :caption: Various
   :maxdepth: 1

   content/various/*



.. toctree:: 
   :maxdepth: 1
   :caption: References

   References<references>
   Glossary<content/glossary>