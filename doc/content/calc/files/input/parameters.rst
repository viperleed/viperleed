.. include:: /substitutions.rst

.. _parameters:

PARAMETERS
==========

The |PARAMETERS| file lists parameters defining what to run, how to
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
An exception are parameters like :ref:`SITEDEF`, which can be defined for
different elements.

.. toctree::
    :maxdepth: 1

    ../../param_name
    ../../param_topics
    ../../param_section

The |PARAMETERS| file given as input may be automatically edited by |calc|.
At the end of each |calc| execution, the file given as input for that run
is renamed to :file:`PARAMETERS_ori`, while the (potentially) edited file
is copied to the root directory (from |OUT|) as a new |PARAMETERS| file.
This ensures that further invocations of |calc| will automatically use
the output of previous executions as an input. You can manually call the
|bookkeeper| utility after a specific |calc| run if this behavior is not
desirable. See the :ref:`bookkeeper` page for more details.

.. versionchanged:: 0.13.0
    In earlier versions of |calc| the |PARAMETERS| files given as inputs
    were renamed to :file:`PARAMETERS_ori_<timestamp>` before editing.
    The ``_<timestamp>`` was dropped in ``v0.13.0``. The original,
    unedited version is stored in :file:`SUPP/original_inputs`.
    The edited |PARAMETERS| file can be found in |OUT|.
