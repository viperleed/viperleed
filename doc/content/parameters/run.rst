.. include:: /substitutions.rst

.. _run:

===
RUN
===

RUN defines what parts of the TensErLEED package should be run.

**Default:** 0-3

**Allowed values:** list of integers and/or integer ranges, with integers
0-3, or single integers up to 6. Instead of the integers, also the aliases
in the table can be used:

+-----+------------------------------------------------------------------------------------------+-----------------------------------------------+
| RUN | Action                                                                                   | Aliases                                       |
+=====+==========================================================================================+===============================================+
| 0   | Initialization                                                                           | ``initialization, initialisation, init, ini`` |
+-----+------------------------------------------------------------------------------------------+-----------------------------------------------+
| 1   | Reference calculation                                                                    | ``refcalc, ref, fd``                          |
+-----+------------------------------------------------------------------------------------------+-----------------------------------------------+
| 11  | |R-Factor| calculation (automatically inserted after reference calculation,              | ``rfactor, rfac``                             |
|     | if experimental beam file is found)                                                      |                                               |
+-----+------------------------------------------------------------------------------------------+-----------------------------------------------+
| 12  | |R-Factor| calculation based on :ref:`Superpos<super_pos>` results                       | ``rfacsuper``                                 |
|     | (automatically inserted after superpos calculation,                                      |                                               |
|     | if experimental beam file is found)                                                      |                                               |
+-----+------------------------------------------------------------------------------------------+-----------------------------------------------+
| 2   | DeltaAmplitudes                                                                          | ``deltaamplitudes, deltas, delta, del``       |
+-----+------------------------------------------------------------------------------------------+-----------------------------------------------+
| 3   | Search                                                                                   | ``search``                                    |
+-----+------------------------------------------------------------------------------------------+-----------------------------------------------+
| 31  | :ref:`Superpos<super_pos>`  calculation (automatically inserted after search)            | ``superpos, super, sup``                      |
+-----+------------------------------------------------------------------------------------------+-----------------------------------------------+
| 4   | :ref:`Domain calculation<domain_calculation>`                                            | ``domains, domain, dom``                      |
+-----+------------------------------------------------------------------------------------------+-----------------------------------------------+
| 5   | :ref:`Error calculation<error_calculation>`                                              | ``error, err``                                |
+-----+------------------------------------------------------------------------------------------+-----------------------------------------------+
| 6   | :ref:`Full-dynamic optimization<Fdoptimization>`                                         | ``opt, optimize, fdopt``                      |
+-----+------------------------------------------------------------------------------------------+-----------------------------------------------+

Syntax
------

::

   RUN = 0-3 1
   RUN = 0
   RUN = ini     ! same as RUN = 0
   RUN = 1       ! will automatically be changed to '0 1'
   RUN = 0 3 1   ! runs initialization, search, and a final reference calculation. Requires pre-existing deltas.
   RUN = ref-search ref   ! same as 1-3 1

Initialization will always be run, even if the leading 0 is not there. Indices
1–3 are meant to be performed in that order, and any other (without data from
a previous run) will fail.

Index 4 is a helper index for :ref:`domain calculations<domain_calculation>`,
and should only be used together with the :ref:`DOMAIN`  parameter.
If index 4 is specified in RUN, the input files for the different domains
will be checked for compatibility; if Tensor files exist for all domains
and are compatible, the calculation will proceed straight to DeltaAmplitudes
and Search. Otherwise, reference calculations will be run as needed for the
domains. Note that you can also specify domains by the :ref:`DOMAIN` parameter
and use RUN as normal, without setting it to 4. In that case, the program will
try to execute the specified segments, and stop if that's not possible (e.g.,
if no reference calculation was specified, but the Tensors are missing or
incompatible).
