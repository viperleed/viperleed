.. include:: /substitutions.rst

.. _run:

===
RUN
===

RUN defines what parts of the TensErLEED package should be run.

**Default:** 0-3

**Allowed values:** list of integers and/or integer ranges, with integers
0-3, or single integers up to 6. Instead of the integers, also the section
names in the table, or their aliases, can be used. Aliases are not case
sensitive, and the underscores in them are optional (e.g., ``fd_opt`` and
``fdopt`` are equivalent). Longer and shorter forms of the same words are
also accepted (e.g., ``calc`` or ``calculation``, ``optimize`` or
``optimise``). Names that start with an alias are also accepted (e.g.,
``superposition``), and are interpreted as the section with the longest
such alias, provided they consist only of letters, digits and
underscores. Spaces cannot be used inside a name, as they separate items
in RUN. Hyphens and colons separate ranges. A hyphen can be used inside
a name only where it replaces the underscore before ``calc``,
``calculation``, ``amplitudes``, ``opt`` or ``optimization`` (e.g.,
``delta-amplitudes``, ``ref-calc`` or ``fd-opt``). Any other name with
a hyphen or a colon is interpreted as a range.

+-----+------------------------------------------------------------------------------------------+----------------------------------------------------+
| RUN | Action                                                                                   | Aliases (examples)                                 |
+=====+==========================================================================================+====================================================+
| 0   | Initialization                                                                           | ``init``, ``ini``, ``initialize``,                 |
|     |                                                                                          | ``initialization``                                 |
+-----+------------------------------------------------------------------------------------------+----------------------------------------------------+
| 1   | Reference calculation                                                                    | ``ref``, ``fd``, ``refcalc``,                      |
|     |                                                                                          | ``reference_calculation``, ``full_dynamic``,       |
|     |                                                                                          | ``fd_calc``                                        |
+-----+------------------------------------------------------------------------------------------+----------------------------------------------------+
| 11  | |R-Factor| calculation (automatically inserted after reference calculation,              | ``rfac``, ``rfactor``, ``r_factor``                |
|     | if experimental beam file is found)                                                      |                                                    |
+-----+------------------------------------------------------------------------------------------+----------------------------------------------------+
| 12  | |R-Factor| calculation based on :ref:`Superpos<super_pos>` results                       | ``rfacsuper``, ``rfac_super``,                     |
|     | (automatically inserted after superpos calculation,                                      | ``rfactor_superpos``                               |
|     | if experimental beam file is found)                                                      |                                                    |
+-----+------------------------------------------------------------------------------------------+----------------------------------------------------+
| 2   | DeltaAmplitudes                                                                          | ``del``, ``delta``, ``deltas``,                    |
|     |                                                                                          | ``delta_calculation``, ``deltaamplitudes``         |
+-----+------------------------------------------------------------------------------------------+----------------------------------------------------+
| 3   | Search                                                                                   | ``search``, ``structure_search``,                  |
|     |                                                                                          | ``struct_search``                                  |
+-----+------------------------------------------------------------------------------------------+----------------------------------------------------+
| 31  | :ref:`Superpos<super_pos>`  calculation (automatically inserted after search)            | ``sup``, ``super``, ``superpos``                   |
+-----+------------------------------------------------------------------------------------------+----------------------------------------------------+
| 4   | :ref:`Domain calculation<domain_calculation>`                                            | ``dom``, ``domain``, ``domains``,                  |
|     |                                                                                          | ``domain_calculation``                             |
+-----+------------------------------------------------------------------------------------------+----------------------------------------------------+
| 5   | :ref:`Error calculation<error_calculation>`                                              | ``err``, ``error``, ``errors``,                    |
|     |                                                                                          | ``error_calculation``                              |
+-----+------------------------------------------------------------------------------------------+----------------------------------------------------+
| 6   | :ref:`Full-dynamic optimization<Fdoptimization>`                                         | ``opt``, ``optimize``, ``fdopt``, ``fd_opt``,      |
|     |                                                                                          | ``full_dynamic_optimization``                      |
+-----+------------------------------------------------------------------------------------------+----------------------------------------------------+

.. versionchanged:: 0.15.0
   More aliases are accepted (e.g., ``fd_opt`` or ``structure_search``).
   A name that starts with an alias is interpreted as the section with
   the longest matching alias. Before, it was interpreted as the first
   section with a matching alias, so that, e.g., ``fd_opt`` was a
   reference calculation. Names with hyphens are ranges, except where
   the hyphen replaces the underscore before ``calc``, ``calculation``,
   ``amplitudes``, ``opt`` or ``optimization``.

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
