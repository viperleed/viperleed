.. _searchparsinfo:

searchpars.info
===============

The searchpars.info file is generated to help with interpreting raw 
TensErLEED output of the search, i.e. the 
:ref:`control.chem<controlchem>` and :ref:`sdtl` files.
The search handles parameters under variation as a long list of indices 
that can vary over a given range, but gives those parameters rather 
unhelpful labels.
The searchpars.info file lists parameters in the same order as they 
occur in the :ref:`control.chem<controlchem>`  and :ref:`sdtl` 
files, and gives the following information on each parameter:

-   The number (at what point it occurs in the TensErLEED output).
-   The label of the parameter in the TensErLEED output.
-   Which atom is affected by the parameter 
    (same atom numbering as in :ref:`POSCAR`).
-   Which element of that atom is affected by the parameter 
    (only relevant if :ref:`ELEMENT_MIX` is used).
-   What is being varied by the parameter, i.e. **geo**\ metry, 
    **vib**\ rational amplitude, or site **occ**\ upation.
-   The number of steps in the parameter's displacements range.
-   Potential constraints, i.e. whether the parameter is linked to any 
    other parameters, or fixed to a specific index. 
    A constraint of i.e. "#3" means that the parameter is linked to 
    parameter number 3, as listed in the searchpars.info file.
