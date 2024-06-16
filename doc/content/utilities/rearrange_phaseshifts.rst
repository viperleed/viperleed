.. _rearrange_phaseshifts:

=====================
rearrange_phaseshifts
=====================

``modify_phaseshifts`` is a simple :term:`Python` utility script packaged with
ViPErLEED. It reads a :ref:`PHASESHIFTS` file and allows the user
to copy and rearrange the blocks. This can be helpful if it is desirable to
reuse a custom :ref:`PHASESHIFTS` file for structures with the
same chemical elements but different stoichiometries or different
:ref:`site definitions<sitedef>`.

Usage
=====

To use the rearrange_phaseshifts utility, call it from the command line as

.. code-block:: bash

    $ viperleed util rearrange_phaseshifts

The utility expects to be called in the directory containing the
:ref:`PHASESHIFTS` file. Once called, the script will
instruct and prompt the user on how to proceed to reorder the file.
