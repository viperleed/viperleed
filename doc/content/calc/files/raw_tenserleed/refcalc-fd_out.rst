.. _fd_out:

refcalc-fd.out and refcalc-amp.out
==================================

TensErLEED writes the output of the :ref:`ref-calc` for each energy step, i.e.,
the diffraction intensity :math:`I` and complex amplitude :math:`A` of each
requested beam (index :math:`(h|k)`) at energy :math:`E` to two files called
``fd.out`` and ``amp.out``.

In a first step, ViPErLEED takes these outputs and combines them in the files
``refcalc-fd.out`` and ``refcalc-amp.out``. These files thus contain the raw,
unformatted output of the reference calculation. They are provided primarily
for compatibility and debugging purposes.

The same data is provided in a more readable and useful format in the
:ref:`THEOBEAMS.csv<THEOBEAMS>` file.
