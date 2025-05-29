.. _experiment_symmetry:

experiment_symmetry.ini
=======================

The :file:`experiment_symmetry.ini` file is generated during the
:ref:`initialization`, and contains basic information on the unit
cell and symmetry of the structure, on the energy range used, and
on the incidence direction of the electron beam. It serves as an
input file to the GUI, which can display a simulated LEED pattern
for the structure, and export a "spot-pattern file" for the Spot
Tracker of the :ref:`imagej_plugins`.


Changes
-------

.. versionchanged:: 0.14.0
   The name and extension of the file changed from
   :file:`PatternInfo.tlm` to :file:`experiment_symmetry.ini`.
