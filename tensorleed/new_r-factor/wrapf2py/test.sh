#!/bin/bash

gfortran -c -g -Og -llapack ../interpolation/interpolation.f90
cp ../interpolation/interpolation.mod ./interpolation.mod
cp ../interpolation/interpolation.o ./interpolation.o
cp ../rfactor.f90 ./rfactor.f90
gfortran -c -g -Og rfactor.f90


f2py -m rfactor rfactor.f90 -h rfactor.pyf --overwrite-signature --debug-capi only: prepare_beams r_pendry_beam_y r_pendry_beamset_y
f2py -c --f90flags="-g -Og" rfactor.pyf rfactor.f90 -I interpolation.o

ipython -i rf_test.py