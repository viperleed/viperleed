#!/bin/bash

gfortran -c -g -Og -llapack ../interpolation/interpolation.f90
cp ../interpolation/interpolation.mod ./interpolation.mod
cp ../interpolation/interpolation.o ./interpolation.o
cp ../rfactor.f90 ./rfactor.f90
gfortran -c -g -Og rfactor.f90


f2py -m rfactor rfactor.f90 -h rfactor.pyf --overwrite-signature --debug-capi only: prepare_beams
f2py -c --verbose rfactor.pyf rfactor.f90 -I interpolation.o

python3 rf_test.py