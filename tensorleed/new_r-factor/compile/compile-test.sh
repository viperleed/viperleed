#!/bin/bash

cp ../wrapf2py/rfactor.o .
cp ../wrapf2py/interpolation.o .
cp ../wrapf2py/r_factor_new.mod .
cp ../wrapf2py/interpolation.mod .

gfortran -Og -o rfacsb.o -c rfacsb.f -llapack -lpthread -lblas -w -fallow-argument-mismatch -g -fbacktrace -fcheck=all
gfortran -Og -o main.o -c rfactor_legacy.f -llapack -lpthread -lblas -w -fallow-argument-mismatch -g -fbacktrace -fcheck=all
gfortran -Og -o rfactor-test rfacsb.o main.o rfactor.o interpolation.o -llapack -lpthread -lblas -w -fallow-argument-mismatch -g -fbacktrace -fcheck=all -Wuninitialized
#rm *.o

python3 test.py
