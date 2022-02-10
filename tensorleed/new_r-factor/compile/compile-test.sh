#!/bin/bash

cp ../wrapf2py/rfactor.o .
cp ../wrapf2py/interpolation.o .
cp ../wrapf2py/r_factor_new.mod .
cp ../wrapf2py/interpolation.mod .

gfortran -Og -o rfacsb.o -c rfacsb.f -llapack -lpthread -lblas -fallow-argument-mismatch -g -fbacktrace -fcheck=all -std=legacy -static
gfortran -Og -o main.o -c rfactor_legacy.f -llapack -lpthread -lblas -fallow-argument-mismatch -g -fbacktrace -fcheck=all -Wuninitialized -std=legacy -cpp -DNEW
gfortran -Og -o rfactor-test rfacsb.o main.o rfactor.o interpolation.o -llapack -lpthread -lblas -fallow-argument-mismatch -g -fbacktrace -fcheck=all -Wuninitialized -cpp -DNEW
#rm *.o

python3 test.py
