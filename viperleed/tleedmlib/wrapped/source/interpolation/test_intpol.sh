#!/bin/bash
#gfortran -Og fortran-utils-master/src/lapack.f90 fortran-utils-master/src/types.f90 fortran-utils-master/src/utils.f90 fortran-utils-master/src/splines.f90 fortran-utils-master/src/intpol_test_splines.f90 -fcheck=all -fbacktrace -Wall
rm *.mod
#gfortran -cpp -c bspline_kinds_module.F90 bspline_sub_module.f90 bspline_oo_module.f90 bspline_module.f90
gfortran -cpp TensErLEED_intpol.f90   lapack.f90 types.f90 utils.f90 splines.f90 interpolation.f90 intpol_test.f90 -llapack -o intpol_test.out -fcheck=all -fbacktrace -Ofast
./intpol_test.out

rm *.mod

echo "testscript finished"
