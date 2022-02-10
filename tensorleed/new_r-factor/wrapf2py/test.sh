#!/bin/bash

gfortran -c -g -Og -llapack ../interpolation/interpolation.f90 -fbacktrace -fcheck=all -Wuninitialized -static
#cp ../interpolation/interpolation.mod ./interpolation.mod
#cp ../interpolation/interpolation.o ./interpolation.o
#f2py -m interpolation ../interpolation/interpolation.f90 -h interpolation.pyf --overwrite-signature --debug-capi only: get_array_sizes de_Boor_size single_calc_spline pre_eval_input_grid calc_spline_with_pre_eval single_interpolate_coeffs_to_grid get_intervals calc_deBoor eval_bspline_fast
#f2py -c --f90flags="-g -Og -llapack -fbacktrace -fcheck=all" interpolation.pyf ../interpolation/interpolation.f90

cp ../rfactor.f90 ./rfactor.f90
gfortran -c -g -Og rfactor.f90 -fbacktrace -fcheck=all -Wuninitialized -static

# Subroutines must be given in lowercase letters here
#f2py -m rfactor rfactor.f90 -h rfactor.pyf --overwrite-signature --debug-capi only: prepare_beams r_pendry_beam_y r_pendry_beamset_y r_pendry_beamset_v0r_opt_on_grid parabola_lsq_fit parabola parabola_r_squared
#f2py -c --f90flags="-g -Og -fbacktrace -fcheck=all" rfactor.pyf rfactor.f90 -I interpolation.o

#ipython -i rf_test.py