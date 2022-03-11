
comp := gfortran

GFOPTFLAGS := -g -Og
GFDEBUGFLAGS := -fbacktrace -fcheck=all -Wuninitialized -static -fPIC
GFLAPACK := -llapack

# Target are the .so python module executables
all: rfactor.so interpolation.so
	cp *.so ../../tleedmlib/wrapped
	echo "Success: Modules copied to tleedmlib/wrapped"

rfactor.so: interpolation.o interpolation.mod rfactor.o rfactor.mod rfactor.pyf rfactor/rfactor.f90 interpolation/interpolation.f90
	f2py -c --f90flags="$(GFDEBUGFLAGS)" rfactor.pyf rfactor/rfactor.f90 -I interpolation/interpolation.f90 $(GFLAPACK)
	mv rfactor.cpython**.so rfactor.so

	
interpolation.so: interpolation.pyf interpolation.o interpolation.mod
	f2py -c --f90flags="$(GFDEBUGFLAGS)" interpolation.pyf interpolation/interpolation.f90 --opt="$(GFOPTFLAGS)" --debug-capi $(GFLAPACK)
	mv interpolation.cpython**.so interpolation.so

rfactor.pyf: rfactor/rfactor.f90
	f2py -m rfactor rfactor/rfactor.f90 -h rfactor.pyf --overwrite-signature --debug-capi only: prepare_beams r_pendry_beam_y r_pendry_beamset_y r_beamset_v0r_opt_on_grid parabola_lsq_fit parabola parabola_r_squared r_beamtype_grouping r2_beam_intensity r2_beamset_intensity apply_beamset_shift trapez_integration_const_dx tenser_intsum r_pendry_beam_y_wrong

interpolation.pyf: interpolation/interpolation.f90
	f2py -m interpolation interpolation/interpolation.f90 -h interpolation.pyf --overwrite-signature --debug-capi only: get_array_sizes de_Boor_size single_calc_spline pre_eval_input_grid calc_spline_with_pre_eval single_interpolate_coeffs_to_grid get_intervals calc_deBoor eval_bspline_fast

rfactor.o rfactor.mod: interpolation.o interpolation.mod
	$(comp) -c rfactor/rfactor.f90 $(GFLAPACK) $(GFOPTFLAGS) $(GFDEBUGFLAGS)

interpolation.o interpolation.mod: interpolation/interpolation.f90
	$(comp) -c interpolation/interpolation.f90 $(GFLAPACK) $(GFOPTFLAGS) $(GFDEBUGFLAGS)


.PHONY: clean
clean:
	rm *.o *.mod *.pyf *.so