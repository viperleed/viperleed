
gcomp := gfortran

GFOPTFLAGS := -Ofast -fno-finite-math-only

all: beamgen.v1.7 EEASiSSS.x

beamgen.v1.7: beamgen_source/beamgen.v1.7.f
	$(gcomp) beamgen_source/beamgen.v1.7.f -o beamgen.v1.7 $(GFOPTFLAGS)


EEASiSSS.x: eeasisss_code/modified/imported_routines.f90 eeasisss_code/modified/eeasisss.f90
	$(gcomp) eeasisss_code/modified/imported_routines.f90 eeasisss_code/modified/eeasisss.f90 -o EEASiSSS.x $(GFOPTFLAGS)
	rm *.mod