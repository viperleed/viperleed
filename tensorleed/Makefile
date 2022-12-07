gcc_fcomp := gfortran
gcc_ccomp := gcc

intel_fcomp := ifort
intel_comp := icc

GCC_FFLAGS := -O3 -fno-finite-math-only
GCC_CFLAGS := -O3 -fno-finite-math-only

INTEL_FFLAGS := -O3
INTEL_CFLAGS := -O3


intel: beamgen_source/beamgen.v1.7.f eeasisss_code/modified/imported_routines.f90 eeasisss_code/modified/eeasisss.f90
	$(intel_fcomp) beamgen_source/beamgen.v1.7.f -o beamgen.v1.7 $(INTEL_FFLAGS)
	$(intel_fcomp) eeasisss_code/modified/imported_routines.f90 eeasisss_code/modified/eeasisss.f90 -o EEASiSSS.x $(INTEL_FFLAGS)


gcc: beamgen_source/beamgen.v1.7.f eeasisss_code/modified/imported_routines.f90 eeasisss_code/modified/eeasisss.f90
	$(gcc_fcomp) beamgen_source/beamgen.v1.7.f -o beamgen.v1.7 $(GCC_FFLAGS)
	$(gcc_fcomp) eeasisss_code/modified/imported_routines.f90 eeasisss_code/modified/eeasisss.f90 -o EEASiSSS.x $(GCC_FFLAGS)
        
clean: 
	rm beamgen.v1.7 EEASiSSS.x