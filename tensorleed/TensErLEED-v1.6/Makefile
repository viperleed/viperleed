gcc_mpi_ccomp := mpicc
gcc_ccomp := gcc

intel_ccomp := icc
intel_mpi_ccomp := mpiicc

GCC_CFLAGS := -O3 

INTEL_CFLAGS := -O3

intel: lib/random_.c
	$(intel_mpi_ccomp) lib/random_.c -c -o lib/MPIrandom_.o $(INTEL_CFLAGS)
	$(intel_ccomp) lib/random_.c -c -o lib/random_.o $(INTEL_CFLAGS)

gcc: lib/random_.c
	$(gcc_mpi_ccomp) lib/random_.c -c -o lib/MPIrandom_.o $(GCC_CFLAGS)
	$(gcc_ccomp) lib/random_.c -c -o lib/random_.o $(GCC_CFLAGS)

clean:
	rm lib/random_.o lib/MPIrandom_.o