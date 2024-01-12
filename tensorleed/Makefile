gcc_fcomp := gfortran
gcc_ccomp := gcc

intel_fcomp := ifort
intel_comp := icc

GCC_FFLAGS := -O3 -fno-finite-math-only
GCC_CFLAGS := -O3 -fno-finite-math-only

INTEL_FFLAGS := -O3
INTEL_CFLAGS := -O3

# Set the sha256 command to be used for checksum verification
SHA256_CMD = sha256sum
ifeq ($(shell uname), Darwin)
        SHA256_CMD = openssl sha256 -r
endif


EEASISS_IMPORTED_ROUTINES_SHA256 := a3de67cdad38a5ed5980aa4be5f2d86a8362f2e5d1b68c12524eb26cd89926a5
EEASISS_SHA256 := 919ebb06e72dff8f914b392d66fec6ffc2413b8f0f714db4b6682a7562784e3e


intel: eeasisss_code/modified/imported_routines.f90 eeasisss_code/modified/eeasisss.f90 correct_sha256
	$(intel_fcomp) eeasisss_code/modified/imported_routines.f90 eeasisss_code/modified/eeasisss.f90 -o eeasisss $(INTEL_FFLAGS)


gcc: eeasisss_code/modified/imported_routines.f90 eeasisss_code/modified/eeasisss.f90 correct_sha256
	$(gcc_fcomp) eeasisss_code/modified/imported_routines.f90 eeasisss_code/modified/eeasisss.f90 -o eeasisss $(GCC_FFLAGS)

correct_sha256:
	@IMPORTED_ROUTINES_ACTUAL_SHA256=`cat eeasisss_code/modified/imported_routines.f90 | tr -d '\r' | $(SHA256_CMD)`; \
		case "$$IMPORTED_ROUTINES_ACTUAL_SHA256 " in \
			($(EEASISS_IMPORTED_ROUTINES_SHA256)\ *) : ok ;; \
			(*) echo eeasisss_code/modified/imported_routines.f90 checksum mismatch, expected=\"$(EEASISS_IMPORTED_ROUTINES_SHA256)\" actual=\"$$IMPORTED_ROUTINES_ACTUAL_SHA256\"; \
			exit 1 ;; \
		esac
	@EEASISS_ACTUAL_SHA256=`cat eeasisss_code/modified/eeasisss.f90 | tr -d '\r' | $(SHA256_CMD)`; \
		case "$$EEASISS_ACTUAL_SHA256 " in \
			($(EEASISS_SHA256)\ *) : ok ;; \
			(*) echo eeasisss_code/modified/eeasisss.f90 checksum mismatch, expected=\"$(EEASISS_SHA256)\" actual=\"$$EEASISS_ACTUAL_SHA256\"; \
			exit 1 ;; \
		esac

clean: 
	rm eeasisss