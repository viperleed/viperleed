comp := gfortran


BUILD_DIR := ./build
SRC_DIR := .

# Set the sha256 command to be used for checksum verification
SHA256_CMD = sha256sum
ifeq ($(shell uname), Darwin)
        SHA256_CMD = openssl sha256 -r
endif

SMOOTHING_SHA256 := 4df3d0ee921802b85de42108e3d749fa00a1a0263119f15f9e10c168a1ee2996
INTERPOLATION_SHA256 := 403ebddec583fd88950ce923b91811e2b7ba9193c8a76e2d27b3ba6063d1400f
RFACTOR_SHA256 := c0d9fbae6ece536e17f1dae8d9cef35c4349744469ae6460a6cee7429fad4b9b

# Make sure that -fno-finite-math-only is enabled when using aggressive optimization, otherwise checks WILL fail! 
GFOPTFLAGS := -g -Ofast -fno-finite-math-only -fsignaling-nans
GFDEBUGFLAGS := -fbacktrace -fcheck=all -Wuninitialized -static -fPIC
GFLAPACK := -llapack

# Target are the .so python module executables
all: rfactor.so interpolation.so MS_smoothing.so
	echo "Success: Modules saved to tleedmlib/wrapped"

rfactor.so: $(BUILD_DIR)/interpolation.o $(BUILD_DIR)/interpolation.mod $(BUILD_DIR)/rfactor.o $(BUILD_DIR)/rfactor.mod $(BUILD_DIR)/rfactor.pyf $(SRC_DIR)/rfactor/rfactor.f90 $(SRC_DIR)/interpolation/interpolation.f90 correct_sha256
	f2py -c --opt="$(GFOPTFLAGS)" --f90flags="$(GFDEBUGFLAGS) -I$(BUILD_DIR)" $(BUILD_DIR)/rfactor.pyf $(SRC_DIR)/rfactor/rfactor.f90 -I $(SRC_DIR)/interpolation/interpolation.f90 $(GFLAPACK)
	mv rfactor.cpython**.so rfactor.so

	
interpolation.so: $(BUILD_DIR)/interpolation.pyf $(BUILD_DIR)/interpolation.o $(BUILD_DIR)/interpolation.mod correct_sha256
	f2py -c --opt="$(GFOPTFLAGS)" --f90flags="$(GFDEBUGFLAGS) -I$(BUILD_DIR)" $(BUILD_DIR)/interpolation.pyf $(SRC_DIR)/interpolation/interpolation.f90 --debug-capi $(GFLAPACK)
	mv interpolation.cpython**.so interpolation.so

$(BUILD_DIR)/rfactor.pyf: $(SRC_DIR)/rfactor/rfactor.f90 correct_sha256
	mkdir -p $(BUILD_DIR)/
	f2py -m rfactor $(SRC_DIR)/rfactor/rfactor.f90 -h $(BUILD_DIR)/rfactor.pyf --overwrite-signature --debug-capi only: prepare_beams r_pendry_beam_y r_pendry_beamset_y r_beamset_v0r_opt_on_grid parabola_lsq_fit parabola parabola_r_squared r_beamtype_grouping r2_beam_intensity r2_beamset_intensity apply_beamset_shift trapez_integration_const_dx tenser_intsum r_pendry_beam_y_wrong test_compilation pendry_y_beamset pendry_y alloc_beams_arrays
$(BUILD_DIR)/interpolation.pyf: $(SRC_DIR)/interpolation/interpolation.f90 correct_sha256
	mkdir -p $(BUILD_DIR)/
	f2py -m interpolation $(SRC_DIR)/interpolation/interpolation.f90 -h $(BUILD_DIR)/interpolation.pyf --overwrite-signature --debug-capi only: get_array_sizes de_Boor_size single_calc_spline pre_eval_input_grid calc_spline_with_pre_eval single_interpolate_coeffs_to_grid get_intervals calc_deBoor eval_bspline_fast

$(BUILD_DIR)/rfactor.o $(BUILD_DIR)/rfactor.mod: $(BUILD_DIR)/interpolation.o $(BUILD_DIR)/interpolation.mod correct_sha256
	mkdir -p $(BUILD_DIR)/
	$(comp) -c $(SRC_DIR)/rfactor/rfactor.f90 $(GFLAPACK) $(GFOPTFLAGS) $(GFDEBUGFLAGS) -I$(BUILD_DIR) -J$(BUILD_DIR) -o $(BUILD_DIR)/interpolation.o

$(BUILD_DIR)/interpolation.o $(BUILD_DIR)/interpolation.mod: $(SRC_DIR)/interpolation/interpolation.f90 correct_sha256
	mkdir -p $(BUILD_DIR)/
	$(comp) -c $(SRC_DIR)/interpolation/interpolation.f90 $(GFLAPACK) $(GFOPTFLAGS) $(GFDEBUGFLAGS) -I$(BUILD_DIR) -J$(BUILD_DIR) -o $(BUILD_DIR)/rfactor.o
	
MS_smoothing.so: $(BUILD_DIR)/MS_smoothing.pyf correct_sha256
	f2py -c $(BUILD_DIR)/MS_smoothing.pyf $(SRC_DIR)/smoothing/MS_smoothing.f90 --f90flags="$(GFOPTFLAGS)" --opt="$(GFOPTFLAGS)"
	mv MS_smoothing.cpython**.so MS_smoothing.so
	
$(BUILD_DIR)/MS_smoothing.pyf: correct_sha256
	mkdir -p $(BUILD_DIR)/
	f2py -m MS_smoothing $(SRC_DIR)/smoothing/MS_smoothing.f90 -h $(BUILD_DIR)/MS_smoothing.pyf --overwrite-signature --debug-capi

correct_sha256:
	SMOOTHING_ACTUAL_SHA=`$(SHA256_CMD) $(SRC_DIR)/smoothing/MS_smoothing.f90`; \
		case "$$SMOOTHING_ACTUAL_SHA " in \
			($(SMOOTHING_SHA256)\ *) : ok ;; \
			(*) echo $(SRC_DIR)/smoothing/MS_smoothing.f90 checksum mismatch, expected=\"$(SMOOTHING_SHA256)\" actual=\"$$SMOOTHING_ACTUAL_SHA\"; \
			exit 1 ;; \
		esac
	INTERPOLATION_ACTUAL_SHA=`$(SHA256_CMD) $(SRC_DIR)/interpolation/interpolation.f90`; \
		case "$$INTERPOLATION_ACTUAL_SHA " in \
			($(INTERPOLATION_SHA256)\ *) : ok ;; \
			(*) echo $(SRC_DIR)/interpolation/interpolation.f90 checksum mismatch, expected=\"$(INTERPOLATION_SHA256)\" actual=\"$$INTERPOLATION_ACTUAL_SHA\"; \
			exit 1 ;; \
		esac
	RFACTOR_ACTUAL_SHA=`$(SHA256_CMD) $(SRC_DIR)/rfactor/rfactor.f90`; \
		case "$$RFACTOR_ACTUAL_SHA " in \
			($(RFACTOR_SHA256)\ *) : ok ;; \
			(*) echo $(SRC_DIR)/rfactor/rfactor.f90 checksum mismatch, expected=\"$(RFACTOR_SHA256)\" actual=\"$$RFACTOR_ACTUAL_SHA\"; \
			exit 1 ;; \
		esac


.PHONY: clean
clean:
	rm -r $(BUILD_DIR)/