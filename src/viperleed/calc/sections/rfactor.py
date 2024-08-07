"""Section Rfactor."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-11'
__license__ = 'GPLv3+'

import logging
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

from viperleed.calc import DEFAULT_WORK
from viperleed.calc.files import iorfactor
from viperleed.calc.files import iotensors
from viperleed.calc.files.iorefcalc import readFdOut
from viperleed.calc.lib import leedbase
from viperleed.calc.lib.checksums import validate_multiple_files

try:
    from viperleed.calc.wrapped.rfactor import r_factor_new as rf
except ImportError:
    _HAS_NEW_R_FACTOR = False
else:
    from viperleed.calc.wrapped.error_codes import check_ierr
    _HAS_NEW_R_FACTOR = True

logger = logging.getLogger(__name__)


def rfactor(sl, rp, index, for_error=False, only_vary=None):                    # TODO: Parameters __doc__
    """Runs the r-factor calculation for either the reference calculation
    (index 11) or the superpos (index 12)."""
    if rp.THEO_ENERGIES.n_energies < 2:
        logger.info("Only one theoretical energy found: Cannot "
                    "calculate a meaningful R-Factor. Stopping...")
        return []
    if for_error and rp.best_v0r is None:
        # Cannot proceed with planned section
        err_msg = ("Cannot calculate R-factor for error without a stored"
                   "V0r value. Execute a normal R-factor calculation first.")
        rp.setHaltingLevel(3)
        logger.error(err_msg)
        raise RuntimeError(err_msg)

    # Both refcalc (11) and superpos should work fine with new R-factor
    if index == 11:
        name = "refcalc"
    else:
        name = "superpos"
    _need_spectra = (index == 11 and not len(rp.refcalc_fdout)
                     or index == 12 and not len(rp.superpos_specout))
    if _need_spectra:
        _fetch_and_check_spectra(rp, index, name)
    theobeams = rp.theobeams[name]
    expbeams = rp.expbeams

    # Branch off for new R factor calculation
    if not rp.R_FACTOR_LEGACY and _HAS_NEW_R_FACTOR:                            # TODO: we may want to raise some errors if we want to use it but it's not there
        return run_new_rfactor(sl, rp, for_error, name, theobeams, expbeams)
    return run_legacy_rfactor(sl, rp, for_error, name, theobeams, index, only_vary)


def _fetch_and_check_spectra(rp, index, name):
    """Look for full-dynamic spectra, then check consistency."""
    if index == 11:
        fn = Path("refcalc-fd.out")
    else:
        fn = Path("superpos-spec.out")
    directory = None
    path = None
    if fn.is_file():
        directory = DEFAULT_WORK
        path = fn
    elif ("OUT" / fn).is_file():
        directory = "OUT"
        path = "OUT" / fn
    elif index == 11:
        # try getting from Tensors
        iotensors.getTensors(rp.TENSOR_INDEX, required=False)
        directory = Path(f"Tensors_{rp.TENSOR_INDEX:03d}")
        if ("Tensors" / directory / fn).is_file():
            path = "Tensors" / directory / fn

    if path:
        logger.warning("R-factor calculation was called without stored "
                       f"spectrum data. Reading from file {fn} in {directory} "
                       "folder...")
    else:
        logger.error("Cannot execute R-factor calculation: no stored "
                     f"spectrum data and no {fn} file was found.")
        raise RuntimeError("No spectrum data found")
    try:
        theobeams, theospec = readFdOut(readfile=path)
    except Exception:
        logger.error(f"Failed to read {path}")
        raise

    # if we haven't returned, then it was read. check data vs ivbeams:
    if index == 11:
        rp.refcalc_fdout = theospec
    else:
        rp.superpos_specout = theospec
    rp.theobeams[name] = theobeams

    eps = 1e-3
    consistent = (len(rp.ivbeams) == len(theobeams)
                  and all(any(ivbeam.isEqual(theobeam, eps=eps)
                              for theobeam in theobeams)
                              for ivbeam in rp.ivbeams))
    if not consistent:
        logger.error("The list of beams read from IVBEAMS is not equivalent "
                     f"to the list of beams in {path}. R-Factor calculation "
                     "cannot proceed.")
        raise ValueError("Contradiction in beam sets")


def run_new_rfactor(sl, rp, for_error, name, theobeams, expbeams):
    logger.debug("Using new R-factor calculation. This is still experimental!")
    which_r = rp.R_FACTOR_TYPE

    if which_r == 1:
        n_derivs = 1
    elif which_r == 2:
        n_derivs = 0
    else:
        check_ierr(701, logger)

    (_, theo_range,
     iv_shift,
     intpol_step) = iorfactor.prepare_rfactor_energy_ranges(
        rp, theobeams, for_error,
        n_expand=(rp.INTPOL_DEG - 1) // 2
        )
    out_grid = np.arange(theo_range.min,
                         theo_range.max + 0.1 * intpol_step,
                         intpol_step)

    # find correspondence experimental to theoretical beams:
    beamcorr = leedbase.getBeamCorrespondence(sl, rp)
    # integer & fractional beams
    iorf = []
    for (i, beam) in enumerate(rp.expbeams):
        if beam.hk[0] % 1.0 != 0.0 or beam.hk[1] % 1.0 != 0.0:
            iorf.append(2)
        else:
            iorf.append(1)
    iorf.extend([0] * (len(rp.ivbeams) - len(rp.expbeams)))

    n_beams = len(beamcorr) - beamcorr.count(-1)  # -1 indicates no correspondence
    (exp_grid, exp_id_start,
     exp_n_E_beams,
     exp_intensities_in) = iorfactor.beamlist_to_array(rp.expbeams)
    corr_theobeams = []
    for i in beamcorr:
        if i != -1:
            corr_theobeams.append(theobeams[i])
    (theo_grid,
     theo_id_start,
     theo_n_E_beams,
     theo_intensities_in) = iorfactor.beamlist_to_array(corr_theobeams)

    deg = rp.INTPOL_DEG
    v0i = rp.V0_IMAG

    # Should we do reordering averaging etc. using the subroutine? - unused right now... # TODO
    averaging_scheme = np.int32(np.arange(n_beams) + 1)  # Fortran index

    # create buffer arrays
    (
        exp_e_start_beams_out,
        exp_n_e_beams_out,
        exp_intpol_int,
        exp_intpol_deriv_y,
        ierr,
    ) = rf.alloc_beams_arrays(n_beams, len(out_grid))
    check_ierr(ierr, logger)

    # interpolation and derivative
    ierr = rf.prepare_beams(
        e_grid_in=exp_grid,
        intensities=exp_intensities_in,
        e_start_beams=exp_id_start + 1,  # Fortran indices
        n_e_beams=exp_n_E_beams,
        deg=deg,
        n_derivs=n_derivs,
        e_grid_out=out_grid,
        e_start_beams_out=exp_e_start_beams_out,
        n_e_beams_out=exp_n_e_beams_out,
        intpol_intensity=exp_intpol_int,
        intpol_derivative=exp_intpol_deriv_y,
    )
    # check for erros
    check_ierr(ierr, logger)

    # Same for theoretical beams now
    # create buffer arrays
    (
        theo_e_start_beams_out,
        theo_n_e_beams_out,
        theo_intpol_int,
        theo_intpol_deriv_y,
        ierr,
    ) = rf.alloc_beams_arrays(n_beams, len(out_grid))
    check_ierr(ierr, logger)

    # interpolation and derivative
    ierr = rf.prepare_beams(
        e_grid_in=theo_grid,
        intensities=theo_intensities_in,
        e_start_beams=theo_id_start + 1,  # Fortran indices
        n_e_beams=theo_n_E_beams,
        deg=deg,
        n_derivs=n_derivs,
        e_grid_out=out_grid,
        e_start_beams_out=theo_e_start_beams_out,
        n_e_beams_out=theo_n_e_beams_out,
        intpol_intensity=theo_intpol_int,
        intpol_derivative=theo_intpol_deriv_y,
    )
    # check for erros
    check_ierr(ierr, logger)

    # For Factor R2 we use intesities for further calculations, otherwise the Y function
    if which_r == 1:
        # Make Y-functions
        rf.pendry_y_beamset(
            exp_intpol_int,
            exp_intpol_deriv_y,  # derivative is overwritten by the Y function
            exp_e_start_beams_out,
            exp_n_e_beams_out,
            v0i,
        )

        rf.pendry_y_beamset(
            theo_intpol_int,
            theo_intpol_deriv_y,  # derivative is overwritten by the Y function
            theo_e_start_beams_out,
            theo_n_e_beams_out,
            v0i,
        )

        exp_data = exp_intpol_deriv_y
        theo_data = theo_intpol_deriv_y

    elif which_r == 2:
        exp_data = exp_intpol_int
        theo_data = theo_intpol_int

        # optimize V0r and calculate R factor

        ## settings for V0r optimization
    bounds = iv_shift.min, iv_shift.max
    v0r_range = np.array(
        [round(bound / intpol_step) for bound in bounds], dtype="int32"
    )
    v0r_center = int((v0r_range[0] + v0r_range[1]) / 2)
    start_guess = np.array([v0r_range[0], v0r_center, v0r_range[1]], dtype="int32")

    # TODO: make below into User Parameters
    fast_search = False
    tol_r = (1 - 5e-2,)
    tol_r_2 = (1 - 5e-2,)
    max_fit_range = 6

    (
        best_v0r_step,
        best_v0r,
        best_rfac,
        _,  # n_v0r_evaluated
        R_beams,
        numerators,
        denominators,
        n_overlapp_beams,
        ierr,
    ) = rf.r_beamset_v0r_opt_on_grid(
        which_r,
        v0r_range,
        start_guess,
        fast_search,
        intpol_step,
        exp_data,
        theo_data,
        exp_e_start_beams_out,
        theo_e_start_beams_out,
        exp_n_e_beams_out,
        theo_n_e_beams_out,
        tol_r=tol_r,
        tol_r_2=tol_r_2,
        max_fit_range=max_fit_range,
    )

    check_ierr(ierr, logger)

    logger.info(
        "With inner potential shift of {:.2f} eV: "
        "R = {:.4f}\n".format(best_v0r, best_rfac)
    )
    rp.best_v0r = best_v0r

    # seperate into integer and fractional beams
    n_groups = 2
    grouping = np.int32(iorf)  # 1 is integer, 2 is fractional

    # Call grouping routine
    (r_factor_groups, _, ierr) = rf.r_beamtype_grouping(  # n_overlapp_groups unused
        which_r, R_beams, numerators, denominators, n_overlapp_beams, n_groups, grouping
    )
    # error 903 and 904 are acceptable - that just means there are no integer/fractional beams

    check_ierr(ierr, logger)

    (r_int, r_frac) = list(r_factor_groups)

    if np.isclose(best_rfac, 0):
        logger.error(
            "R-Factor reported as zero. This means "
            "something went wrong in the reference "
            "calculation or in the R-factor calculation."
        )
        rp.setHaltingLevel(2)
    elif best_rfac > 2:
        logger.error(
            f"R-Factor reported as {best_rfac}. This means "
            "something went wrong in the reference "
            "calculation or in the R-factor calculation."
        )
        rp.setHaltingLevel(2)
    else:
        rp.last_R = best_rfac
        rp.stored_R[name] = (best_rfac, r_int, r_frac)

    if rp.PLOT_IV["plot"]:
        # Plot R-factor
        outname = "Rfactor_plots_{}.pdf".format(name)
        aname = "Rfactor_analysis_{}.pdf".format(name)
        if rp.PLOT_IV["overbar"]:
            labelstyle = "overbar"
        else:
            labelstyle = "minus"
        # labelwidth = max([beam.getLabel(style=labelstyle)[1] for beam in rp.expbeams]) # TODO currently unused? why?

        labels = [beam.label for beam in expbeams]

        # Apply v0r shift to arrays for plotting

        shifted_intpol_exp = np.copy(exp_intpol_int)
        shifted_intpol_theo = np.copy(theo_intpol_int)
        shifted_y_exp = np.copy(exp_data)
        shifted_y_theo = np.copy(theo_data)
        shifted_E_start_exp = np.copy(exp_e_start_beams_out)
        shifted_E_start_theo = np.copy(theo_e_start_beams_out)
        shifted_n_exp = np.copy(exp_n_e_beams_out)
        shifted_n_theo = np.copy(theo_n_e_beams_out)

        """
            np.savetxt("viper_exp_start.csv", exp_e_start_beams_out, delimiter = ',')
            np.savetxt("viper_exp_n_beams.csv", exp_n_e_beams_out, delimiter = ',')
            np.savetxt("viper_exp_y.csv", exp_yfunc, delimiter = ',')
            np.savetxt("viper_exp_intensity.csv", exp_intpol_intensity, delimiter = ',')
            np.savetxt("viper_theo_start.csv", theo_e_start_beams_out, delimiter = ',')
            np.savetxt("viper_theo_n_beams.csv", theo_n_e_beams_out, delimiter = ',')
            np.savetxt("viper_theo_y.csv", theo_yfunc, delimiter = ',')
            np.savetxt("viper_theo_intensity.csv", theo_intpol_intensity, delimiter = ',')
            np.savetxt("viper_theo_energies.csv", theo_grid, delimiter = ',')
            np.savetxt("viper_exp_energies.csv", exp_grid, delimiter = ',')
            """

        ierr = rf.apply_beamset_shift(
            shifted_intpol_exp,
            shifted_E_start_exp,
            shifted_n_exp,
            shifted_intpol_theo,
            shifted_E_start_theo,
            shifted_n_theo,
            shift=best_v0r_step,
            fill_outside=1,
        )
        check_ierr(ierr, logger)
        ierr = rf.apply_beamset_shift(
            shifted_y_exp,
            np.copy(exp_e_start_beams_out),
            np.copy(exp_n_e_beams_out),
            shifted_y_theo,
            np.copy(theo_e_start_beams_out),
            np.copy(theo_n_e_beams_out),
            shift=best_v0r_step,
            fill_outside=1,
        )
        check_ierr(ierr, logger)

        # shifted_E_start_exp must be now same as shifted_E_start_theo -> no need to pass both
        # same for n_e_beams
        iorfactor.writeRfactorPdf_new(
            n_beams,
            labels,
            R_beams,
            out_grid,
            shifted_E_start_exp,
            shifted_n_exp,
            shifted_intpol_exp,
            shifted_intpol_theo,
            shifted_y_exp,
            shifted_y_theo,
            outName=outname,
            analysisFile=aname,
            v0i=rp.V0_IMAG,
            formatting=rp.PLOT_IV,
        )
    rfaclist = list(R_beams)
    return rfaclist


def run_legacy_rfactor(sl, rp, for_error, name, theobeams, index, only_vary):

    if index == 11:
        theospec = rp.refcalc_fdout
    elif index == 12:
        theospec = rp.superpos_specout
    # WEXPEL before PARAM, to make sure number of exp. beams is correct
    try:
        iorfactor.writeWEXPEL(sl, rp, theobeams, for_error=for_error)
    except Exception:
        logger.error("Exception during writeWEXPEL: ")
        raise
    try:
        iorfactor.writeRfactPARAM(rp, theobeams,
                                  for_error=for_error,
                                  only_vary=only_vary)
    except Exception:
        logger.error("Exception during writeRfactPARAM: ")
        raise

    # get fortran files and compile
    try:
        tl_source = rp.get_tenserleed_directory()
        tl_path = tl_source.path
        libpath = tl_path / "lib"
        libname = next(libpath.glob("rfacsb*"))                                 # StopIteration??
        srcpath = tl_path / "src"
        srcname = next(srcpath.glob("rfactor.*"))                               # StopIteration??
        shutil.copy2(libname, libname.name)  # Copy here from source
        shutil.copy2(srcname, srcname.name)  # Copy here from source
    except Exception:
        logger.error(
            "Error getting TensErLEED files for r-factor calculation:"
            )
        raise
    if rp.SUPPRESS_EXECUTION:
        logger.warning("SUPPRESS_EXECUTION parameter is on. R-factor "
                       "calculation will not proceed. Stopping...")
        rp.setHaltingLevel(3)
        return []
    # validate checksums
    if not rp.TL_IGNORE_CHECKSUM:
        files_to_check = (Path(libpath) / libname,
                          Path(srcpath) / srcname)
        validate_multiple_files(files_to_check, logger,
                                "R-factor", rp.TL_VERSION)

    logger.info("Compiling fortran input files...")
    rfacname = f"rfactor-{rp.timestamp}"
    if rp.FORTRAN_COMP[0] == "":
        rp.getFortranComp()
    # compile
    _pre, _post = rp.FORTRAN_COMP
    ctasks = [
        (f"{_pre} -o rfacsb.o -c", libname.name, _post),
        (f"{_pre} -o main.o -c", srcname.name, _post),
        (f"{_pre} -o {rfacname}", "rfacsb.o main.o", _post)
        ]

    # Compile
    compile_log = Path("compile-rfactor.log")
    try:
        leedbase.fortran_compile_batch(ctasks, logname=compile_log)
    except Exception:
        logger.error("Error compiling fortran files: ", exc_info=True)
        # move compile log to compile_logs directory
        leedbase.copy_compile_log(rp, compile_log, "rfactor-compile")
        raise

    # move log file to supp
    try:
        shutil.move(compile_log, "compile_logs" / compile_log)
    except OSError:
        logger.warning(f"Could not move {compile_log} to SUPP")
    # run
    rfaclogname = Path(rfacname).with_suffix(".log")
    logger.info(
        "Starting R-factor calculation...\n"
        f"R-factor log will be written to file {rfaclogname}"
        )
    try:
        with rfaclogname.open("w") as log:
            subprocess.run(
                Path(rfacname).resolve(),
                input=theospec,
                encoding="ascii",
                stdout=log,
                stderr=log,
                )
    except Exception:
        logger.error("Error during R-factor calculation. "
                     "Also check R-factor log file.")
        raise
    logger.info("Finished R-factor calculation. Processing files...")
    # rename and move files
    try:
        os.rename("WEXPEL", "rfactor-WEXPEL")
    except OSError:
        logger.warning("Failed to rename R-factor input "
                       "file WEXPEL to rfactor-WEXPEL")
    try:
        os.rename("PARAM", "rfactor-PARAM")
    except OSError:
        logger.warning("Failed to rename R-factor input "
                       "file PARAM to rfactor-PARAM")
    if not Path("ROUT").is_file():
        logger.error("No ROUT file was found after R-Factor calculation!")
        rp.setHaltingLevel(2)
        return []

    # read output
    if for_error:
        try:
            rfaclist = iorfactor.readROUTSHORT()
        except Exception:
            logger.error("Error reading ROUTSHORT file",
                         exc_info=rp.is_debug_mode)
            rp.setHaltingLevel(2)
        return rfaclist

    try:
        (rfac, r_int, r_frac), v0rshift, rfaclist = iorfactor.readROUT()
    except Exception:  # TODO catch correct exception
        logger.error("Error reading ROUT file", exc_info=rp.is_debug_mode)
        rp.setHaltingLevel(2)
        return []

    logger.info("With inner potential shift of {:.2f} eV: "
                "R = {:.4f}\n".format(v0rshift, rfac))
    rp.best_v0r = v0rshift
    dir_list = [Path(), Path("OUT")]
    for dir_name in dir_list:
        for f_name in dir_name.glob(f"R_OUT_{rp.timestamp}*"):
            if not f_name.is_file():
                continue
            try:  # delete old R_OUT files
                f_name.unlink()
            except Exception:
                pass
    if rfac <= 0.00001:
        logger.warning(
            "ROUT reports R-Factor as zero. This may mean "
            "something went wrong in the reference "
            "calculation or in the R-factor calculation. "
            "If you are comparing with pseudo-experiment data, "
            "you can ignore this warning."
            )

    f_name = f"R_OUT_{name}_{rp.timestamp}_R={rfac:.4f}"
    rp.last_R = rfac
    rp.stored_R[name] = (rfac, r_int, r_frac)
    try:
        os.rename("ROUT", f_name)
    except OSError:
        logger.warning("Failed to rename R-factor "
                       f"output file ROUT to {f_name}")
    if len(rfaclist) != len(rp.expbeams):
        logger.warning("Failed to read R-Factors per "
                       "beam from R-factor output file ROUT.")
        rfaclist = [-1] * len(rp.expbeams)

    if rp.PLOT_IV["plot"]:
        outname = f"Rfactor_plots_{name}.pdf"
        aname = f"Rfactor_analysis_{name}.pdf"
        labelstyle = "overbar" if rp.PLOT_IV["overbar"] else "minus"
        labelwidth = max(beam.getLabel(style=labelstyle)[1]
                         for beam in rp.expbeams)
        try:
            iorfactor.writeRfactorPdf(
                [(b.getLabel(lwidth=labelwidth, style=labelstyle)[0], r)
                 for b, r in zip(rp.expbeams, rfaclist)],
                outName=outname,
                analysisFile=aname,
                v0i=rp.V0_IMAG,
                formatting=rp.PLOT_IV,
                )
        except Exception:                                                       # TODO catch correct exception
            logger.warning("Error plotting R-factors.", exc_info=True)
    return rfaclist
