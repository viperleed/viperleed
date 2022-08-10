"""
Created on Aug 11 2020

@author: Florian Kraushofer, Alexander M. Imre

Tensor LEED Manager section Rfactor
"""

import os
import logging
import shutil
import subprocess
from tabnanny import check
import numpy as np

from viperleed.tleedmlib.files.iorefcalc import readFdOut
from viperleed.tleedmlib.leedbase import fortran_compile_batch, getTLEEDdir, getTensors
import viperleed.tleedmlib.files.iorfactor as io

from viperleed.tleedmlib.wrapped.rfactor import r_factor_new as rf
from viperleed.tleedmlib.wrapped.error_codes import error_codes, check_ierr


logger = logging.getLogger("tleedm.rfactor")


def rfactor(sl, rp, index, for_error=False, only_vary=None):
    """Runs the r-factor calculation for either the reference calculation
    (index 11) or the superpos (index 12)."""
    if int((rp.THEO_ENERGIES[1] - rp.THEO_ENERGIES[0]) / rp.THEO_ENERGIES[2]) + 1 < 2:
        logger.info(
            "Only one theoretical energy found: Cannot calculate "
            "a meaningful R-Factor. Stopping..."
        )
        return []
    if for_error and rp.best_v0r is None:
        logger.error(
            "Cannot calculate R-factor for error without a stored"
            "V0r value. Execute a normal R-factor calculation first."
        )
    # Both refcalc (11) and superpos should work fine with new R-factor
    if index == 11:
        name = "refcalc"
    else:
        name = "superpos"
    if (index == 11 and len(rp.refcalc_fdout) == 0) or (
        index == 12 and len(rp.superpos_specout) == 0
    ):
        if index == 11:
            fn = "refcalc-fd.out"
        else:
            fn = "superpos-spec.out"
        if os.path.isfile(os.path.join(".", fn)):
            logger.warning(
                f"R-factor calculation was called without "
                f"stored spectrum data. Reading from file {fn} in work folder..."
            )
            path = os.path.join(".", fn)
        elif os.path.isfile(os.path.join(".", "OUT", fn)):
            logger.warning(
                "R-factor calculation was called without "
                "stored spectrum data. Reading from file " + fn + " in OUT folder..."
            )
            path = os.path.join(".", "OUT", fn)
        else:
            path = ""
            if index == 11:
                # try getting from Tensors
                getTensors(rp.TENSOR_INDEX, required=False)
                dn = "Tensors_" + str(rp.TENSOR_INDEX).zfill(3)
                if os.path.isfile(os.path.join(".", "Tensors", dn, fn)):
                    logger.warning(
                        "R-factor calculation was called without "
                        "stored spectrum data. Reading from file "
                        + fn
                        + " in "
                        + dn
                        + " folder..."
                    )
                    path = os.path.join(".", "Tensors", dn, fn)
            if not path:
                logger.error(
                    "Cannot execute R-factor calculation: no stored "
                    "spectrum data and no " + fn + " file was "
                    "found."
                )
                raise RuntimeError("No spectrum data found")
        try:
            theobeams, theospec = readFdOut(readfile=path)
            if index == 11:
                rp.refcalc_fdout = theospec
            else:
                rp.superpos_specout = theospec
            rp.theobeams[name] = theobeams
        except Exception:
            logger.error("Failed to read " + path)
            raise
        # if we haven't returned, then it was read. check data vs ivbeams:
        eq = True
        eps = 1e-3
        if len(rp.ivbeams) != len(theobeams):
            eq = False
        else:
            eq = all(
                rp.ivbeams[i].isEqual(theobeams[i], eps=eps)
                for i in range(0, len(rp.ivbeams))
            )
        if not eq:
            logger.error(
                "The list of beams read from IVBEAMS is not "
                "equivalent to the list of beams in " + path + ". R-Factor "
                "calculation cannot proceed."
            )
            raise ValueError("Contradiction in beam sets")
    theobeams = rp.theobeams[name]
    expbeams = rp.expbeams

    # Branch off for new R factor calculation
    if not rp.R_FACTOR_LEGACY:
        return run_new_rfactor(sl, rp, for_error, name, theobeams, expbeams)

    return run_legacy_rfactor(sl, rp, for_error, name, theobeams, index, only_vary)


def run_new_rfactor(sl, rp, for_error, name, theobeams, expbeams):
    logger.debug("Using new R-factor calculation. This is still experimental!")
    
    deg, v0i, n_beams, which_r, n_derivs, real_iv_shift = get_general_rfactor_params(rp, for_error)
    
    
    theo_energies = []
    for b in theobeams:
        theo_energies.extend([k for k in b.intens if k not in theo_energies])
    theo_energies.sort()

    exp_energies = get_exp_energies(rp)
    intpol_step = min(
        exp_energies[1] - exp_energies[0], theo_energies[1] - theo_energies[0]
    )

    out_grid= get_rfactor_energy_grid(rp, real_iv_shift, intpol_step)

    if rp.IV_SHIFT_RANGE[2] > 0:
        intpol_step = min(intpol_step, rp.IV_SHIFT_RANGE[2])

    

    # find correspondence experimental to theoretical beams:
    beamcorr = io.getBeamCorrespondence(sl, rp)
    
    # integer & fractional beams
    iorf = []
    for (i, beam) in enumerate(rp.expbeams):
        if beam.hk[0] % 1.0 != 0.0 or beam.hk[1] % 1.0 != 0.0:
            iorf.append(2)
        else:
            iorf.append(1)
    iorf.extend([0] * (len(rp.ivbeams) - len(rp.expbeams)))

    # -1 indicates no correspondence
    n_beams = len(beamcorr) - beamcorr.count(-1)
    exp_grid, exp_id_start, exp_n_E_beams, exp_intensities_in = io.beamlist_to_array(
        rp.expbeams
    )
    corr_theobeams = []
    for i in beamcorr:
        if i != -1:
            corr_theobeams.append(theobeams[i])
    
    # interpolate & get derivative for experimental beams
    (
        exp_e_start_beams_out,
        exp_n_e_beams_out,
        exp_intpol_int,
        exp_intpol_deriv_y
    ) = perform_full_intpol_deriv(deg, n_derivs, out_grid, n_beams, exp_grid, exp_id_start, exp_n_E_beams, exp_intensities_in)
    
    # Same for theoretical beams now
    (
        theo_grid,
        theo_id_start,
        theo_n_E_beams,
        theo_intensities_in,
    ) = io.beamlist_to_array(corr_theobeams)
    
    # interpolate & get derivative for theory (Note we only do this once in rfactor section -- in search this is done more efficiently)
    (
        theo_e_start_beams_out,
        theo_n_e_beams_out,
        theo_intpol_int,
        theo_intpol_deriv_y
    ) = perform_full_intpol_deriv(deg, n_derivs, out_grid, n_beams, theo_grid, theo_id_start, theo_n_E_beams, theo_intensities_in)
    
    # Beams are now interpolated & derivatives known - continue with Y-function & R-factor
    
    
    # For Factor R2 we use intesities for further calculations, otherwise the Y function
    if which_r == 1:
        # Make Y-functions
        rf.pendry_y_beamset(
            exp_intpol_int,
            exp_intpol_deriv_y,  # buffer array with derivative is overwritten by the Y function
            exp_e_start_beams_out,
            exp_n_e_beams_out,
            v0i,
        )

        rf.pendry_y_beamset(
            theo_intpol_int,
            theo_intpol_deriv_y,  # buffer array with derivative is overwritten by the Y function
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
    v0r_range = np.array(
        [int(bound / intpol_step) for bound in rp.IV_SHIFT_RANGE[:2]], dtype="int32"
    )
    v0r_center = int((v0r_range[0] + v0r_range[1]) / 2)
    start_guess = np.array([v0r_range[0], v0r_center, v0r_range[1]], dtype="int32")

    # TODO: make below into User Parameters?
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
        io.writeRfactorPdf_new(
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

def get_rfactor_energy_grid(rp, real_iv_shift, intpol_step):
    
    exp_energies = get_exp_energies(rp)
    
    # extend energy range if they are close together
    if abs(min(exp_energies) - rp.THEO_ENERGIES[0]) < abs(real_iv_shift[0]):
        min_en = max(min(exp_energies), rp.THEO_ENERGIES[0]) + real_iv_shift[0]
    else:
        min_en = max(min(exp_energies), rp.THEO_ENERGIES[0])
    if abs(max(exp_energies) - rp.THEO_ENERGIES[1]) < abs(real_iv_shift[1]):
        max_en = (
            min(max(exp_energies), rp.THEO_ENERGIES[1]) + real_iv_shift[1]
        )  # TODO: should this be + or - ? I think + ...
    else:
        max_en = min(max(exp_energies), rp.THEO_ENERGIES[1])
        
    # generate energy grid
    energy_grid = np.arange(min_en, max_en + intpol_step, intpol_step)
    
    return energy_grid

def get_exp_energies(rp):
    exp_energies = []
    for b in rp.expbeams:
        exp_energies.extend([k for k in b.intens if k not in exp_energies])
    exp_energies.sort()
    return exp_energies


def get_general_rfactor_params(rp, for_error):
    deg = rp.INTPOL_DEG
    v0i = rp.V0_IMAG
    n_beams = len(rp.ivbeams)
    which_r = rp.R_FACTOR_TYPE
    if which_r == 1:
        n_derivs = 1
    elif which_r == 2:
        n_derivs = 0
    else:
        check_ierr(701, logger)
        
    if not for_error:
        real_iv_shift = rp.IV_SHIFT_RANGE[:2]
    else:
        real_iv_shift = [rp.best_v0r] * 2
        
    return deg, v0i, n_beams, which_r, n_derivs, real_iv_shift




def perform_full_intpol_deriv(deg, n_derivs, out_grid, n_beams, grid, id_start, n_E_beams, intensities_in):
    """Performs full (non-optimized i.e slower but more convenient) interpolation and derivative caculation on a beams array.

    Args:
        deg (int): interpolation degree
        n_derivs (_type_): _description_
        out_grid (_type_): _description_
        n_beams (_type_): _description_
        grid (_type_): _description_
        id_start (_type_): _description_
        n_E_beams (_type_): _description_
        intensities_in (_type_): _description_

    Returns:
        _type_: _description_
    """
    # Should we do reordering averaging etc. using the subroutine? - unused right now... # TODO
    averaging_scheme = np.int32(np.arange(n_beams) + 1)  # Fortran index

    # create buffer arrays
    (
        e_start_beams_out,
        n_e_beams_out,
        intpol_int,
        intpol_deriv_y,
        ierr,
    ) = rf.alloc_beams_arrays(n_beams, len(out_grid))
    check_ierr(ierr, logger)

    # interpolation and derivative
    ierr = rf.prepare_beams(
        e_grid_in=grid,
        intensities=intensities_in,
        e_start_beams=id_start + 1,  # Fortran indices
        n_e_beams=n_E_beams,
        deg=deg,
        n_derivs=n_derivs,
        e_grid_out=out_grid,
        e_start_beams_out=e_start_beams_out,
        n_e_beams_out=n_e_beams_out,
        intpol_intensity=intpol_int,
        intpol_derivative=intpol_deriv_y,
    )
    # check for errors
    check_ierr(ierr, logger)
    return e_start_beams_out,n_e_beams_out,intpol_int,intpol_deriv_y


def run_legacy_rfactor(sl, rp, for_error, name, theobeams, index, only_vary):

    if index == 11:
        theospec = rp.refcalc_fdout
    elif index == 12:
        theospec = rp.superpos_specout
    # WEXPEL before PARAM, to make sure number of exp. beams is correct
    try:
        io.writeWEXPEL(sl, rp, theobeams, for_error=for_error)
    except Exception:
        logger.error("Exception during writeWEXPEL: ")
        raise
    try:
        io.writeRfactPARAM(rp, theobeams, for_error=for_error, only_vary=only_vary)
    except Exception:
        logger.error("Exception during writeRfactPARAM: ")
        raise
    # get fortran files and compile
    try:
        tldir = getTLEEDdir(home=rp.sourcedir, version=rp.TL_VERSION)
        if not tldir:
            raise RuntimeError("TensErLEED code not found.")
        libpath = os.path.join(tldir, "lib")
        libname = [f for f in os.listdir(libpath) if f.startswith("rfacsb")][0]
        shutil.copy2(os.path.join(libpath, libname), libname)
        srcpath = os.path.join(tldir, "src")
        srcname = [f for f in os.listdir(srcpath) if f.startswith("rfactor.")][0]
        shutil.copy2(os.path.join(srcpath, srcname), srcname)
    except Exception:
        logger.error("Error getting TensErLEED files for r-factor " "calculation: ")
        raise
    if rp.SUPPRESS_EXECUTION:
        logger.warning(
            "SUPPRESS_EXECUTION parameter is on. R-factor "
            "calculation will not proceed. Stopping..."
        )
        rp.setHaltingLevel(3)
        return []
    logger.info("Compiling fortran input files...")
    rfacname = "rfactor-" + rp.timestamp
    if rp.FORTRAN_COMP[0] == "":
        rp.getFortranComp()
    # compile
    ctasks = [
        (rp.FORTRAN_COMP[0] + " -o " + oname + " -c", fname, rp.FORTRAN_COMP[1])
        for (fname, oname) in [(libname, "rfacsb.o"), (srcname, "main.o")]
    ]
    ctasks.append(
        (rp.FORTRAN_COMP[0] + " -o " + rfacname, "rfacsb.o main.o", rp.FORTRAN_COMP[1])
    )

    # Compile
    compile_log = "compile-rfactor.log"
    try:
        fortran_compile_batch(ctasks, logname=compile_log)
    except Exception:
        logger.error("Error compiling fortran files: ", exc_info=True)
        raise
    else:
        # move log file to supp
        try:
            shutil.move(compile_log, os.path.join("compile_logs", compile_log))
        except OSError:
            logger.warning("Could not move compile-rfactor.log to SUPP")
    # run
    rfaclogname = rfacname + ".log"
    logger.info(
        "Starting R-factor calculation...\n"
        "R-factor log will be written to file " + rfaclogname
    )
    try:
        with open(rfaclogname, "w") as log:
            subprocess.run(
                os.path.join(".", rfacname),
                input=theospec,
                encoding="ascii",
                stdout=log,
                stderr=log,
            )
    except Exception:
        logger.error(
            "Error during R-factor calculation. Also check " "R-factor log file."
        )
        raise
    logger.info("Finished R-factor calculation. Processing files...")
    # rename and move files
    try:
        os.rename("WEXPEL", "rfactor-WEXPEL")
    except OSError:
        logger.warning(
            "Failed to rename R-factor input file WEXPEL to " "rfactor-WEXPEL"
        )
    try:
        os.rename("PARAM", "rfactor-PARAM")
    except OSError:
        logger.warning("Failed to rename R-factor input file PARAM to " "rfactor-PARAM")
    if not os.path.isfile(os.path.join(".", "ROUT")):
        logger.error("No ROUT file was found after R-Factor calculation!")
        rp.setHaltingLevel(2)
        return []

    # read output
    if for_error:
        try:
            rfaclist = io.readROUTSHORT()
        except Exception:
            logger.error("Error reading ROUTSHORT file", exc_info=rp.LOG_DEBUG)
            rp.setHaltingLevel(2)
        return rfaclist

    try:
        (rfac, r_int, r_frac), v0rshift, rfaclist = io.readROUT()
    except Exception:  # TODO catch correct exception
        logger.error("Error reading ROUT file", exc_info=rp.LOG_DEBUG)
        rp.setHaltingLevel(2)
    else:
        logger.info(
            "With inner potential shift of {:.2f} eV: "
            "R = {:.4f}\n".format(v0rshift, rfac)
        )
        rp.best_v0r = v0rshift
        dir_list = ["."]
        if os.path.isdir("OUT"):
            dir_list.append("OUT")
        for dir_name in dir_list:
            for f_name in [
                f
                for f in os.listdir(dir_name)
                if f.startswith("R_OUT_" + rp.timestamp)
                and os.path.isfile(os.path.join(dir_name, f))
            ]:
                try:  # delete old R_OUT files
                    os.remove(os.path.join(dir_name, f_name))
                except Exception:
                    pass
        if rfac == 0:
            logger.error(
                "ROUT reports R-Factor as zero. This means "
                "something went wrong in the reference "
                "calculation or in the R-factor calculation."
            )
            rp.setHaltingLevel(2)
            f_name = "R_OUT_" + name + "_" + rp.timestamp
        else:
            f_name = "R_OUT_" + name + "_" + rp.timestamp + "_R={:.4f}".format(rfac)
            rp.last_R = rfac
            rp.stored_R[name] = (rfac, r_int, r_frac)
        try:
            os.rename("ROUT", f_name)
        except OSError:
            logger.warning("Failed to rename R-factor output file " "ROUT to " + f_name)
        if len(rfaclist) != len(rp.expbeams):
            logger.warning(
                "Failed to read R-Factors per beam from " "R-factor output file ROUT."
            )
            rfaclist = [-1] * len(rp.expbeams)

        if rp.PLOT_IV["plot"]:
            outname = "Rfactor_plots_{}.pdf".format(name)
            aname = "Rfactor_analysis_{}.pdf".format(name)
            if rp.PLOT_IV["overbar"]:
                labelstyle = "overbar"
            else:
                labelstyle = "minus"
            labelwidth = max(
                [beam.getLabel(style=labelstyle)[1] for beam in rp.expbeams]
            )
            try:
                io.writeRfactorPdf(
                    [
                        (
                            b.getLabel(lwidth=labelwidth, style=labelstyle)[0],
                            rfaclist[i],
                        )
                        for (i, b) in enumerate(rp.expbeams)
                    ],
                    outName=outname,
                    analysisFile=aname,
                    v0i=rp.V0_IMAG,
                    formatting=rp.PLOT_IV,
                )
            except Exception:  # TODO catch correct excpetion
                logger.warning("Error plotting R-factors.", exc_info=True)
    return rfaclist
