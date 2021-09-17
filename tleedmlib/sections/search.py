# -*- coding: utf-8 -*-

"""
Created on Aug 11 2020

@author: Florian Kraushofer

Tensor LEED Manager section Search
"""

import os
import sys
import logging
import shutil
import subprocess
import time
from timeit import default_timer as timer
import numpy as np
import signal
import re
import copy
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
import scipy

import viperleed.tleedmlib.files.iosearch as io
import viperleed.tleedmlib as tl
# from tleedmlib.polynomialfeatures_no_interaction import PolyFeatNoMix
from viperleed.tleedmlib.leedbase import fortran_compile_batch
from viperleed.tleedmlib.files.parameters import updatePARAMETERS
from viperleed.tleedmlib.files.displacements import readDISPLACEMENTS_block
from viperleed.tleedmlib.files.searchpdf import (
    writeSearchProgressPdf, writeSearchReportPdf)

logger = logging.getLogger("tleedm.search")


def processSearchResults(sl, rp, final=True):
    """
    Looks at the final block in the SD.TL file and gets the data from the
    best structure in into the slab.

    Parameters
    ----------
    sl : Slab
        The Slab object to be modified
    rp : Rparams
        Run parameters.
    final : bool, optional
        Defines whether this is the final interpretation that *must* work and
        should go into the log, or if this will be repeated anyway. The
        default is True.

    Returns
    -------
    None

    """
    # get the last block from SD.TL:
    try:
        lines = io.readSDTL_end()
    except FileNotFoundError:
        logger.error("Could not process Search results: SD.TL file not "
                     "found.")
        raise
    except Exception:
        logger.error("Failed to get last block from SD.TL file.")
        raise
    # read the block
    sdtlContent = io.readSDTL_blocks("\n".join(lines),
                                     whichR=rp.SEARCH_BEAMS, logInfo=final)
    if not sdtlContent:
        if final:
            logger.error("No data found in SD.TL file!")
            rp.setHaltingLevel(2)
        raise RuntimeError("No data in SD.TL")
    (generation, rfacs, configs) = sdtlContent[0]
    # collect equal entries
    pops = []  # tuples (rfactor, list of parameter indices)
    popcount = []   # how often a given population is there
    dparlist = []   # for control.chem
    for i in range(0, len(rfacs)):
        if (rfacs[i], configs[i]) in pops:
            popcount[pops.index((rfacs[i], configs[i]))] += 1
        else:
            pops.append((rfacs[i], configs[i]))
            popcount.append(1)
        dparlist.append(configs[i])
    writeControlChem = False
    if not os.path.isfile("control.chem"):
        writeControlChem = True
    else:
        # check file, which generation
        try:
            with open("control.chem", "r") as rf:
                rf.readline()  # first line is empty.
                s = rf.readline()
                s = s.split("No.")[-1].split(":")[0]
                if int(s) != generation:
                    writeControlChem = True
                if len(rf.readlines()) < rp.SEARCH_POPULATION:
                    writeControlChem = True
        except Exception:
            # try overwriting, just to be safe
            writeControlChem = True
    # write backup file
    output = ("\nParameters of generation No." + str(generation).rjust(6)
              + ":\n")
    if rp.domainParams:
        astep = rp.DOMAIN_STEP
    else:
        astep = 100
    if rp.TL_VERSION < 1.7:
        ctrl_width = 3
    else:
        ctrl_width = 4
    for dpars in dparlist:
        ol = ""
        areapars = []
        for (percent, l) in dpars:
            for ind in l:
                ol += str(ind).rjust(ctrl_width)
            areapars.append(int(percent/astep)+1)
        for ap in areapars:
            ol += str(ap).rjust(ctrl_width)
        output += ol + "\n"
    rp.controlChemBackup = output
    if writeControlChem:
        try:
            with open("control.chem", "w") as wf:
                wf.write(output)
        except Exception:
            logger.error("Failed to write control.chem")
            rp.setHaltingLevel(1)

    # info for log:
    maxpop = max(popcount)
    if final:
        if len(pops) == 1:
            logger.info("All trial structures converged to the same "
                        "configuration (R = {:.4f})".format(pops[0][0]))
        elif any([v > rp.SEARCH_POPULATION/2 for v in popcount]):
            info = ("The search outcome is dominated by one configuration "
                    "(population {} of {}, R = {:.4f})\n".format(
                        maxpop, rp.SEARCH_POPULATION,
                        pops[popcount.index(maxpop)][0]))
        else:
            info = ("The search outcome is not dominated by any "
                    "configuration. The best result has population {} (of {})"
                    ", R = {:.4f}\n".format(popcount[0], rp.SEARCH_POPULATION,
                                            pops[0][0]))
        if len(pops) > 1:
            info += ("The best configurations are:\nPOP       R |")
            if rp.domainParams:
                info += " area |"
            info += " PARAMETERS\n"
            for i in range(0, min(5, len(pops))):
                for (j, (percent, pars)) in enumerate(pops[i][1]):
                    if j == 0:
                        info += "{:>3}  {:.4f} |".format(popcount[i],
                                                         pops[i][0])
                    else:
                        info += "            |"
                    if rp.domainParams:
                        info += " {:>3}% |".format(percent)
                    for v in pars:
                        info += "{:>3}".format(v)
                    info += "\n"
            logger.info(info)
    # now writeSearchOutput:
    if not rp.domainParams:
        # the order here matters, as the final write operation will store the
        #  new states (for later searches)
        if any(sp.parabolaFit["min"] for sp in rp.searchpars):
            parab_inds = list(pops[0][1][0][1])
            for i, sp in enumerate(rp.searchpars):
                if sp.parabolaFit["min"] is not None:
                    parab_inds[i] = sp.parabolaFit["min"]
            io.writeSearchOutput(sl, rp, parab_inds, silent=True,
                                 suffix="_parabola")
        io.writeSearchOutput(sl, rp, pops[0][1][0][1], silent=(not final))
    else:
        home = os.getcwd()
        for (i, dp) in enumerate(rp.domainParams):
            try:
                os.chdir(dp.workdir)
                if any(sp.parabolaFit["min"] for sp in rp.searchpars):
                    parab_inds = list(pops[0][1][i][1])
                    for j, sp in enumerate(dp.rp.searchpars):
                        if sp.parabolaFit["min"] is not None:
                            parab_inds[j] = sp.parabolaFit["min"]
                    io.writeSearchOutput(dp.sl, dp.rp, parab_inds, silent=True,
                                         suffix="_parabola")
                io.writeSearchOutput(dp.sl, dp.rp, pops[0][1][i][1],
                                     silent=(not final))
            except Exception:
                logger.error("Error while writing search output for domain {}"
                             .format(dp.name), exc_info=rp.LOG_DEBUG)
                rp.setHaltingLevel(2)
                raise
            finally:
                os.chdir(home)
    return


def parabolaFit(rp, datafiles, r_best, x0=None, max_configs=0, **kwargs):
    """
    Performs a parabola fit to all rfactor/configuration data in r_configs.
    Stores the result in the rp searchpar objects. Keyword arguments are read
    from rp by default.

    Parameters
    ----------
    rp : Rparams
        The Rparams object containing runtime information
    datafiles : list of str
        The files to read data from
    r_best : float
        Best known R-factor so far
    x0 : numpy.array, optional
        A starting guess for the minimum
    max_configs : int
        0 to read all, else specifies the maximum number of configurations to
        read in. Will sort configurations and read in the ones with lowest
        R-factors first, discard the rest.
    Keyword Arguments:
        localize : real, optional
            If not zero, r_configs will be reduced to only use points close to
            the current best result. The value determines the region to be used
            ascompared to the full data, i.e. 0.5 would use only half of the
            space in each dimesion.
        mincurv : float, optional
            The minimum curvature that the N-dimensional paraboloid is
            required to have along the axis of a given parameter, in order for
            the minimum and curvature to be saved for that parameter.
        type : str, optional
            Defines the regression model to use in the fit. Allowed values are
            LinearRegression (default), Ridge, Lasso, and ElasticNet. Note that
            for the fit, parameters are centered around the current best
            configuration. Therefore, all models except LinearRegression put
            some penalty not only on curvature, but also on the parameter
            offset from the best-R configuration.
            Setting 'type' to 'none' turns off parabola fit entirely.
        alpha : float, optional
            The alpha value used by Ridge, Lasso and ElasticNet regression. The
            default is 1e-5.

    Returns
    -------
    np.array
        The result vector in parameter space, consisting of the coordinates of
        the parabola minimum for those dimensions where the parabola curvature
        is greater than mincurv, and the values from the configuration with
        the lowest R-factor for all dimensions. Can be passed back as x0 in
        the next iteration.
    float
        Predicted minimum R-factor.

    """

    def optimizerHelper(array, func):
        return func(array.reshape(1, -1))

    def castToMatrix(features, dim):
        """Casts the weights of PolynomialFeatures of degree 2 into a square
        symmetric 2D array, discarding linear features."""
        m = np.zeros((dim, dim))
        start = dim+1
        for i in range(0, dim):
            m[i][i:] = features[start:start + dim - i]
            start += dim - i
        m = 0.5*(m + m.T)
        return m

    if not datafiles:
        return None, None
    if "type" not in kwargs:
        which_regression = rp.PARABOLA_FIT["type"]
    else:
        which_regression = kwargs["type"]
    if which_regression == "none":
        return None, None
    for a in ("localize", "mincurv", "alpha", "type"):
        if a not in kwargs:
            kwargs[a] = rp.PARABOLA_FIT[a]

    # read data
    r_cutoff = 1.0
    if datafiles:
        varR = np.sqrt(8*np.abs(rp.V0_IMAG) / rp.total_energy_range())*r_best
        rc = np.array((io.readDataChem(
            rp, datafiles,
            cutoff=r_best + r_cutoff * varR,
            max_configs=max_configs)), dtype=object)
    if len(rc) < 100*rp.indyPars:
        return None, None

    rfacs, configs = rc[:, 0].astype(float), rc[:, 1]
    localizeFactor = kwargs["localize"]
    if localizeFactor == 0:
        localizeFactor = 1  # no localization
    # reduce to independent parameters during fit
    if not rp.domainParams:
        sps = [sp for sp in rp.searchpars if sp.el != "vac" and
               sp.mode != "dom" and sp.steps*localizeFactor >= 3 and
               sp.linkedTo is None and sp.restrictTo is None]
        indep_pars = np.array([*np.array([*np.array(configs, dtype=object)],
                                         dtype=object)
                               .reshape(-1, 1, 2)[:, 0, 1]])
        indep_pars = np.delete(indep_pars, [i for i in
                                            range(len(rp.searchpars)-1)
                                            if rp.searchpars[i] not in sps], 1)
    else:
        # first the percentages:
        sps = [sp for sp in rp.searchpars if sp.mode == "dom"]
        reshaped = (np.array([*np.array(configs, dtype=object)], dtype=object)
                    .reshape(-1, len(rp.domainParams), 2))
        indep_pars = reshaped[:, :, 0].astype(int)  # contains the percentages
        # then the 'real' parameters:
        for (j, dp) in enumerate(rp.domainParams):
            new_sps = [sp for sp in dp.rp.searchpars if
                       sp.el != "vac" and sp.steps*localizeFactor >= 3 and
                       sp.linkedTo is None and sp.restrictTo is None and
                       sp.mode != "dom"]
            new_ip = np.array([*reshaped[:, j, 1]])
            new_ip = np.delete(new_ip, [i for i in range(len(dp.rp.searchpars))
                               if dp.rp.searchpars[i] not in new_sps], 1)
            sps.extend(new_sps)
            indep_pars = np.append(indep_pars, new_ip, axis=1)
    indep_pars = indep_pars.astype(float)
    best_config = np.copy(indep_pars[np.argmin(rfacs)])
    sps_original = sps[:]

    which_regression = kwargs["type"].lower()
    alpha = kwargs["alpha"]
    if which_regression == 'lasso':
        polyreg = make_pipeline(PolynomialFeatures(degree=2),
                                Lasso(alpha=alpha, normalize=True))
    elif which_regression == 'ridge':
        polyreg = make_pipeline(PolynomialFeatures(degree=2),
                                Ridge(alpha=alpha, normalize=True))
    elif which_regression == 'elasticnet':
        polyreg = make_pipeline(PolynomialFeatures(degree=2),
                                ElasticNet(alpha=alpha, normalize=True))
    else:
        if which_regression not in ('linearregression', 'linear'):
            logger.warning("Regression model {} not found, parabola fit "
                           "defaulting to linear regression."
                           .format(which_regression))
        which_regression = "linearregression"
        polyreg = make_pipeline(PolynomialFeatures(degree=2),
                                LinearRegression())
    xmin = np.copy(indep_pars[np.argmin(rfacs), :])
    rr = np.sqrt(8*np.abs(rp.V0_IMAG) / rp.total_energy_range())
    ip_tofit = np.copy(indep_pars)
    rf_tofit = np.copy(rfacs)
    # throw out high R-factors - TODO: perhaps also throw out highest X% ?
    to_del = np.where(rf_tofit > min(rf_tofit) + 3*rr)
    ip_tofit = np.delete(ip_tofit, to_del, axis=0)
    rf_tofit = np.delete(rf_tofit, to_del, axis=0)
    # center on best config; important because offset from the
    #  best configuration may also be penalized in regression
    ip_tofit = ip_tofit - xmin
    # fit
    polyreg.fit(ip_tofit, rf_tofit)
    # now find minimum within bounds (still centered around xmin)
    bounds = np.array([(1, sp.steps) for sp in sps]) - xmin.reshape((-1, 1))
    if x0 is None:  # x0 is the starting guess
        x0 = np.copy(xmin)
    # else:
    #     x0 = np.delete(x0, deletedPars)
    parab_min = scipy.optimize.minimize(
        optimizerHelper, x0, args=(polyreg.predict,), method='L-BFGS-B',
        bounds=bounds)
    predictR = parab_min.fun[0]
    # find errors
    m = castToMatrix(polyreg.named_steps[which_regression].coef_, len(sps))
    w, v = np.linalg.eig(m)
    # error along main axes
    d = np.copy(np.diag(m))
    d[d <= 0] = np.nan
    with np.errstate(invalid="ignore"):
        err_unco = 2 * np.sqrt(rr * predictR / d)
        err_unco[err_unco < 0] = np.nan
    err_unco = np.sqrt(err_unco)
    # error along eigenvectors
    w2 = np.copy(w)
    w2[w2 == 0] = 1e-100  # to avoid divide-by-zero; dealt with below
    with np.errstate(invalid="ignore"):
        err_ev = 2 * np.sqrt(rr * predictR / w2)
        err_ev[err_ev < 0] = 0  # dealt with below
    err_ev = np.sqrt(err_ev)
    # correlated error
    err_co = np.dot(v**2, err_ev)
    for i in range(len(err_co)):
        if any(w[j] <= 0 and v[i, j] >= 0.1 for j in range(len(w))):
            err_co[i] = np.nan
    # !!! TODO: maybe discard parameters and re-fit?

    # deletedPars = []
    # while True:
    #     ip_tofit = np.copy(indep_pars)
    #     rf_tofit = np.copy(rfacs)
    #     xmin = np.copy(ip_tofit[np.argmin(rfacs), :])
    #     if kwargs["localize"] != 0:
    #         # discard points that are far from global min in any dimension
    #         base = np.copy(xmin)  # parameter vector of best conf
    #         for i in range(len(base)):
    #             r = sps[i].steps*localizeFactor*0.5
    #             base[i] = max(base[i], 1 + r)
    #             base[i] = min(base[i], sps[i].steps - r)
    #         dist_norm = np.array([1/(sp.steps-1) for sp in sps])
    #         dist = np.abs(dist_norm*(ip_tofit - base))
    #         maxdist = np.max(dist, 1)
    #         # dist = np.linalg.norm(dist_norm*(indep_pars - base), axis=1)
    #         # cutoff = np.percentile(dist, 30)
    #         # cutoff = np.max(dist) * 0.75
    #         cutoff = 0.5*localizeFactor
    #         # weights = [1 if d <= cutoff else 0.1 for d in dist]
    #         to_del = [i for i in range(len(maxdist)) if maxdist[i] > cutoff]
    #         ip_tofit = np.delete(ip_tofit, to_del, axis=0)
    #         rf_tofit = np.delete(rf_tofit, to_del, axis=0)
    #     # center on best config; important because offset from the
    #     #  best configuration may also be penalized in regression
    #     ip_tofit = ip_tofit - xmin
    #     # check curvature of the parabolas per parameter
    #     polyreg.fit(ip_tofit, rf_tofit)
    #     curvs = [(polyreg.named_steps[which_regression]
    #               .coef_[polyreg.named_steps['polynomialfeatures']
    #               .get_feature_names().index('x{}^2'.format(i))])
    #              * sp.steps**2
    #              for i, sp in enumerate(sps)]
    #     if min(curvs) >= kwargs["mincurv"]:
    #         if len(rf_tofit) < 50*len(sps):  # too few points for good fit
    #             # logger.debug("{} points - too few for fit of {} dimensions"
    #             #               .format(len(rf_tofit), len(sps)))
    #             for sp in rp.searchpars:
    #                 sp.parabolaFit = {"curv": None, "min": None}
    #             return None, None
    #         break  # all parabola dimensions now have reasonable shape
    #     if len(sps) == 1:
    #         # logger.debug("Found no parameter with sufficient curvature")
    #         sps[0].parabolaFit = {"curv": None, "min": None}
    #         return None, None  # found no parameter with parabola shape
    #     # throw out the parameters with lowest curvature; repeat.
    #     #   discard: min. 1, max. all with low curv.
    #     k = min(max(1, int(len(curvs)*0.2)),
    #             sum(c < kwargs["mincurv"] for c in curvs))
    #     ind = np.argpartition(curvs, k)[:k]
    #     for i in ind:
    #         sps[i].parabolaFit = {"curv": None, "min": None}
    #         deletedPars.append(sps_original.index(sps[i]))
    #     sps = [sps[i] for i in range(len(sps)) if i not in ind]
    #     indep_pars = np.delete(indep_pars, ind, axis=1)
    # for (i, sp) in enumerate(sps):
    #     sp.parabolaFit["curv"] = curvs[i]

    # logger.debug("Parabola fit: {}/{} parameters were fit with {} points "
    #              "({:.4f} s)".format(len(sps), len(sps_original),
    #                                  len(rf_tofit), (timer() - starttime)))

    for (i, sp) in enumerate(sps):
        sp.parabolaFit["min"] = parab_min.x[i] + xmin[i]
        sp.parabolaFit["err_co"] = err_co[i]
        sp.parabolaFit["err_unco"] = err_unco[i]
        best_config[sps_original.index(sp)] = parab_min.x[i] + xmin[i]
    for sp in [sp for sp in rp.searchpars if type(sp.linkedTo) == tl.SearchPar
               or type(sp.restrictTo) == tl.SearchPar]:
        if type(sp.linkedTo) == tl.SearchPar:
            sp.parabolaFit = copy.copy(sp.linkedTo.parabolaFit)
        elif type(sp.restrictTo) == tl.SearchPar:
            sp.parabolaFit = copy.copy(sp.restrictTo.parabolaFit)
    return (best_config, predictR)


def search(sl, rp):
    """
    Generates input for and runs the search, then interprets output.

    Parameters
    ----------
    sl : Slab
        Slab object containing atom, site and layer information.
    rp : Rparams
        The run parameters.

    Raises
    ------
    RuntimeError
        Raised if execution cannot continue.

    Returns
    -------
    None.

    """

    def kill_process(proc, default_pgid=None):
        """Cleanly kill the mpirun subprocess and its children. If the process
        is not alive any more (and therefore the pgid cannot be determined),
        will instead try to terminate the default_pgid, if passed."""
        # determine pgid
        try:
            pgid = os.getpgid(proc.pid)
        except ProcessLookupError:
            pgid = default_pgid
        # kill main process
        try:
            proc.kill()
            proc.wait()
        except ProcessLookupError:
            pass   # already dead
        if pgid is not None:
            # kill children
            try:
                os.killpg(pgid, signal.SIGTERM)
                if os.name == "nt":
                    os.waitpid(pgid)
                else:
                    os.waitpid(-pgid, 0)
            except (ProcessLookupError, ChildProcessError):
                pass  # already dead or no children
        return

    rp.searchResultConfig = None
    if rp.domainParams:
        initToDo = [(dp.rp, dp.sl, dp.workdir) for dp in rp.domainParams]
    else:
        initToDo = [(rp, sl, ".")]
    for (rpt, slt, path) in initToDo:
        # read DISPLACEMENTS block
        if not rpt.disp_block_read:
            readDISPLACEMENTS_block(rpt, slt,
                                    rpt.disp_blocks[rpt.search_index])
            rpt.disp_block_read = True
        # get Deltas
        if 2 not in rpt.runHistory:
            if "Tensors" in rpt.manifest:
                logger.error("New tensors were calculated, but no new delta "
                             "files were generated. Cannot execute search.")
                raise RuntimeError("Delta calculations was not run for "
                                   "current tensors.")
            tl.leedbase.getDeltas(rpt.TENSOR_INDEX, basedir=path,
                                  targetdir=path, required=True)
    rp.updateCores()
    # generate rf.info
    try:
        rfinfo = io.writeRfInfo(sl, rp, filename="rf.info")
    except Exception:
        logger.error("Error generating search input file rf.info")
        raise
    # generate PARAM and search.steu
    #   needs to go AFTER rf.info, as writeRfInfo may remove expbeams!
    try:
        io.generateSearchInput(sl, rp)
    except Exception:
        logger.error("Error generating search input")
        raise
    if rp.indyPars == 0:  # never for calculations with domains # !!! CHECK
        logger.info("Found nothing to vary in search. Will proceed "
                    "directly to writing output and starting SUPERPOS.")
        rp.searchResultConfig = [(100, [1] * (len(rp.searchpars))-1)]
        for (i, sp) in enumerate(rp.searchpars):
            if type(sp.restrictTo) == int:
                rp.searchResultsConfig[i] = sp.restrictTo
            elif type(sp.restrictTo) == tl.SearchPar:
                rp.searchResultsConfig[i] = (rp.searchpars.index(
                                                    sp.restrictTo) + 1)
            elif type(sp.linkedTo) == tl.SearchPar:
                rp.searchResultsConfig[i] = (rp.searchpars.index(
                                                    sp.linkedTo) + 1)
        io.writeSearchOutput(sl, rp)
        return None
    if rp.SUPPRESS_EXECUTION:
        logger.warning("SUPPRESS_EXECUTION parameter is on. Search "
                       "will not proceed. Stopping...")
        rp.setHaltingLevel(3)
        return None
    # check for mpirun, decide whether to use parallelization
    usempi = True
    if rp.N_CORES == 1:
        logger.warning(
            "The N_CORES parameter is set to 1. The search will be run "
            "without multiprocessing. This will be much slower!")
        usempi = False

    if usempi and shutil.which("mpirun", os.X_OK) is None:
        usempi = False
        logger.warning(
            "mpirun is not present. Search will be compiled and executed "
            "without parallelization. This will be much slower!")
        if rp.FORTRAN_COMP[0] == "":
            try:
                rp.getFortranComp()
            except Exception:
                logger.error("No fortran compiler found, cancelling...")
                raise RuntimeError("Fortran compile error")
    if usempi:
        if rp.FORTRAN_COMP_MPI[0] == "":
            try:
                rp.getFortranMpiComp()
            except Exception:
                logger.error("No fortran mpi compiler found, cancelling...")
                raise RuntimeError("Fortran compile error")
    else:
        if rp.FORTRAN_COMP[0] == "":
            try:
                rp.getFortranComp()
            except Exception:
                logger.error("No fortran compiler found, cancelling...")
                raise RuntimeError("Fortran compile error")
    # get fortran files
    try:
        tldir = tl.leedbase.getTLEEDdir(home=rp.sourcedir,
                                        version=rp.TL_VERSION)
        if not tldir:
            raise RuntimeError("TensErLEED code not found.")
        srcpath = os.path.join(tldir, 'src')
        if usempi:
            srcname = [f for f in os.listdir(srcpath)
                       if f.startswith('search.mpi')][0]
        else:
            srcname = [f for f in os.listdir(srcpath)
                       if f.startswith('search') and 'mpi' not in f][0]
        shutil.copy2(os.path.join(srcpath, srcname), srcname)
        libpath = os.path.join(tldir, 'lib')
        if usempi:
            libname = [f for f in os.listdir(libpath)
                       if f.startswith('lib.search.mpi')][0]
        else:
            libname = [f for f in os.listdir(libpath)
                       if f.startswith('lib.search') and 'mpi' not in f][0]
        shutil.copy2(os.path.join(libpath, libname), libname)
        hashing_files = [f for f in os.listdir(libpath)
                         if f.startswith('intarr_hashing')
                         and f.endswith('.f90')]
        if hashing_files:
            hashname = hashing_files[0]
            shutil.copy2(os.path.join(libpath, hashname), hashname)
        else:
            hashname = ""
        if usempi:  # these are short C scripts - use pre-compiled versions
            randnamefrom = "MPIrandom_.o"
        else:
            randnamefrom = "random_.o"
        randname = "random_.o"
        shutil.copy2(os.path.join(libpath, randnamefrom), randname)
        globalname = "GLOBAL"
        shutil.copy2(os.path.join(srcpath, globalname), globalname)
    except Exception:
        logger.error("Error getting TensErLEED files for search: ")
        raise
    # compile fortran files
    searchname = "search-"+rp.timestamp
    if usempi:
        fcomp = rp.FORTRAN_COMP_MPI
    else:
        fcomp = rp.FORTRAN_COMP
    logger.info("Compiling fortran input files...")
    # compile
    ctasks = [(fcomp[0]+" -o lib.search.o -c", libname, fcomp[1])]
    if hashname:
        ctasks.append((fcomp[0]+" -c", hashname, fcomp[1]))
    ctasks.append((fcomp[0]+" -o restrict.o -c", "restrict.f", fcomp[1]))
    format_tag = ""
    if any([f.endswith('.f90') for f in (libname, srcname, hashname)]):
        format_tag = "-fixed"
        if any([s in fcomp[0] for s in ("gfortran", "mpifort")]):
            # assume that mpifort also uses gfortran
            format_tag = "--fixed-form"  #  different formatting string
    ctasks.append((fcomp[0]+" -o search.o -c "+format_tag, srcname, fcomp[1]))
    to_link = "search.o random_.o lib.search.o restrict.o"
    if hashname:
        to_link += " intarr_hashing.o"
    ctasks.append((fcomp[0] + " -o " + searchname, to_link, fcomp[1]))
    try:
        fortran_compile_batch(ctasks)
    except Exception:
        logger.error("Error compiling fortran files: ", exc_info=True)
        raise
    logger.debug("Compiled fortran files successfully")
    # run
    if rp.LOG_SEARCH:
        searchlogname = searchname+".log"
        logger.info("Search log will be written to file "+searchlogname)
        if rp.TL_VERSION > 1.6:
            rp.manifest.append(searchlogname)
    # if there is an old SD.TL file, it needs to be removed
    if os.path.isfile("SD.TL"):
        try:
            os.remove("SD.TL")
        except Exception:
            logger.warning("Failed to delete old SD.TL file. This may "
                           "cause errors in the interpretation of search "
                           "progress.")
    # same for old data.chem
    for fn in [f for f in os.listdir() if re.match(r'data\d+\.chem$',
                                                   f.lower())]:
        try:
            os.remove(fn)
        except Exception:
            logger.warning("Failed to delete old {} file. This may cause "
                           "errors in the parabola fit.".format(fn))
    # get config size
    config_size = (
        sys.getsizeof((1., tuple()))
        + sys.getsizeof(tuple([tuple()] * max(1, len(rp.domainParams))))
        + sys.getsizeof((1, tuple())))
    if rp.domainParams:
        config_size += sum(
            [sys.getsizeof(tuple([0] * len(dp.rp.searchpars)))
             + sys.getsizeof(1) * len(dp.rp.searchpars)
             for dp in rp.domainParams])
    else:
        config_size += (sys.getsizeof(tuple([0] * len(rp.searchpars)))
                        + sys.getsizeof(1) * len(rp.searchpars))
    max_memory_parab = 1e9   # maximum memory in bytes to use on parabola fit
    max_read_configs = int(max_memory_parab / config_size)
    # start search process
    repeat = True
    first = True
    genOffset = 0
    gens = []  # generation numbers in SD.TL, but continuous if search restarts
    markers = []
    rfaclist = []
    parab_x0 = None     # starting guess for parabola
    rfac_predict = []  # tuples (gen, r) from parabola fit
    realLastConfig = {"all": [], "best": [], "dec": []}
    realLastConfigGen = {"all": 0, "best": 0, "dec": 0}
    convergedConfig = {"all": None, "best": None, "dec": None}
    lastconfig = None
    rp.searchMaxGenInit = rp.SEARCH_MAX_GEN
    absstarttime = timer()
    tried_repeat = False        # if SD.TL is not written, try restarting
    pgid = None
    while repeat:
        if first:
            logger.info("Starting search. See files Search-progress.pdf "
                        "and SD.TL for progress information.")
            first = False
        repeat = False
        interrupted = False
        proc = None
        if usempi:
            command = ["mpirun", "-n", str(rp.N_CORES),
                       os.path.join(".", searchname)]
        else:
            command = os.path.join('.', searchname)
        try:
            if not rp.LOG_SEARCH:
                proc = subprocess.Popen(
                    command, encoding="ascii",
                    stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT,
                    preexec_fn=os.setsid)
            else:
                logExists = os.path.isfile(searchlogname)
                with open(searchlogname, "a") as log:
                    if logExists:
                        log.write("\n\n-------\nRESTARTING\n-------\n\n")
                    proc = subprocess.Popen(
                        command, encoding="ascii",
                        stdout=log, stderr=log,
                        preexec_fn=os.setsid)
            pgid = os.getpgid(proc.pid)
        except Exception:
            logger.error("Error starting search. Check SD.TL file.")
            raise
        if proc is None:
            logger.error("Error starting search subprocess... Stopping.")
            raise RuntimeError("Error running search")
        # FEED INPUT
        try:
            proc.communicate(input=rfinfo, timeout=0.1)
        except subprocess.TimeoutExpired:
            pass  # started successfully; monitoring below
        except Exception:
            logger.error("Error starting search. Check SD.TL file.")
        # MONITOR SEARCH
        searchStartTime = timer()
        printt = searchStartTime
        filepos = 0
        timestep = 1  # time step to check files
        # !!! evaluation time could be higher - keep low only for debugging
        evaluationTime = 60  # how often should SD.TL be evaluated
        lastEval = 0  # last evaluation time (s), counting from searchStartTime
        comment = ""
        sdtlGenNum = 0
        gaussianWidthOri = rp.GAUSSIAN_WIDTH
        check_datafiles = False
        try:
            while proc.poll() is None:
                time.sleep(timestep)
                updatePARAMETERS(rp)
                # check convergence criteria
                stop = False
                checkrepeat = True
                if rp.STOP is True:
                    stop = True
                    checkrepeat = False
                    logger.info("Search stopped by STOP command.")
                    if not os.path.isfile("SD.TL"):
                        # try saving by waiting for SD.TL to be created...
                        logger.warning("SD.TL file not found. Trying to "
                                       "wait, maximum 5 minutes...")
                        i = 0
                        while not os.path.isfile("SD.TL") and not i >= 300:
                            time.sleep(1)
                            i += 1
                if rp.GAUSSIAN_WIDTH != gaussianWidthOri:
                    stop = True
                    repeat = True
                    comment = ("GAUSSIAN_WIDTH = {}"
                               .format(rp.GAUSSIAN_WIDTH))
                    logger.info("GAUSSIAN_WIDTH parameter changed. "
                                "Search will restart.")
                t = timer() - searchStartTime
                if (t - lastEval > evaluationTime) or stop:
                    # evaluate
                    lastEval = t
                    newData = []
                    if os.path.isfile("SD.TL"):
                        filepos, content = io.readSDTL_next(offset=filepos)
                        if content:
                            newData = io.readSDTL_blocks(
                                content, whichR=rp.SEARCH_BEAMS)
                    elif t >= 900 and rp.HALTING < 3:
                        stop = True
                        if tried_repeat:
                            logger.warning(
                                "No SD.TL file was written for 15 minutes "
                                "after restarting the search. Search will "
                                "stop.")
                            repeat = False
                            rp.setHaltingLevel(2)
                        else:
                            repeat = True
                            tried_repeat = True
                            logger.warning(
                                "No SD.TL file was written for 15 minutes "
                                "after the search started. Trying to restart. "
                                "You can suppress this behaviour by setting "
                                "the HALTING parameter to 3.")
                    for (gen, rfacs, configs) in newData:
                        gens.append(gen + genOffset)
                        sdtlGenNum = gen
                        rfaclist.append(np.array(rfacs))
                        if gen % 1000 == 0:
                            speed = 1000*(timer() - absstarttime)/gens[-1]
                            logger.debug(
                                "R = {:.4f} (Generation {}, {:.1f} s "
                                "since last, {:.1f} s/kG overall)".format(
                                    min(rfacs), gens[-1], timer() - printt,
                                    speed))
                            printt = timer()
                            check_datafiles = True
                        dgen = {}
                        for k in ["dec", "best", "all"]:
                            dgen[k] = gens[-1] - realLastConfigGen[k]
                        if configs != realLastConfig["all"]:
                            realLastConfig["all"] = configs
                            realLastConfigGen["all"] = gens[-1]
                        if (configs[:(rp.SEARCH_POPULATION // 10 + 1)]
                                != realLastConfig["dec"]):
                            realLastConfig["dec"] = configs[:(
                                            rp.SEARCH_POPULATION // 10 + 1)]
                            realLastConfigGen["dec"] = gens[-1]
                        if configs[0] != realLastConfig["best"]:
                            realLastConfig["best"] = configs[0]
                            realLastConfigGen["best"] = gens[-1]
                        for k in ["dec", "best", "all"]:
                            if (rp.SEARCH_MAX_DGEN[k] > 0
                                    and len(gens) > 1 and not stop
                                    and rp.GAUSSIAN_WIDTH_SCALING != 1
                                    and (dgen[k] >= rp.SEARCH_MAX_DGEN[k])):
                                stop = True
                                o = {"all": "all structures",
                                     "best": "best structure",
                                     "dec": "best 10% of structures"}
                                logger.info(
                                    "Search convergence criterion reached: "
                                    "max. generations without change ({}): "
                                    "{}/{}.".format(
                                        o[k], dgen[k],
                                        int(rp.SEARCH_MAX_DGEN[k])))
                                break
                    if len(newData) > 0:
                        lastconfig = newData[-1][2]
                    if (check_datafiles and not (stop and repeat)
                            and rp.PARABOLA_FIT["type"] != "none"):
                        check_datafiles = False
                        datafiles = [f for f in os.listdir()
                                     if re.match(r'data\d+\.chem$', f.lower())]
                        if rfac_predict and (rfac_predict[-1][1] >
                                             rfaclist[-1][0]):
                            # if current prediction is bad, reset x0
                            parab_x0 = None
                        try:
                            parab_x0, predictR = parabolaFit(
                                rp, datafiles, min(rfacs), x0=parab_x0,
                                max_configs=max_read_configs)
                            if predictR is not None:
                                if not rfac_predict:
                                    logger.debug("Starting parabola fits "
                                                 "to R-factor data.")
                                rfac_predict.append((gens[-1], predictR))
                        except KeyboardInterrupt:
                            raise
                        except Exception:
                            logger.warning("Parabolic fit of R-factor "
                                           "data failed",
                                           exc_info=rp.LOG_DEBUG)
                    if len(gens) > 1:
                        try:
                            writeSearchProgressPdf(
                                rp, gens, rfaclist, lastconfig,
                                markers=markers, rfac_predict=rfac_predict)
                        except KeyboardInterrupt:
                            raise
                        except Exception:
                            logger.warning("Error writing Search-progress.pdf",
                                           exc_info=rp.LOG_DEBUG)
                        try:
                            writeSearchReportPdf(rp)
                        except KeyboardInterrupt:
                            raise
                        except Exception:
                            logger.warning("Error writing Search-report.pdf",
                                           exc_info=rp.LOG_DEBUG)
                    if (len(gens) > 1 and os.path.isfile("SD.TL")
                            and (repeat or checkrepeat or not stop)):
                        try:
                            processSearchResults(sl, rp, final=False)
                        except Exception as e:
                            logger.warning("Failed to update POSCAR_OUT "
                                           "and VIBROCC_OUT: " + str(e))
                if stop:
                    logger.info("Stopping search...")
                    kill_process(proc)
                    if (not repeat and not rp.GAUSSIAN_WIDTH_SCALING == 1
                            and checkrepeat):
                        repeat = True
                        block = False
                        for k in [k for k in ["dec", "best", "all"]
                                  if rp.SEARCH_MAX_DGEN[k] > 0]:
                            if convergedConfig[k] != realLastConfig[k]:
                                convergedConfig[k] = realLastConfig[k][:]
                                block = True
                            else:
                                repeat = False
                                o = {"all": "any structure",
                                     "best": "best structure",
                                     "dec": "best 10% of structures"}
                                logger.info(
                                    "Convergence reached: No improvement to "
                                    + o[k] + " since changing GAUSSIAN_WIDTH.")
                        if not repeat and block:
                            repeat = True
                            logger.info("Multiple convergence "
                                        "criteria are defined, not all are "
                                        "met. Search continues.")
                        if repeat:
                            if rp.GAUSSIAN_WIDTH == 0.0001:
                                logger.info(
                                    "GAUSSIAN_WIDTH cannot be reduced "
                                    "further, continuing search...")
                            else:
                                rp.GAUSSIAN_WIDTH *= rp.GAUSSIAN_WIDTH_SCALING
                                if rp.GAUSSIAN_WIDTH < 0.0001:
                                    rp.GAUSSIAN_WIDTH = 0.0001
                                logger.info(
                                    "Reducing GAUSSIAN_WIDTH parameter to {} "
                                    "and restarting search..."
                                    .format(round(rp.GAUSSIAN_WIDTH, 4)))
                                comment = ("GAUSSIAN_WIDTH = {}".format(
                                    round(rp.GAUSSIAN_WIDTH, 4)))
                            for k in ["dec", "best", "all"]:
                                rp.SEARCH_MAX_DGEN[k] *= (
                                            rp.SEARCH_MAX_DGEN_SCALING[k])
                                realLastConfigGen[k] = gens[-1]
        except KeyboardInterrupt:
            if not os.path.isfile("SD.TL"):
                # try saving by waiting for SD.TL to be created...
                logger.warning("SD.TL file not found. Trying to wait, "
                               "interrupt again to override...")
                try:
                    i = 0
                    while not os.path.isfile("SD.TL") and not i >= 60:
                        time.sleep(1)
                        i += 1
                except KeyboardInterrupt:
                    pass   # user insisted, give up
            interrupted = True
            rp.STOP = True
            kill_process(proc)
            logger.warning("Search interrupted by user. Attempting "
                           "analysis of results...")
        except Exception:
            logger.error("Error during search. Check SD.TL file.")
            kill_process(proc)
            raise
        if repeat:
            rp.SEARCH_START = "control"
            if gens:
                genOffset = gens[-1]
                rp.SEARCH_MAX_GEN -= sdtlGenNum
                markers.append((genOffset, comment))
            try:
                io.generateSearchInput(sl, rp, steuOnly=True,
                                       cull=True, info=False)
            except Exception:
                logger.error("Error re-generating search input")
                raise
            if os.path.isfile("SD.TL"):
                try:
                    os.remove("SD.TL")
                except Exception:
                    logger.warning("Failed to delete old SD.TL file. "
                                   "This may cause errors in the "
                                   "interpretation of search progress.")
    if proc is not None:
        try:     # should generally not be necessary, but just to make sure
            kill_process(proc, default_pgid=pgid)
        except Exception:
            pass
    if not interrupted:
        logger.info("Finished search. Processing files...")
    else:
        logger.info("Processing files...")
    # final parabola fit  # !!! double check later
    if rfac_predict and (rfac_predict[-1][1] > rfaclist[-1][0]):
        parab_x0 = None
    datafiles = [f for f in os.listdir()
                 if re.match(r'data\d+\.chem$', f.lower())]
    if len(rfaclist) > 0:
        parab_x0, predictR = parabolaFit(rp, datafiles, min(rfacs),
                                         x0=parab_x0,
                                         max_configs=max_read_configs)
        if predictR is not None:
            rfac_predict.append((gens[-1], predictR))
    # write pdf one more time
    if len(gens) > 1:
        try:
            writeSearchProgressPdf(rp, gens, rfaclist, lastconfig,
                                   markers=markers, rfac_predict=rfac_predict)
        except Exception:
            logger.warning("Error writing Search-progress.pdf",
                           exc_info=True)
        try:
            writeSearchReportPdf(rp)
        except Exception:
            logger.warning("Error writing Search-report.pdf",
                           exc_info=True)
    # process SD.TL to get POSCAR_OUT, VIBROCC_OUT
    try:
        processSearchResults(sl, rp)
    except FileNotFoundError:
        logger.error("Cannot interpret search results without SD.TL file.")
        rp.setHaltingLevel(2)
    except Exception:
        logger.error("Error processing search results: ", exc_info=True)
        raise
    # if deltas were copied from domain folders, clean them up
    if rp.domainParams:
        rgx = re.compile(r'D\d+_DEL_')
        for file in [f for f in os.listdir() if rgx.match(f)]:
            try:
                os.remove(file)
            except Exception:
                logger.warning('Failed to deleted redundant domain delta file '
                               + file)
    # process files
    try:
        os.rename('PARAM', 'search-PARAM')
    except Exception:
        logger.warning("Failed to rename search input file PARAM to "
                       "search-PARAM")
    try:
        os.rename('rf.info', 'search-rf.info')
    except Exception:
        logger.warning("Failed to rename search input file rf.info to "
                       "search-rf.info")
    if lastconfig is not None:
        rp.searchResultConfig = lastconfig
    return
