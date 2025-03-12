"""Section Search."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-11'
__license__ = 'GPLv3+'

from collections import Counter
import copy
import logging
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time

import numpy as np
import scipy
from sklearn.linear_model import ElasticNet
from sklearn.linear_model import Lasso
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures

from viperleed.calc.classes.searchpar import SearchPar
from viperleed.calc.files import iosearch
from viperleed.calc.files import parameters
from viperleed.calc.files import searchpdf
from viperleed.calc.files.displacements import readDISPLACEMENTS_block
from viperleed.calc.lib import leedbase
from viperleed.calc.lib.checksums import validate_multiple_files
from viperleed.calc.lib.time_utils import ExecutionTimer
from viperleed.calc.lib.time_utils import ExpiringOnCountTimer
from viperleed.calc.lib.time_utils import ExpiringTimerWithDeadline
from viperleed.calc.lib.version import Version

logger = logging.getLogger(__name__)


class SearchError(Exception):
    """Base class for exceptions in this module."""

    base_message = ''
    detailed_message = ''
    log_records = ()   # Strings in the log that signal this exception
    _matcher = all     # Should all or any of the log records match?

    def __init__(self, *args, **kwargs):
        """Initialize exception instance."""
        message = (self.base_message
                   + (' ' if self.base_message else '')
                   + self.detailed_message)
        if message:
            args = message, *args
        super().__init__(*args, **kwargs)

    @classmethod
    def matches(cls, log_contents):
        """Return whether this exception matches `log_contents`."""
        return cls._matcher(rec in log_contents for rec in cls.log_records)


class SearchParameterError(SearchError):
    """A value in PARAMETERS is inappropriate for the search."""


class TensErLEEDSearchError(SearchError):
    """A known, fatal TensErLEED error occurred during the search."""

    base_message = ('TensErLEED error encountered in '
                    'search. Execution cannot proceed.')


class InconsistentV0ImagError(SearchParameterError, TensErLEEDSearchError):
    """The V0i value in PARAMETERS does not match the one in the Deltas."""

    detailed_message = (
        'TensErLEED search stopped because stored Delta files were calculated '
        'with a different imaginary inner potential (optical potential) than '
        'requested in PARAMETERS. Make sure to use the same value of V0_IMAG '
        'for delta amplitudes and structure search.'
        )
    log_records = ('Average optical potential value in rf.info is incorrect:',)


class MaxIntensitiesError(SearchParameterError, TensErLEEDSearchError):
    """Search failed with "too high intensity" error."""

    detailed_message = (
        'TensErLEED stopped due to unreasonably high beam amplitudes.\n'
        'This may be caused by convergence problems due to a too small '
        'value of LMAX. Alternatively, this error may be caused by '
        'scatterers with very small distances as a result of a very '
        'large range used in DISPLACEMENTS. Check your input files and '
        'consider increasing LMAX or decreasing the DISPLACEMENTS ranges.'
        )
    log_records = ('MAX. INTENS. IN THEOR. BEAM',
                   'IS SUSPECT',
                   '****** STOP PROGRAM ******')


class NotEnoughSlotsError(SearchParameterError):
    """Could not spawn the requested N_CORES Open MPI processes."""

    detailed_message = (
        'GNU mpirun failed to allocate the number of processes '
        'requested via N_CORES. Try reducing the N_CORES parameter.'
        )
    log_records = ('There are not enough slots available in the system',)


class ProcessKilledError(SearchError):
    """The mpirun process died. Typically, it required too much memory."""

    detailed_message = (
        'The search was stopped because the MPI process was killed by the '
        'system. This is most likely because the process ran out of available '
        'memory. Consider reducing the number of MPI processes by setting a '
        'smaller value for N_CORES or limiting the parameter space in the '
        'DISPLACEMENTS file.'
        )
    _matcher = any  # Two different error messages for Intel and GNU
    log_records = (
        'APPLICATION TERMINATED WITH THE EXIT STRING: Killed (signal 9)',
        '=   KILLED BY SIGNAL: 9 (Killed)',
        )


class SigbusError(TensErLEEDSearchError):
    """Tried to write to a non-existing memory address."""

    detailed_message = (
        'TensErLEED stopped due to a SIGBUS signal.\nThis may be caused by a '
        'compiler error. If you are using gfortran, consider lowering the '
        'optimization level to \'-O1\' or \'-O0\' using the FORTRAN_COMP '
        'parameter.'
        )
    log_records = ('received signal SIGBUS',
                   'undefined portion of a memory object')


def processSearchResults(sl, rp, search_log_path, final=True):
    """Read the best structure from the last block of 'SD.TL' into a slab.

    Parameters
    ----------
    sl : Slab
        The Slab object to be modified
    rp : Rparams
        Run parameters.
    search_log_path : Pathlike or None
        Path to the search log. Will be checked for TensErLEED errors
        messages. If None, assumes search was not logged and skips
        checks.
    final : bool, optional
        Defines whether this is the final interpretation that *must* work and
        should go into the log, or if this will be repeated anyway. The
        default is True.

    Returns
    -------
    None.
    """
    # Read search log if available, and check for errors
    _check_search_log(search_log_path)

    # get the last block from SD.TL:
    wait_time = 5 if not final else 20
    try:
        sdtl_content = iosearch.repeat_fetch_SDTL_last_block(
            which_beams=rp.SEARCH_BEAMS,
            expected_params=rp.SEARCH_POPULATION,
            final=final, wait_time=wait_time
            )
    except iosearch.SearchIOEmptyFileError as err:
        if final:
            logger.error("No data found in SD.TL file!")
            rp.setHaltingLevel(2)
        raise RuntimeError("No data in SD.TL") from err
    except iosearch.SearchIORaceConditionError:
        return  # if we can't find a block, we also can't store it.

    (generation, rfactors, configs), *_ = sdtl_content

    assert len(rfactors) == len(configs)                                        # TODO: catch and complain. Previous version would pick only the first len(rfactors) configs, but the two should hopefully be the same length.
    writeControlChem = False                                                    # TODO: I don't think we should do this at all. We might mess with the fortran subprocess that is trying to write to the file.
    if not os.path.isfile("control.chem"):
        writeControlChem = True
    _write_control_chem(rp, generation, configs)

    # Collect populations by (r-factor, parameter indices) pairs
    populations = Counter(zip(rfactors, configs))
    if final:
        _write_info_to_log(rp, populations)


    # Pick only the best configuration (possibly for multiple
    # domains): it is the one with the highest R factor, i.e.,
    # the first one. best_doms is a list of pairs (domain_area,
    # domain_parameters) for each domain under variation
    best_doms = configs[0]
    if rp.domainParams:
        doms_info = ((dp.sl, dp.rp, dp.workdir, dp.name)
                     for dp in rp.domainParams)
    else:
        doms_info = ((sl, rp, Path.cwd(), ""),)

    # Finally write out the best structures
    home = Path.cwd()
    for (*dom_info, work, name), (_, dom_config) in zip(doms_info, best_doms):
        os.chdir(work)
        try:
            _store_and_write_best_structure(rp, *dom_info,
                                            dom_config, final)
        except Exception:                                                       # TODO: catch better
            _err = "Error while writing search output"
            if name:
                _err += f" for domain {name}"
            logger.error(_err, exc_info=rp.is_debug_mode)
            rp.setHaltingLevel(2)
            raise
        finally:
            os.chdir(home)


def _write_control_chem(rp, generation, configs):
    """Store information in configs in rp, and write it to file if needed."""
    _ctrl_chem_p = Path("control.chem")
    _write_file = not _ctrl_chem_p.is_file()
    if not _write_file:
        # check file, which generation
        try:
            with _ctrl_chem_p.open("r", encoding="utf-8") as control_chem:
                _, header, *other_lines = control_chem
        except (OSError, ValueError):
            # try overwriting, just to be safe
            _write_file = True
        else:
            try:
                gen_nr = int(header.split("No.")[-1].split(":")[0])
            except (IndexError, ValueError):
                _write_file = True
            else:
                _write_file = (gen_nr != generation
                               or len(other_lines) < rp.SEARCH_POPULATION)

    # Prepare lines to store as backup, and perhaps write to file
    lines = ["\n",  # Starts with one empty line
             f"Parameters of generation No.{generation:>6d}:\n"]
    astep = rp.DOMAIN_STEP if rp.domainParams else 100
    if rp.TL_VERSION < Version('1.7.0'):
        ctrl_width = 3
    else:
        ctrl_width = 4
    for dom_pars in configs:
        areapars = (int(percent/astep) + 1 for percent, _ in dom_pars)
        lines.append(
            "".join(f"{i:>{ctrl_width}d}" for _, par in dom_pars for i in par)
            + "".join(f"{a:>{ctrl_width}d}" for a in areapars) + "\n"
            )
    rp.controlChemBackup = "".join(lines)
    if not _write_file:
        return
    try:
        with _ctrl_chem_p.open("w", encoding="utf-8") as control_chem:
            control_chem.writelines(lines)                                      # TODO: I don't think we should do this at all. We might mess with the fortran subprocess that is trying to write to the file.
    except Exception:
        logger.error("Failed to write control.chem")
        rp.setHaltingLevel(1)


def _write_info_to_log(rp, pops):
    """Write information about populations to the current logger."""
    (most_common_r, _), n_most_common = pops.most_common(1)[0]
    if n_most_common == rp.SEARCH_POPULATION:
        logger.info("All trial structures converged to the same "
                    f"configuration (R = {most_common_r:.4f})")
        return

    (best_r, _), n_best = next(iter(pops.items()))
    # if any(v > rp.SEARCH_POPULATION / 2 for v in popcount):                   # TODO: used to be like this, but I think the current one is the same.
    if n_most_common > rp.SEARCH_POPULATION / 2:
        info = ("The search outcome is dominated by one configuration "
                f"(population {n_most_common} of {rp.SEARCH_POPULATION}, "
                f"R = {most_common_r:.4f})\n")
    else:
        info = ("The search outcome is not dominated by any "
                "configuration. The best result has population "
                f"{n_best} (of {rp.SEARCH_POPULATION}), "
                f"R = {best_r:.4f}\n")
    # if len(pops) > 1:                                                         # TODO: what follows used to be done only in this case, but I think it should already be covered by the first "if n_most_common == rp.SEARCH_POPULATION"
    info += ("The best configurations are:\nPOP       R |")
    if rp.domainParams:
        info += " area |"
    info += " PARAMETERS\n"
    for (pop_r, pop_pars), n_pops in pops.most_common(5):
        info += f"{n_pops:>3}  {pop_r:.4f} |"
        for j, (percent, pars) in enumerate(pop_pars):
            if j:
                info += "            |"
            if rp.domainParams:
                info += f" {percent:>3}% |"
            info += "".join(f"{v:>3}" for v in pars) + "\n"
    logger.info(info)


def _store_and_write_best_structure(rp, dom_slab, dom_rp, best_config, final):
    """Modify dom_slab to the parameters in best_config. Write it out.

    If there is a predicted optimum for the parameters (via parabola
    fit), POSCAR_OUT and VIBROCC_OUT files are also written for this
    predicted optimum. However, `dom_slab` will always contain the
    configuration passed via `best_config` at the end of this call.

    Parameters
    ----------
    rp : Rparams
        The PARAMETERS for the whole calculation.
    dom_slab : Slab
        The slab of the domain to be modified, and used to write
        POSCAR_OUT and VIBROCC_OUT files.
    dom_rp : Rparams
        The PARAMETERS for this specific domain to be written out.
    best_config : Sequence
        The indices of the optimal parameters to be used for this
        domain's configuration.

    Returns
    -------
    None.
    """
    # The order here matters, as the final write operation
    # will store the new states (e.g., for later searches)
    # in dom_slab, its atoms, and its sites.
    if any(sp.parabolaFit["min"] for sp in rp.searchpars):
        parab_inds = list(best_config)
        for j, sp in enumerate(dom_rp.searchpars):
            if sp.parabolaFit["min"] is not None:
                parab_inds[j] = sp.parabolaFit["min"]
        iosearch.writeSearchOutput(dom_slab, dom_rp, parab_inds,
                                   silent=True, suffix="_parabola")
    iosearch.writeSearchOutput(dom_slab, dom_rp, best_config, silent=not final)


def _check_search_log(search_log_path):
    """Check the file at `search_log_path` for known fatal errors.

    Parameters
    ----------
    search_log_path : str or Path or None
        Path to the log file. No checking is performed if None.

    Raises
    ------
    SearchError
        If error messages are found in the log file at `search_log_path`
        that correspond to a fatal error (from TensErLEED or other
        source).
    """
    if not search_log_path:
        return

    search_log_path = Path(search_log_path)
    try:
        log_contents = search_log_path.read_text(encoding='utf-8')
    except OSError:
        logger.error(f'Could not read search log file {search_log_path}. '
                     'Execution will continue, but errors related '
                     'to the search may not be reported properly.')
        return
    _known_errors = (
        MaxIntensitiesError,
        SigbusError,
        InconsistentV0ImagError,
        ProcessKilledError,
        NotEnoughSlotsError,
        )
    try:
        raise next(e for e in _known_errors if e.matches(log_contents))
    except StopIteration:  # No error
        pass


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
        varR = np.sqrt(8*abs(rp.V0_IMAG) / rp.total_energy_range())*r_best
        rc = np.array((iosearch.readDataChem(
            rp, datafiles,
            cutoff=r_best + r_cutoff * varR,
            max_configs=max_configs
            )), dtype=object)
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
        indep_pars = np.array([*np.array([*np.array(configs, dtype=object)],    # TODO: pylint complains about too-many-function-args. I don't understand the line
                                         dtype=object)
                               .reshape(-1, 1, 2)[:, 0, 1]])
        indep_pars = np.delete(indep_pars, [i for i in
                                            range(len(rp.searchpars)-1)
                                            if rp.searchpars[i] not in sps], 1)
    else:
        # first the percentages:
        sps = [sp for sp in rp.searchpars if sp.mode == "dom"]
        reshaped = (np.array([*np.array(configs, dtype=object)], dtype=object)  # TODO: pylint complains about too-many-function-args. I don't understand the line
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
                                Lasso(alpha=alpha, normalize=True))             # TODO: incorrect keyword normalize
    elif which_regression == 'ridge':
        polyreg = make_pipeline(PolynomialFeatures(degree=2),
                                Ridge(alpha=alpha, normalize=True))             # TODO: incorrect keyword normalize
    elif which_regression == 'elasticnet':
        polyreg = make_pipeline(PolynomialFeatures(degree=2),
                                ElasticNet(alpha=alpha, normalize=True))        # TODO: incorrect keyword normalize
    else:
        if which_regression not in ('linearregression', 'linear'):
            logger.warning("Regression model {} not found, parabola fit "
                           "defaulting to linear regression."
                           .format(which_regression))
        which_regression = "linearregression"
        polyreg = make_pipeline(PolynomialFeatures(degree=2),
                                LinearRegression())
    xmin = np.copy(indep_pars[np.argmin(rfacs), :])
    rr = np.sqrt(8*abs(rp.V0_IMAG) / rp.total_energy_range())
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
    #         dist = abs(dist_norm*(ip_tofit - base))
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
    for sp in [sp for sp in rp.searchpars if isinstance(sp.linkedTo, SearchPar)
               or isinstance(sp.restrictTo, SearchPar)]:
        if isinstance(sp.linkedTo, SearchPar):
            sp.parabolaFit = copy.copy(sp.linkedTo.parabolaFit)
        elif isinstance(sp.restrictTo, SearchPar):
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
            pgid = os.getpgid(proc.pid)                                         # TODO: only UNIX!
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
                os.killpg(pgid, signal.SIGTERM)                                 # TODO: only UNIX
                if os.name == "nt":
                    os.waitpid(pgid)
                else:
                    os.waitpid(-pgid, 0)
            except (ProcessLookupError, ChildProcessError):
                pass  # already dead or no children

    rp.searchResultConfig = None
    if rp.domainParams:
        initToDo = [(dp.rp, dp.sl, dp.workdir) for dp in rp.domainParams]
    else:
        initToDo = [(rp, sl, rp.paths.work)]
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
            leedbase.getDeltas(rpt.TENSOR_INDEX, basedir=path,
                               targetdir=path, required=True)
    rp.updateCores()
    # generate rf.info
    try:
        rf_info_path = rp.paths.work / "rf.info"
        rf_info_content = iosearch.writeRfInfo(sl, rp, file_path=rf_info_path)
    except Exception:
        logger.error("Error generating search input file rf.info")
        raise
    # generate PARAM and search.steu
    #   needs to go AFTER rf.info, as writeRfInfo may remove expbeams!
    try:
        iosearch.generateSearchInput(sl, rp)
    except Exception:
        logger.error("Error generating search input")
        raise
    # delete old control.chem if present - not earlier because can be input!
    try:
        os.remove("control.chem")
    except FileNotFoundError:
        pass
    if not rp.indyPars:  # never for calculations with domains                  # !!! CHECK
        logger.info("Found nothing to vary in search. Will proceed "
                    "directly to writing output and starting SUPERPOS.")
        this_domain_pars = [1] * len(rp.searchpars)                             # TODO: Used to be "[1] * (len) - 1" --> TypeError
        for (i, sp) in enumerate(rp.searchpars):
            if isinstance(sp.restrictTo, int):
                this_domain_pars[i] = sp.restrictTo
            elif isinstance(sp.restrictTo, SearchPar):
                this_domain_pars[i] = rp.searchpars.index(sp.restrictTo) + 1
            elif isinstance(sp.linkedTo, SearchPar):
                this_domain_pars[i] = rp.searchpars.index(sp.linkedTo) + 1
        rp.searchResultConfig = ((100, this_domain_pars),)                      # TODO: Used to assign straight to rp.searchResultConfig[i] -> IndexError as there was only a single element.
        iosearch.writeSearchOutput(sl, rp)
        return None

    if rp.SUPPRESS_EXECUTION:
        logger.warning("SUPPRESS_EXECUTION parameter is on. Search "
                       "will not proceed. Stopping...")
        rp.setHaltingLevel(3)
        return None

    # Decide whether to use parallelization
    usempi = rp.N_CORES > 1 and shutil.which('mpirun')
    if not usempi:
        reason = (
            f'The N_CORES parameter is set to {rp.N_CORES}' if rp.N_CORES <= 1
            else 'mpirun is not present'
            )
        logger.warning(f'{reason}. The search will be run without '
                       'parallelization. This will be much slower!')

    _find_compiler = None
    if usempi and not rp.FORTRAN_COMP_MPI[0]:
        _find_compiler = rp.getFortranMpiComp
    elif not rp.FORTRAN_COMP[0]:
        _find_compiler = rp.getFortranComp
    if _find_compiler:
        try:
            _find_compiler()
        except Exception as exc:
            _mpi = 'mpi ' if usempi else ''
            logger.error(f'No fortran {_mpi}compiler found, cancelling...')
            raise FileNotFoundError('Fortran compile error') from exc

    # get fortran files
    try:
        tl_source = rp.get_tenserleed_directory()
        tl_path = tl_source.path
        srcpath = tl_path / 'src'
        if usempi:
            src_file = next(srcpath.glob('search.mpi*'), None)
        else:
            src_files = (f for f in srcpath.glob('search*')
                         if 'mpi' not in f.name)
            src_file = next(src_files, None)
        if src_file is None:
            raise FileNotFoundError('No Fortran source for '
                                    f'search in {tl_path}')
        shutil.copy2(src_file, src_file.name)
        libpath = tl_path / 'lib'
        libpattern = 'lib.search'
        if usempi and rp.TL_VERSION <= Version('1.7.3'):
            libpattern += '.mpi'
        lib_file = next(libpath.glob(libpattern + '*'), None)
        if lib_file is None:
            raise FileNotFoundError(f'File {libpattern}.f not found.')

        # copy to work dir
        shutil.copy2(lib_file, lib_file.name)
        hashing_file = next(libpath.glob("intarr_hashing*.f90"), None)
        if hashing_file:
            hashname = hashing_file.name
            shutil.copy2(hashing_file, hashname)
        else:
            hashname = ""

        randname = 'MPIrandom_.o' if usempi else 'random_.o'
        if rp.TL_VERSION <= Version('1.7.3'):
            # Try to copy randomizer lib object file
            # These are short C scripts - use pre-compiled versions
            try:
                shutil.copy2(libpath / randname, randname)
            except FileNotFoundError:
                logger.error(f'Could not find required {randname} file. '
                             'You may have forgotten to compile random_.c.')
                raise

        globalname = "GLOBAL"
        shutil.copy2(srcpath / globalname, globalname)
    except Exception:
        logger.error("Error getting TensErLEED files for search: ")
        raise

    # Validate TensErLEED input files
    if not rp.TL_IGNORE_CHECKSUM:
        files_to_check = (lib_file,
                          src_file,
                          srcpath / globalname,
                          hashing_file)
        validate_multiple_files(files_to_check, logger,
                                "search", rp.TL_VERSION)

    # Prepare to compile fortran files
    searchname = f'search-{rp.timestamp}'
    fcomp = rp.FORTRAN_COMP_MPI if usempi else rp.FORTRAN_COMP
    logger.info('Compiling fortran input files...')

    # compile task could be inherited from general CompileTask (issue #43)
    ctasks = [(f"{fcomp[0]} -o lib.search.o -c", lib_file.name, fcomp[1])]
    if hashname:
        ctasks.append((f"{fcomp[0]} -c", hashname, fcomp[1]))
    ctasks.append((f"{fcomp[0]} -o restrict.o -c", "restrict.f", fcomp[1]))
    _fixed_format = any(f.endswith('.f90')
                        for f in (lib_file.name, src_file.name, hashname))
    is_gfortran = any(s in fcomp[0]
                      # Assume that mpifort uses gfortran
                      for s in ('gfortran', 'mpifort'))
    format_tag = ''
    if _fixed_format:
        format_tag =  '--fixed-form' if is_gfortran else '-fixed'
    ctasks.append(
        (f"{fcomp[0]} -o search.o -c {format_tag}", src_file.name, fcomp[1])
        )
    to_link = "search.o lib.search.o restrict.o"
    if rp.TL_VERSION <= Version('1.7.3'):
        to_link += f" {randname}"
    if hashname:
        to_link += " intarr_hashing.o"
    ctasks.append((f"{fcomp[0]} -o {searchname}", to_link, fcomp[1]))
    compile_log = "compile-search.log"
    try:
        leedbase.fortran_compile_batch(ctasks, logname=compile_log)
    except Exception:
        leedbase.copy_compile_log(rp, Path(compile_log),
                                  log_name="search-compile")
        logger.error("Error compiling fortran files: ", exc_info=True)
        raise
    logger.debug("Compiled fortran files successfully")
    # run
    if rp.LOG_SEARCH:
        search_log_path = (rp.paths.work / searchname).with_suffix(".log")
        logger.info(f"Search log will be written to file {search_log_path}.")
    else:
        search_log_path = None
    # if there is an old SD.TL file, it needs to be removed
    if os.path.isfile("SD.TL"):
        try:
            os.remove("SD.TL")
        except OSError:
            logger.warning("Failed to delete old SD.TL file. This may "
                           "cause errors in the interpretation of search "
                           "progress.")
    # same for old data.chem
    _datachem_re = re.compile(r'data\d+\.chem$', re.IGNORECASE)
    for file in Path().glob("*"):
        if not _datachem_re.match(file.name):
            continue
        try:
            file.unlink()
        except OSError:
            logger.warning(f"Failed to delete old {file} file. This "
                           "may cause errors in the parabola fit.")
    # get config size
    config_size = (
        sys.getsizeof((1., tuple()))
        + sys.getsizeof((tuple(),) * max(1, len(rp.domainParams)))
        + sys.getsizeof((1, tuple()))
        )
    if rp.domainParams:
        config_size += sum(
            sys.getsizeof((0,) * len(dp.rp.searchpars))
            + sys.getsizeof(1) * len(dp.rp.searchpars)
            for dp in rp.domainParams
            )
    else:
        config_size += (sys.getsizeof((0,) * len(rp.searchpars))
                        + sys.getsizeof(1) * len(rp.searchpars))
    max_memory_parab = 1e9   # maximum memory in bytes to use on parabola fit
    max_read_configs = int(max_memory_parab / config_size)
    # start search process
    repeat = True
    genOffset = 0
    gens = []  # generation numbers in SD.TL, but continuous if search restarts
    timestamps = []  # timestamps of generations
    markers = []
    rfaclist = []
    parab_x0 = None    # starting guess for parabola                            # TODO: would be nicer to incorporate it into a ParabolaFit class
    rfac_predict = []  # tuples (gen, r) from parabola fit                      # TODO: would be nicer to incorporate it into a ParabolaFit class
    realLastConfig = {"all": [], "best": [], "dec": []}                         # TODO: would be nicer with a SearchConfigurationTracker class that has an .update
    realLastConfigGen = {"all": 0, "best": 0, "dec": 0}
    convergedConfig = {"all": None, "best": None, "dec": None}
    lastconfig = None
    rp.searchMaxGenInit = rp.SEARCH_MAX_GEN
    since_started = ExecutionTimer()
    since_last_debug = ExpiringOnCountTimer(
            interval=rp.output_interval,
            count_start=0,
            )
    tried_repeat = False  # if SD.TL is not written, try restarting
    pgid = None
    logger.info("Starting search. See files Search-progress.pdf "
                "and SD.TL for progress information.")

    # Prepare the command to be run via subprocess
    executable = os.path.join('.', searchname)
    if usempi:
        command = ['mpirun', '-n', str(rp.N_CORES)]
    if usempi and is_gfortran:
        # Assume we're using OpenMPI: we need to specify the use of all
        # CPU threads explicitly, otherwise OpenMPI will use only the
        # physical cores. However, our auto-detected N_CORES counts all
        # the logical cores, not only the physical ones.
        command.append('--use-hwthread-cpus')
    command.append(executable)

    while repeat:                                                               # TODO: all this mess would be nicer to handle with a state machine approach. This would at least help readability on the various ways things are handled
        repeat = False
        interrupted = False
        proc = None
        # if LOG_SEARCH -> log search
        if search_log_path:
            log_exists = search_log_path.is_file()
            search_log_f = search_log_path.open("a")
            if log_exists:  # log file existed before
                search_log_f.write("\n\n-------\nRESTARTING\n-------\n\n")
        else:
            search_log_f = subprocess.DEVNULL
        # NB: log file may be open and must be closed!
        # Create search process
        logger.debug(f'Starting search process with command "{" ".join(command)}".')
        try:
            proc = subprocess.Popen(
                command, encoding="ascii",
                stdout=search_log_f,  # if LOG_SEARCH is False, this is DEVNULL
                stderr=search_log_f,  # same as above
                preexec_fn=os.setsid                                            # TODO: setsid only POSIX; os comments suggest NOT TO USE THIS as it is unsafe against deadlocks. Suggestion is to use start_new_session keyword [UNIX only!] instead of setsid. For a Windows solution see https://stackoverflow.com/questions/47016723. Probably even better: do not use python subprocess, but QProcess (which can run with its own event loop, and has a .kill(), or perhaps .terminate()). QProcess may have some issues with OpenMPI v<=1.7 due to a bug there (see www.qtcentre.org/threads/19636-Qprocess-and-mpi-not-finishing).
                )
        except OSError:  # This should not fail unless the shell is very broken.
            logger.error("Error starting search. Check SD.TL file.")
            if search_log_path:
                search_log_f.close()
            raise
        else:
            pgid = os.getpgid(proc.pid)                                         # TODO: getpgid only POSIX
        if proc is None:
            logger.error("Error starting search subprocess... Stopping.")
            if search_log_path:
                search_log_f.close()
            raise RuntimeError("Error running search."
                               f'Could not start process using "{command}".')
        # FEED INPUT
        try:
            proc.communicate(input=rf_info_content, timeout=0.2)
        except subprocess.TimeoutExpired:
            pass  # started successfully; monitoring below
        except (OSError, subprocess.SubprocessError):
            logger.error("Error feeding input to search process. "
                         "Check files SD.TL and rf.info.")

        # MONITOR SEARCH
        search_eval_timer = ExpiringTimerWithDeadline(                          # TODO: would be nicer with a QTimer, or even a QFileSystemWatcher
            interval=rp.searchEvalTime,  # How often to read SD.TL
            deadline=900,   # Max 15 min without an SD.TL produced
            )
        since_last_debug.synchronize_with(search_eval_timer)
        filepos = 0
        timestep = 1  # time step to check files
        # !!! evaluation time could be higher - keep low only for debugging; TODO
        comment = ""
        sdtlGenNum = 0
        gaussianWidthOri = rp.GAUSSIAN_WIDTH
        check_datafiles = False
        try:
            while proc.poll() is None:  # proc is running
                time.sleep(timestep)
                parameters.update(rp)                                           # TODO: Would be way nicer with a QFileSystemWatcher
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
                        while not os.path.isfile("SD.TL") and i < 300:
                            time.sleep(1)
                            i += 1
                if rp.GAUSSIAN_WIDTH != gaussianWidthOri:
                    stop = True
                    repeat = True
                    comment = f"GAUSSIAN_WIDTH = {rp.GAUSSIAN_WIDTH}"
                    logger.info("GAUSSIAN_WIDTH parameter changed. "
                                "Search will restart.")
                if search_eval_timer.has_expired() or stop:
                    # evaluate
                    newData = []
                    if os.path.isfile("SD.TL"):
                        filepos, content = iosearch.readSDTL_next(
                            offset=filepos
                            )
                        if content:
                            newData = iosearch.readSDTL_blocks(
                                content, whichR=rp.SEARCH_BEAMS,
                                n_expect=rp.SEARCH_POPULATION
                                )
                    elif (search_eval_timer.has_reached_deadline()              # TODO: nicer with a QTimer timeout
                          and rp.HALTING < 3):
                        stop = True
                        if tried_repeat:
                            logger.warning(
                                "No SD.TL file was written for 15 minutes "
                                "after restarting the search. Search will "
                                "stop."
                                )
                            repeat = False
                            rp.setHaltingLevel(2)
                        else:
                            repeat = True
                            tried_repeat = True
                            logger.warning(
                                "No SD.TL file was written for 15 minutes "
                                "after the search started. Trying to restart. "
                                "You can suppress this behaviour by setting "
                                "the HALTING parameter to 3."
                                )
                    for gen, rfacs, configs in newData:
                        gens.append(gen + genOffset)
                        timestamps.append(time.time())
                        sdtlGenNum = gen
                        rfaclist.append(np.array(rfacs))
                        dgen = {}
                        for k in ["dec", "best", "all"]:                        # TODO: would be nicer with a SearchConfigurationTracker class that has an .update
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
                                    and len(gens) > 1
                                    and not stop
                                    and rp.GAUSSIAN_WIDTH_SCALING != 1
                                    and dgen[k] >= rp.SEARCH_MAX_DGEN[k]):
                                stop = True
                                o = {"all": "all structures",
                                     "best": "best structure",
                                     "dec": "best 10% of structures"}
                                logger.info(
                                    "Search convergence criterion reached: "
                                    f"max. generations without change ({o[k]})"
                                    f": {dgen[k]}/{rp.SEARCH_MAX_DGEN[k]}."
                                    )
                                break
                    # decide to write debug info
                    # will only write once per SD.TL read
                    current_gen = gens[-1] if gens else 0
                    if since_last_debug.count_expired(current_gen):
                        # "speed" is actually the inverse of a
                        # speed, in seconds per 1000 generations
                        speed = 1000 * since_started.how_long() / current_gen
                        logger.debug(
                            f'R = {min(rfacs)} (Generation {current_gen}, '
                            f'{since_last_debug.how_long():.3f} s since '
                            f'gen. {since_last_debug.previous_count}, '
                            f'{speed:.1f} s/kG overall)'
                            )
                        since_last_debug.restart()
                        check_datafiles = True
                    if newData:
                        _, _, lastconfig = newData[-1]
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
                                rp, datafiles, min(rfacs), x0=parab_x0,         # TODO: is rfacs correct here? may be undefined if the loop before does not run
                                max_configs=max_read_configs
                                )
                            if predictR is not None:
                                if not rfac_predict:
                                    logger.debug("Starting parabola fits "
                                                 "to R-factor data.")
                                rfac_predict.append((gens[-1], predictR))
                        except Exception:
                            logger.warning("Parabolic fit of R-factor "
                                           "data failed",
                                           exc_info=rp.is_debug_mode)
                    if len(gens) > 1:
                        try:
                            searchpdf.writeSearchProgressPdf(
                                rp, gens, rfaclist, lastconfig,
                                markers=markers, rfac_predict=rfac_predict,
                                timestamps=timestamps
                                )
                        except Exception:
                            logger.warning("Error writing Search-progress.pdf",
                                           exc_info=rp.is_debug_mode)
                        try:
                            searchpdf.writeSearchReportPdf(rp)
                        except Exception:
                            logger.warning("Error writing Search-report.pdf",
                                           exc_info=rp.is_debug_mode)
                    if (len(gens) > 1 and os.path.isfile("SD.TL")
                            and (repeat or checkrepeat or not stop)):
                        try:
                            processSearchResults(sl, rp, search_log_path,
                                                 final=False)
                        except Exception as exc:                                # TODO: too general
                            logger.warning("Failed to update POSCAR_OUT "
                                           f"and VIBROCC_OUT: {exc}")
                if stop:
                    logger.info("Stopping search...")
                    kill_process(proc, default_pgid=pgid)
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
                            if rp.GAUSSIAN_WIDTH <= 0.0001:
                                logger.info(
                                    "GAUSSIAN_WIDTH cannot be reduced "
                                    "further, continuing search...")
                            else:
                                rp.GAUSSIAN_WIDTH = max(
                                    rp.GAUSSIAN_WIDTH * rp.GAUSSIAN_WIDTH_SCALING,
                                    0.0001
                                )
                                logger.info(
                                    "Reducing GAUSSIAN_WIDTH parameter to {} "
                                    "and restarting search..."
                                    .format(round(rp.GAUSSIAN_WIDTH, 4)))
                                comment = ("GAUSSIAN_WIDTH = {}".format(
                                    round(rp.GAUSSIAN_WIDTH, 4)))
                            for k in ["dec", "best", "all"]:
                                rp.SEARCH_MAX_DGEN[k] = int(
                                            rp.SEARCH_MAX_DGEN[k] *
                                            rp.SEARCH_MAX_DGEN_SCALING[k])
                                realLastConfigGen[k] = gens[-1] if gens else 0
        except KeyboardInterrupt:                                               # TODO: would probably be nicer to install a custom signal handler for SIGINT, CTRL_C_EVENT, and CTRL_BREAK_EVENT, then place back the previous handler when done processing.
            if not os.path.isfile("SD.TL"):
                # try saving by waiting for SD.TL to be created...
                logger.warning("SD.TL file not found. Trying to wait, "
                               "interrupt again to override...")
                try:
                    i = 0
                    while not os.path.isfile("SD.TL") and i < 60:
                        time.sleep(1)
                        i += 1
                except KeyboardInterrupt:
                    pass   # user insisted, give up
            interrupted = True
            rp.STOP = True
            kill_process(proc, default_pgid=pgid)
            logger.warning("Search interrupted by user. Attempting "
                           "analysis of results...")
        except Exception:
            logger.error("Error during search. Check SD.TL file.")
            kill_process(proc, default_pgid=pgid)
            raise
        finally:
            # close open input and log files
            if search_log_path:
                search_log_f.close()
        if repeat:
            rp.SEARCH_START = "control"
            if gens:
                genOffset = gens[-1]
                rp.SEARCH_MAX_GEN -= sdtlGenNum
                markers.append((genOffset, comment))
            try:
                iosearch.generateSearchInput(sl, rp, steuOnly=True,
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
            searchpdf.writeSearchProgressPdf(rp, gens, rfaclist, lastconfig,
                                             markers=markers,
                                             rfac_predict=rfac_predict,
                                             timestamps=timestamps)
        except Exception:
            logger.warning("Error writing Search-progress.pdf",
                           exc_info=True)
        try:
            searchpdf.writeSearchReportPdf(rp)
        except Exception:
            logger.warning("Error writing Search-report.pdf",
                           exc_info=True)
    # process SD.TL to get POSCAR_OUT, VIBROCC_OUT
    try:
        processSearchResults(sl, rp, search_log_path)
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
