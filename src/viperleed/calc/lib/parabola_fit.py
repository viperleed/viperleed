"""Lib parabola_fit (to be deprecated)."""

__authors__ = (
    "Florian Kraushofer (@fkraushofer)",
    "Alexander M. Imre (@amimre)",
)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-07-16"
__license__ = "GPLv3+"

import copy
import numpy as np
import scipy
from sklearn.linear_model import ElasticNet
from sklearn.linear_model import Lasso
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures


from viperleed.calc.files import iosearch
from viperleed.calc.classes.searchpar import SearchPar

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
        start = dim + 1
        for i in range(0, dim):
            m[i][i:] = features[start : start + dim - i]
            start += dim - i
        m = 0.5 * (m + m.T)
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
        varR = np.sqrt(8 * abs(rp.V0_IMAG) / rp.total_energy_range()) * r_best
        rc = np.array(
            (
                iosearch.readDataChem(
                    rp,
                    datafiles,
                    cutoff=r_best + r_cutoff * varR,
                    max_configs=max_configs,
                )
            ),
            dtype=object,
        )
    if len(rc) < 100 * rp.indyPars:
        return None, None

    rfacs, configs = rc[:, 0].astype(float), rc[:, 1]
    localizeFactor = kwargs["localize"]
    if localizeFactor == 0:
        localizeFactor = 1  # no localization
    # reduce to independent parameters during fit
    if not rp.domainParams:
        sps = [
            sp
            for sp in rp.searchpars
            if sp.el != "vac"
            and sp.mode != "dom"
            and sp.steps * localizeFactor >= 3
            and sp.linkedTo is None
            and sp.restrictTo is None
        ]
        indep_pars = np.array(
            [
                *np.array(
                    [
                        *np.array(configs, dtype=object)
                    ],  # TODO: pylint complains about too-many-function-args. I don't understand the line
                    dtype=object,
                ).reshape(-1, 1, 2)[:, 0, 1]
            ]
        )
        indep_pars = np.delete(
            indep_pars,
            [i for i in range(len(rp.searchpars) - 1) if rp.searchpars[i] not in sps],
            1,
        )
    else:
        # first the percentages:
        sps = [sp for sp in rp.searchpars if sp.mode == "dom"]
        reshaped = np.array(
            [*np.array(configs, dtype=object)], dtype=object
        ).reshape(  # TODO: pylint complains about too-many-function-args. I don't understand the line
            -1, len(rp.domainParams), 2
        )
        indep_pars = reshaped[:, :, 0].astype(int)  # contains the percentages
        # then the 'real' parameters:
        for j, dp in enumerate(rp.domainParams):
            new_sps = [
                sp
                for sp in dp.rpars.searchpars
                if sp.el != "vac"
                and sp.steps * localizeFactor >= 3
                and sp.linkedTo is None
                and sp.restrictTo is None
                and sp.mode != "dom"
            ]
            new_ip = np.array([*reshaped[:, j, 1]])
            new_ip = np.delete(
                new_ip,
                [
                    i
                    for i in range(len(dp.rpars.searchpars))
                    if dp.rpars.searchpars[i] not in new_sps
                ],
                1,
            )
            sps.extend(new_sps)
            indep_pars = np.append(indep_pars, new_ip, axis=1)
    indep_pars = indep_pars.astype(float)
    best_config = np.copy(indep_pars[np.argmin(rfacs)])
    sps_original = sps[:]

    which_regression = kwargs["type"].lower()
    alpha = kwargs["alpha"]
    if which_regression == "lasso":
        polyreg = make_pipeline(
            PolynomialFeatures(degree=2), Lasso(alpha=alpha, normalize=True)
        )  # TODO: incorrect keyword normalize
    elif which_regression == "ridge":
        polyreg = make_pipeline(
            PolynomialFeatures(degree=2), Ridge(alpha=alpha, normalize=True)
        )  # TODO: incorrect keyword normalize
    elif which_regression == "elasticnet":
        polyreg = make_pipeline(
            PolynomialFeatures(degree=2), ElasticNet(alpha=alpha, normalize=True)
        )  # TODO: incorrect keyword normalize
    else:
        if which_regression not in ("linearregression", "linear"):
            logger.warning(
                "Regression model {} not found, parabola fit "
                "defaulting to linear regression.".format(which_regression)
            )
        which_regression = "linearregression"
        polyreg = make_pipeline(PolynomialFeatures(degree=2), LinearRegression())
    xmin = np.copy(indep_pars[np.argmin(rfacs), :])
    rr = np.sqrt(8 * abs(rp.V0_IMAG) / rp.total_energy_range())
    ip_tofit = np.copy(indep_pars)
    rf_tofit = np.copy(rfacs)
    # throw out high R-factors - TODO: perhaps also throw out highest X% ?
    to_del = np.where(rf_tofit > min(rf_tofit) + 3 * rr)
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
        optimizerHelper, x0, args=(polyreg.predict,), method="L-BFGS-B", bounds=bounds
    )
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

    for i, sp in enumerate(sps):
        sp.parabolaFit["min"] = parab_min.x[i] + xmin[i]
        sp.parabolaFit["err_co"] = err_co[i]
        sp.parabolaFit["err_unco"] = err_unco[i]
        best_config[sps_original.index(sp)] = parab_min.x[i] + xmin[i]
    for sp in [
        sp
        for sp in rp.searchpars
        if isinstance(sp.linkedTo, SearchPar) or isinstance(sp.restrictTo, SearchPar)
    ]:
        if isinstance(sp.linkedTo, SearchPar):
            sp.parabolaFit = copy.copy(sp.linkedTo.parabolaFit)
        elif isinstance(sp.restrictTo, SearchPar):
            sp.parabolaFit = copy.copy(sp.restrictTo.parabolaFit)
    return (best_config, predictR)
