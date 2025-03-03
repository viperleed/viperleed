"""Functions for writing the SearchProgress.pdf and SearchReport.pdf files."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-19'
__license__ = 'GPLv3+'

import logging

import numpy as np

from viperleed.calc.lib.matplotlib_utils import CAN_PLOT
from viperleed.calc.lib.matplotlib_utils import close_figures
from viperleed.calc.lib.matplotlib_utils import prepare_matplotlib_for_calc
from viperleed.calc.lib.matplotlib_utils import skip_without_matplotlib

if CAN_PLOT:
    prepare_matplotlib_for_calc()
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.markers import MarkerStyle


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # Mute matplotlib debug messages                 # TODO: perhaps nicer to use at_level only in the relevant spots? See also iorfactor and ivplot


@skip_without_matplotlib
def writeSearchProgressPdf(rp, gens, rfacs, lastconfig,
                           outname="Search-progress.pdf",
                           csvname="Search-progress.csv",
                           markers=None,
                           rfac_predict=None):
    """
    Writes a pdf file with reports on R-factor convergence and current
    parameter scatter. Also writes a csv file containing the most basic
    information.

    Parameters
    ----------
    rp : Rparams
        The run parameters
    gens : list of int
        List of generation numbers for which data is stored.
    rfacs : list of (list of float)
        The R-factors for a given generation (corresponding to gens)
    lastconfig : tuple
        Lists (percent, dc) for each domain, with dc the parameter values of
        that domain
    outname : str, optional
        The file name to write to. The default is "Search-progress.pdf".
    csvname : str, optional
        The file name of the csv file to write to. The default is
        "Search-progress.csv".
    markers : list of tuples (int, str), optional
        Pass x-positions in the generations plots to be marked by vertical
        lines as tuples (gen, label), where gen is the position in the plot.
    rfac_predict : list of float, optional
        List of r-factor values determined from the parabola fit, by
        generation. The default is [].

    Returns
    -------
    None.

    """
    markers = markers or []

    figsPerPage = 5
    parsPerFig = 8

    searchname = rp.disp_blocks[rp.search_index][1]

    rfacsMin = np.array([min(rfa) for rfa in rfacs])
    rfacsMax = np.array([max(rfa) for rfa in rfacs])
    rfacsMean = np.array([np.mean(rfa) for rfa in rfacs])
    rp.searchplots[-1] = (searchname, gens, rfacsMin, rfacsMax, rfacsMean)
    # rfacsStd = np.array([np.std(rfa, ddof = 1) for rfa in rfacs])
    rlastunique = [rfacs[-1][0]]
    lastpops = [1]
    for r in rfacs[-1][1:]:
        if r == rlastunique[-1]:
            lastpops[-1] += 1
        else:
            rlastunique.append(r)
            lastpops.append(1)
    allcolors = []
    if rlastunique[-1] == rlastunique[0]:
        colors = [(0., 0., 0.)] * len(rlastunique)  # black
        allcolors = [(0., 0., 0.)] * len(rfacs[-1])
    else:
        colors = []
        for (i, r) in enumerate(rlastunique):
            w = (r - rlastunique[0]) / (rlastunique[-1] - rlastunique[0])
            w = np.sqrt(w)   # stronger scaling towards red
            colors.append((w, 0., 0.))
            for j in range(0, lastpops[i]):
                allcolors.append((w, 0., 0.))
    if (not rp.rfacscatter_all) or (rp.rfacscatter_all[-1][0] != gens[-1]):
        rp.storeRfacScatter([gens[-1]]*len(rlastunique), rlastunique,
                            lastpops, colors)
    deltagens = {"all": {1: 0}, "best": {1: 0}, "dec": {1: 0}}
    previous = {"all": rfacs[0][:],
                "best": min(rfacs[0]),
                "dec": rfacs[0][:(len(rfacs) // 10 + 1)]}
    for i in range(1, len(gens)):
        if min(rfacs[i]) != previous["best"]:
            deltagens["best"][gens[i]] = (gens[i]
                                          - max(deltagens["best"].keys()))
            previous["best"] = min(rfacs[i])
        if list(rfacs[i]) != list(previous["all"]):
            deltagens["all"][gens[i]] = (gens[i]
                                         - max(deltagens["all"].keys()))
            previous["all"] = rfacs[i][:]
        if list(rfacs[i][:(len(rfacs[i]) // 10 + 1)]) != list(previous["dec"]):
            deltagens["dec"][gens[i]] = (gens[i]
                                         - max(deltagens["dec"].keys()))
            previous["dec"] = rfacs[i][:(len(rfacs[i]) // 10 + 1)]

    figsize = (5.8, 8.3)
    figs = []

    # CSV output
    sep = "; "
    width = 12
    titles = ["Generation", "G_Delta_all", "G_Delta_dec", "G_Delta_best",
              "R_min", "R_max", "R_mean"]
    output = ""
    for t in titles:
        output += t.rjust(width)+sep
    output = output[:-len(sep)] + "\n"
    mc = markers[:]
    for i in range(0, len(gens)):
        if mc:
            if gens[i] > mc[0][0]:
                output += mc[0][1] + "\n"   # comment line for marker
                mc.pop(0)
        output += str(gens[i]).rjust(width) + sep
        for k in ["all", "dec", "best"]:
            if gens[i] in deltagens[k]:
                output += str(deltagens[k][gens[i]]).rjust(width) + sep
            else:
                output += "NaN".rjust(width) + sep
        for ls in [rfacsMin, rfacsMax, rfacsMean]:
            output += "{:.4f}".format(ls[i]).rjust(width)+sep
        output = output[:-len(sep)] + "\n"
    try:
        with open(csvname, "w") as wf:
            wf.write(output)
    except KeyboardInterrupt:
        raise
    except Exception:
        logger.warning("Failed to write "+csvname)

    # R-FACTOR AND GENERATION DELTA
    # create figure
    fig, (rfp, dgp) = plt.subplots(2, 1, sharex=True, figsize=figsize)
    dgp.set_xlabel('Generations')
    rfp.set_ylabel('R-Factor')
    dgp.set_ylabel('Generations since last change')
    # plot markers
    part = 0
    for (i, g) in enumerate(gens):
        if g > gens[-1] * 0.2:
            part = len(gens) - i
            break
    part = max(100, part)
    rfmin, rfmax = min(rfacsMin[-part:]), max(rfacsMax[-part:])
    if rfac_predict:
        (pred_x, pred_y) = tuple(zip(*rfac_predict))
        lowbound = max(min(min(pred_y[int(len(pred_y)*0.5):]), rfmin), 0)
    else:
        lowbound = rfmin
    rYrange = [lowbound-(rfmax-lowbound)*0.1, rfmax+(rfmax-rfmin)*0.1]
    if abs(rYrange[1] - rYrange[0]) < 1e-5:
        same = rYrange[0]
        rYrange[0] -= 0.05*same
        rYrange[1] += 0.05*same
    labely = rYrange[0] + (rYrange[1]-rYrange[0])*0.99
    xoff = gens[-1]*0.005
    for (xpos, label) in markers:
        rfp.axvline(x=xpos, lw=0.5, c="black")
        dgp.axvline(x=xpos, lw=0.5, c="black")
        rfp.text(xpos + xoff, labely, label, rotation=-90,
                 verticalalignment="top", size=4)
    # plot data
    if rfac_predict:
        rfp.step(pred_x, pred_y, where="post", color="seagreen",
                 label="Prediction")
    rfp.plot(gens, rfacsMin, '-', color='black', label="Best")
    rfp.fill_between(gens, rfacsMin, rfacsMax, facecolor='grey',
                     alpha=0.2, label="Range")
    rfp.plot(gens, rfacsMean, label="Mean")
    (x, y, s, c) = list(zip(*rp.rfacscatter))
    rfp.scatter(x, y, s=s, c=c)
    scatcol = {"all": "black", "dec": "tab:blue", "best": "tab:red"}
    labels = {"all": "Changes to any", "dec": "Changes in best 10%",
              "best": "Changes to best"}
    sizes = {"all": 4, "dec": 5, "best": 6}
    maxy = 1
    for k in ["all", "dec", "best"]:
        x, y = deltagens[k].keys(), deltagens[k].values()
        dgp.scatter(x, y, s=sizes[k], c=scatcol[k])
        maxdgx, maxdgy = [1], [0]
        for g in deltagens[k].keys():
            if deltagens[k][g] >= maxdgy[-1]:
                maxdgx.append(g)
                maxdgy.append(deltagens[k][g])
        if gens[-1] - max(deltagens[k].keys()) >= maxdgy[-1]:
            maxdgx.append(gens[-1])
            maxdgy.append(gens[-1] - max(deltagens[k].keys()))
        dgp.plot(maxdgx, maxdgy, '-', c=scatcol[k], label=labels[k])
        maxy = max(maxy, max(maxdgy), max(y))
    # layout
    rfp.set_ylim(rYrange)
    dgp.set_ylim([0, maxy*1.1])
    dgp.ticklabel_format(axis="x", style="sci", scilimits=(0, 4))
    fig.tight_layout()
    rfp.legend(loc="lower left")
    dgp.legend(loc="upper left")
    figs.append(fig)

    # SEARCH PARAMETER SCATTER
    fig, axs = plt.subplots(figsPerPage, figsize=figsize, squeeze=True)
    figcount = 0
    labels = {"geo": "GEOMETRY", "vib": "VIBRATION",
              "occ": "OCCUPATION", "dom": "DOMAIN AREAS"}
    offsets = []
    rpToDo = [rp]
    if rp.domainParams:
        rpToDo.extend([dp.rp for dp in rp.domainParams])
    for (k, crp) in enumerate(rpToDo):
        sps = [sp for sp in crp.searchpars if sp.el != "vac" and sp.steps > 1]
        if not rp.domainParams:
            confindex = 0
        elif crp != rp:
            confindex = k-1
        for mode in ["dom", "geo", "vib", "occ"]:
            if rp.domainParams and crp == rp and mode != "dom":
                continue
            if not rp.domainParams and mode == "dom":
                continue
            spm = [sp for sp in sps if sp.mode == mode]
            while len(spm) > 0:
                if figcount >= figsPerPage:
                    fig.tight_layout()
                    figs.append(fig)
                    fig, axs = plt.subplots(figsPerPage, figsize=figsize,
                                            squeeze=True)
                    figcount = 0
                plotpars = spm[:parsPerFig]
                spm = spm[parsPerFig:]
                title = labels[mode]
                addinfo = []
                if mode != "dom" and rp.domainParams:
                    addinfo.append("domain {}"
                                   .format(rp.domainParams[k-1].name))
                if len(crp.disp_blocks) > 1:
                    addinfo.append("search {}".format(searchname[:20]))
                if addinfo:
                    title += " ("
                    while addinfo:
                        title += addinfo.pop(0)
                        if addinfo:
                            title += ", "
                        else:
                            title += ")"
                axs[figcount].set_title(title, pad=8)
                pltpoints = []  # x, y, color, size, alpha
                bestpoints = []
                xlabels = []
                predict = []
                for (i, par) in enumerate(plotpars):
                    vals = []
                    for (j, conf) in enumerate(lastconfig):
                        if mode != "dom":
                            val = ((conf[confindex][1][crp.searchpars
                                                       .index(par)]-1)
                                   / (par.steps-1))
                        else:
                            val = conf[i][0] / 100
                        r, g, b = allcolors[j]
                        if par.linkedTo is None and par.restrictTo is None:
                            alpha = 1.0
                            vals.append(val)
                        else:
                            alpha = 0.5
                        color = (r, g, b, alpha)
                        pltpoints.append((i+1, val, color, 1))
                        if j == 0:
                            bestpoints.append((i+1, val, alpha))
                    if mode == "dom":
                        xlabels.append("#{}\n{}".format(
                            i+1, rp.domainParams[i].name))
                        edgetext = ["0%", "100%"]
                    else:
                        if mode != "occ":
                            el = par.el
                        else:
                            el = par.atom.el
                        xlabels.append(f'#{par.atom.num}\n{el}')
                        edgetext = ["", ""]
                        if isinstance(par.edges[0], (np.floating, float)):
                            edgetext = [str(round(v, 4)) for v in par.edges]
                        elif type(par.edges[0]) == np.ndarray:
                            edgetext = ["[" + ", ".join([str(round(f, 4))
                                                         for f in (
                                                    v * np.array([1, 1, -1]))])
                                        + "]" for v in par.edges]
                        elif type(par.edges[0]) == str:
                            edgetext = par.edges
                    axs[figcount].text(i+1, -0.03, edgetext[0], fontsize=3,
                                       ha="center", va="top")
                    axs[figcount].text(i+1, 1.02, edgetext[1], fontsize=3,
                                       ha="center", va="bottom")
                    if vals:
                        offsets.append(np.std(vals))
                        # mean = np.mean(vals)
                        # offsets.append(np.mean([abs(v - mean) for v in vals])
                        #                * 2)
                    check = ("err_co", "err_unco")
                    if (par.parabolaFit["min"] is not None and
                            not all(np.isnan(par.parabolaFit[s])
                                    for s in check)):
                        # alpha = 0.5
                        # if par.linkedTo is None and par.rstrictTo is None:
                        #     alpha = 1.0
                        a = alpha
                        n = (par.steps - 1)
                        good_fit = True
                        if not any(any(0 < ((par.parabolaFit["min"] - 1)/n
                                            + sign*par.parabolaFit[s]/n) < 1
                                       for sign in (+1, -1)) for s in check):
                            good_fit = False
                            a = 0.5
                        err_unco = par.parabolaFit["err_unco"] / n
                        err_co = par.parabolaFit["err_co"] / n
                        if np.isnan(err_unco):
                            err_unco = 2
                        if np.isnan(err_co):
                            err_co = 2
                        predict.append((i+1,
                                        (par.parabolaFit["min"] - 1) / n,
                                        err_unco,
                                        err_co,
                                        a, good_fit))

                # combine duplicates:
                i = 0
                while i < len(pltpoints):
                    x, y, c, s = pltpoints[i]
                    j = i+1
                    while j < len(pltpoints):
                        if pltpoints[j][0] == x and pltpoints[j][1] == y:
                            s += pltpoints[j][3]
                            if pltpoints[j][2][0] < c[0]:
                                # use the color closer to black
                                c = (pltpoints[j][2][0], c[1], c[2], c[3])
                            pltpoints.pop(j)
                        else:
                            j += 1
                    pltpoints[i] = (x, y, c, s)
                    i += 1
                if predict:
                    m = MarkerStyle("D")
                    m._transform.scale(1.0, 0.5)
                    err_off = 0.08  # error bar offset
                    for alpha in set([p[4] for p in predict]):
                        pred_ok = [p for p in predict
                                   if p[5] and p[4] == alpha]
                        pl1 = tuple(zip(*pred_ok))
                        pl2 = tuple(zip(*[p for p in predict
                                          if p not in pred_ok
                                          and p[4] == alpha]))
                        for plx, c in ((pl1, "seagreen"), (pl2, "red")):
                            if not plx:
                                continue
                            axs[figcount].scatter(plx[0], plx[1], color=c,
                                                  alpha=alpha, marker=m, s=50)
                        # errors
                        for ind in (2, 3):
                            pred_ok = [p for p in predict if
                                       (not np.isnan(p[ind]) and p[ind] < 1)
                                       and p[5] and p[4] == alpha
                                       and any(0 < p[1] + sign*p[ind] < 1
                                               for sign in (+1, -1))]
                            pl1 = tuple(zip(*pred_ok))
                            pl2 = tuple(zip(*[p for p in predict
                                              if p not in pred_ok
                                              and p[5] and p[4] == alpha]))
                            for plx, c in ((pl1, "seagreen"), (pl2, "red")):
                                if not plx:
                                    continue
                                off = -err_off if ind == 2 else err_off
                                px = np.array(plx[0], dtype=float) + off
                                axs[figcount].errorbar(
                                    px, plx[1], plx[ind], capsize=2,
                                    fmt="none", color=c, alpha=alpha,
                                    lw=0.5)
                x, y = [p[0] for p in pltpoints], [p[1] for p in pltpoints]
                c, s = [p[2] for p in pltpoints], [p[3] for p in pltpoints]
                axs[figcount].plot([0, parsPerFig+2], [0.5, 0.5], color='grey',
                                   alpha=0.2)
                axs[figcount].scatter(x, y, s=s, c=c)
                for (px, py, alpha) in bestpoints:
                    axs[figcount].annotate("", (px+0.05, py), (px+0.25, py),
                                           arrowprops=dict(arrowstyle="wedge",
                                                           facecolor="black",
                                                           alpha=alpha))
                    axs[figcount].annotate("", (px-0.05, py), (px-0.25, py),
                                           arrowprops=dict(arrowstyle="wedge",
                                                           facecolor="black",
                                                           alpha=alpha))
                xlabels.extend([""] * (parsPerFig - len(xlabels)))
                axs[figcount].set_xlim([0, parsPerFig+1])
                axs[figcount].set_ylim([0, 1])
                axs[figcount].set_xticks(list(range(1, parsPerFig+1)))
                axs[figcount].set_xticklabels(xlabels)
                axs[figcount].tick_params(which='both', top=False,
                                          bottom=False, left=False,
                                          right=False, labelleft=False,
                                          labelbottom=True, pad=4)
                figcount += 1
    if offsets:
        rp.parScatter[-1].append((gens[-1], np.mean(offsets), max(offsets)))
    for i in range(figcount, figsPerPage):
        axs[i].axis('off')
    if fig not in figs:
        fig.tight_layout()
        figs.append(fig)

    # save
    if searchname in rp.lastParScatterFigs:
        close_figures(plt, *rp.lastParScatterFigs[searchname])
    rp.lastParScatterFigs[searchname] = figs[1:]
    try:
        pdf = PdfPages(outname)
        for fig in figs:
            pdf.savefig(fig)
    except PermissionError:
        logger.warning("Failed to write to " + outname
                       + ": Permission denied.")
    except KeyboardInterrupt:
        raise
    except Exception:
        logger.warning("Failed to write to "+outname)
        raise
    finally:
        try:
            pdf.close()
        except Exception:
            pass
    figures = (f for f in figs
               if searchname not in rp.lastParScatterFigs
               or f not in rp.lastParScatterFigs[searchname])
    close_figures(plt, *figures)


@skip_without_matplotlib
def writeSearchReportPdf(rp, outname="Search-report.pdf"):
    """
    Writes a pdf file with reports on R-factor convergence and parameter
    scatter, collated over the entire run (i.e. potentially multiple searches).

    Parameters
    ----------
    rp : Rparams
        The run parameters
    outname : str, optional
        The file name to write to. The default is "Search-report.pdf".

    Returns
    -------
    None.

    """
    allmin = []
    allmax = []
    allmean = []
    allgens = []
    markers = []
    parScatterLines = []  # list of lists [gens, mean, max] per search
    gencount = 0
    for i in range(0, len(rp.searchplots)):
        (name, gens, rmin, rmax, rmean) = rp.searchplots[i]
        markers.append((gencount, "Search "+name))
        allgens.extend([v + gencount for v in gens])
        allmin.extend(rmin)
        allmax.extend(rmax)
        allmean.extend(rmean)
        if rp.parScatter[i]:
            parScatterLines.append(list(zip(*rp.parScatter[i])))
            parScatterLines[-1][0] = [v + gencount for v in
                                      parScatterLines[-1][0]]
        gencount = allgens[-1]

    figsize = (5.8, 8.3)
    figs = []
    # R-FACTORS AND MEAN SCATTER
    # create figure
    fig, (rfp, msp) = plt.subplots(2, 1, sharex=True, figsize=figsize)
    msp.set_xlabel('Generations')
    rfp.set_ylabel('R-Factor')
    msp.set_ylabel('Parameter scatter')
    # plot markers
    part = 0
    for (i, g) in enumerate(allgens):
        if g > allgens[-1] * 0.2:
            part = len(allgens) - i
            break
    part = max(50, part)

    # Notice that we take max(mean) rather than max(max) in order
    # to zoom-in in the interesting region. Using max(max) would
    # normally obscure this region, as the worst R could be
    # seriously worse than the mean.
    rfmin, rfmax = min(allmin[-part:]), max(allmean[-part:])
    if rfmax <= rfmin:
        rfmin *= 0.95
        rfmax *= 1.05
    rYrange = [rfmin-(rfmax-rfmin)*0.1, rfmax+(rfmax-rfmin)*0.1]
    rYrange[1] = max(rYrange[1], max(allmin)+(rfmax-rfmin)*0.1)

    labely = rYrange[0] + (rYrange[1]-rYrange[0])*0.99
    xoff = allgens[-1]*0.005
    for (xpos, label) in markers:
        rfp.axvline(x=xpos, lw=0.5, c="black")
        msp.axvline(x=xpos, lw=0.5, c="black")
        rfp.text(xpos + xoff, labely, label, rotation=-90,
                 verticalalignment="top", size=4)
    # plot data

    rfp.fill_between(allgens, allmin, allmax, facecolor='grey',
                     alpha=0.2, label="Range")
    rfp.plot(allgens, allmean, label="Mean")
    rfp.plot(allgens, allmin, '-', color='black', label="Best")

    labelled = False
    scattermax = 0
    for (allgens, psmean, psmax) in parScatterLines:
        meanline, = msp.plot(allgens, psmean, '-', color="tab:blue")
        maxline, = msp.plot(allgens, psmax, '-', color="black")
        scattermax = max(scattermax, max(psmax))
        if not labelled:
            meanline.set_label('Mean parameter \u03C3')    # sigma
            maxline.set_label('Highest parameter \u03C3')  # sigma
            labelled = True
    # Avoid matplotlib warnings if scattermax remained == 0. This
    # simply means that the all the individuals in the population
    # converged to the same configuration.
    if scattermax <= 1e-5:
        scattermax = 0.05

    # layout
    rfp.set_ylim(rYrange)
    msp.set_ylim([0, scattermax*1.1])
    msp.ticklabel_format(axis="x", style="sci", scilimits=(0, 4))
    fig.tight_layout()
    rfp.legend(loc="lower left")
    msp.legend(loc="lower left")

    # ADD SCATTERS
    for k in rp.lastParScatterFigs.keys():
        figs.append(rp.lastParScatterFigs[k])
    # save
    try:
        pdf = PdfPages(outname)
        pdf.savefig(fig)
        for k in rp.lastParScatterFigs.keys():
            for f in rp.lastParScatterFigs[k]:
                pdf.savefig(f)
    except PermissionError:
        logger.warning("Failed to write to " + outname
                       + ": Permission denied.")
    except KeyboardInterrupt:
        raise
    except Exception:
        logger.warning("Failed to write to "+outname)
        raise
    finally:
        try:
            pdf.close()
        except Exception:
            pass
    close_figures(plt, fig)
