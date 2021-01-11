# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 13:43:36 2020

@author: Florian Kraushofer

Functions for writing the SearchProgress.pdf and SearchReport.pdf files.
"""

import numpy as np
import logging
from matplotlib.markers import MarkerStyle

logger = logging.getLogger("tleedm.files.searchpdf")
logger.setLevel(logging.INFO)

try:
    import matplotlib
    matplotlib.rcParams.update({'figure.max_open_warning': 0})
    matplotlib.use('Agg')  # !!! check with Michele if this causes conflicts
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
except:
    plotting = False
else:
    plotting = True

def writeSearchProgressPdf(rp, gens, rfacs, lastconfig, 
                           outname = "Search-progress.pdf",
                           csvname = "Search-progress.csv",
                           markers = [],
                           rfac_predict = []):
    global plotting
    if not plotting:
        return None
    
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
        colors = [(0.,0.,0.,1.)] * len(rlastunique)  #black
        allcolors = [(0.,0.,0.,1.)] * len(rfacs[-1])
    else:
        colors = []
        for (i, r) in enumerate(rlastunique):
            w = (r - rlastunique[0]) / (rlastunique[-1] - rlastunique[0])
            w = np.sqrt(w)   # stronger scaling towards red
            colors.append((w, 0., 0.))
            for j in range(0,lastpops[i]):
                allcolors.append((w, 0., 0., 1.))
    if (not rp.rfacscatter_all) or (rp.rfacscatter_all[-1][0] != gens[-1]):
        rp.storeRfacScatter([gens[-1]]*len(rlastunique), rlastunique, 
                            lastpops, colors)
    deltagens = [gens[n] - gens[n-1] for n in range(1, len(gens))]
    maxdgx, maxdgy = [gens[1]], [deltagens[0]]
    for i in range(0, len(deltagens)):
        if deltagens[i] >= maxdgy[-1]:
            maxdgx.append(gens[i+1])
            maxdgy.append(deltagens[i])

    figsize = (5.8, 8.3)
    figs = []
    
    # CSV output
    sep = "; "
    width = 12
    titles = ["Generation", "Gen_Delta", "R_min", "R_max", "R_mean"]
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
        for l in [gens, [0]+deltagens]:
            output += str(l[i]).rjust(width)+sep
        for l in [rfacsMin, rfacsMax, rfacsMean]:
            output += "{:.4f}".format(l[i]).rjust(width)+sep
        output = output[:-len(sep)] + "\n"
    try:
        with open(csvname, "w") as wf:
            wf.write(output)
    except KeyboardInterrupt:
        raise
    except:
        logger.warning("Failed to write "+csvname)
    
    # R-FACTOR AND GENERATION DELTA
    # create figure
    fig, (rfp, dgp) = plt.subplots(2, 1, sharex=True, figsize = figsize)
    dgp.set_xlabel('Generations')
    rfp.set_ylabel('R-Factor')
    dgp.set_ylabel('Generation delta')
    #plot markers
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

    labely = rYrange[0] + (rYrange[1]-rYrange[0])*0.99
    xoff = gens[-1]*0.005
    for (xpos, label) in markers:
        rfp.axvline(x = xpos, lw = 0.5, c = "black")
        dgp.axvline(x = xpos, lw = 0.5, c = "black")
        rfp.text(xpos + xoff, labely, label, rotation=-90, 
                 verticalalignment="top", size=4)
    # plot data
    if rfac_predict:
        rfp.step(pred_x, pred_y, where="post", color="seagreen", 
                 label = "Prediction")
    rfp.plot(gens, rfacsMin, '-', color='black', label = "Best")
    rfp.fill_between(gens, rfacsMin, rfacsMax, facecolor='grey', 
                     alpha=0.2, label="Range")
    rfp.plot(gens, rfacsMean, label = "Mean")
    # rfp.fill_between(gens, rfacsMean - rfacsStd/2, rfacsMean + rfacsStd/2, 
    #                  alpha=0.2, label="Sigma")
    (x,y,s,c) = list(zip(*rp.rfacscatter))
    rfp.scatter(x, y, s=s, c=c)
    # rfp.scatter([gens[-1]]*len(rlastunique), rlastunique, 
    #             s = lastpops, c = colors)
    dgp.scatter(gens[1:], deltagens, s = 4, c='black', label="Every")
    dgp.plot(maxdgx, maxdgy, '-', color='royalblue', label="Max")
    # layout
    rfp.set_ylim(rYrange)
    dgp.set_ylim([0, max(deltagens)*1.1])
    dgp.ticklabel_format(axis="x", style="sci", scilimits=(0,4))
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
                axs[figcount].set_title(title)
                pltpoints = [] # x, y, color, size, alpha
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
                        r, g, b, alpha = allcolors[j]
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
                        xlabels.append("#{}\n{}".format(i+1, 
                                                      rp.domainParams[i].name))
                    else:
                        if mode != "occ":
                            el = par.el
                        else:
                            el = par.atom.el
                        xlabels.append("#{}\n{}".format(par.atom.oriN, el))
                    if vals:
                        offsets.append(np.std(vals))
                        # mean = np.mean(vals)
                        # offsets.append(np.mean([abs(v - mean) for v in vals])
                        #                * 2)
                    if par.parabolaFit["min"] is not None:
                        # alpha = 0.5
                        # if par.linkedTo is None and par.rstrictTo is None:
                        #     alpha = 1.0
                        predict.append((i+1, (par.parabolaFit["min"] - 1)
                                              / (par.steps - 1), alpha))
                        
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
                    pltpoints[i] = (x,y,c,s)
                    i += 1
                if predict:
                    m = MarkerStyle("D")
                    m._transform.scale(1.0, 0.5)
                    for alpha in set([p[2] for p in predict]):
                        (px, py, _) = tuple(zip(*[p for p in predict 
                                                  if p[2] == alpha]))
                        axs[figcount].scatter(px, py, color="seagreen", 
                                              alpha=alpha, marker=m, s=50)
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
                axs[figcount].set_xticks(list(range(1,parsPerFig+1)))
                axs[figcount].set_xticklabels(xlabels)
                axs[figcount].tick_params(axis='y', which='both', left=False, 
                                right=False, labelleft=False)
                figcount += 1
    if offsets:
        rp.parScatter[-1].append((gens[-1],np.mean(offsets),max(offsets)))
    for i in range(figcount, figsPerPage):
        axs[i].axis('off')
    if not fig in figs:
        fig.tight_layout()
        figs.append(fig)
        
    # save
    if searchname in rp.lastParScatterFigs:
        for f in rp.lastParScatterFigs[searchname]:
            try:
                plt.close(f)
            except:
                pass
    rp.lastParScatterFigs[searchname] = figs[1:]
    try:
        pdf = PdfPages(outname)
        for fig in figs:
            pdf.savefig(fig)
    except PermissionError:
        logger.warning("Failed to write to " + outname
                        + ": Permission denied.")
        raise
    except KeyboardInterrupt:
        raise
    except:
        logger.warning("Failed to write to "+outname)
        raise
    finally:
        try:
            pdf.close()
        except:
            pass
    for fig in [f for f in figs if 
                not searchname in rp.lastParScatterFigs or
                not f in rp.lastParScatterFigs[searchname]]:
        try:
            plt.close(fig)
        except:
            pass
    return None

def writeSearchReportPdf(rp, outname = "Search-report.pdf"):
    global plotting
    if not plotting:
        return None
    
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
    fig, (rfp, msp) = plt.subplots(2, 1, sharex=True, figsize = figsize)
    msp.set_xlabel('Generations')
    rfp.set_ylabel('R-Factor')
    msp.set_ylabel('Parameter scatter')
    #plot markers
    part = 0
    for (i, g) in enumerate(allgens):
        if g > allgens[-1] * 0.2:
            part = len(allgens) - i
            break
    part = max(50, part)
    rfmin, rfmax = min(allmin[-part:]), max(allmean[-part:])
    if rfmax == rfmin:
        rfmin *= 0.95
        rfmax *= 1.05
    rYrange = [rfmin-(rfmax-rfmin)*0.1, rfmax+(rfmax-rfmin)*0.1]
    rYrange[1] = min(rYrange[1], 2*max(allmin) - rYrange[0])
    
    labely = rYrange[0] + (rYrange[1]-rYrange[0])*0.99
    xoff = allgens[-1]*0.005
    for (xpos, label) in markers:
        rfp.axvline(x = xpos, lw = 0.5, c = "black")
        msp.axvline(x = xpos, lw = 0.5, c = "black")
        rfp.text(xpos + xoff, labely, label, rotation=-90, 
                 verticalalignment="top", size=4)
    # plot data
    
    rfp.fill_between(allgens, allmin, allmax, facecolor='grey', 
                     alpha=0.2, label="Range")
    rfp.plot(allgens, allmean, label = "Mean")
    rfp.plot(allgens, allmin, '-', color='black', label = "Best")

    labelled = False
    scattermax = 0
    for (allgens, psmean, psmax) in parScatterLines:
        meanline, = msp.plot(allgens, psmean, '-', color="tab:blue")
        maxline, = msp.plot(allgens, psmax, '-', color="black")
        scattermax = max(scattermax, max(psmax))
        if not labelled:
            meanline.set_label('Mean')
            maxline.set_label('Max')
            labelled = True

    # layout
    rfp.set_ylim(rYrange)
    msp.set_ylim([0, scattermax*1.1])
    msp.ticklabel_format(axis="x", style="sci", scilimits=(0,4))
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
        raise
    except KeyboardInterrupt:
        raise
    except:
        logger.warning("Failed to write to "+outname)
        raise
    finally:
        try:
            pdf.close()
        except:
            pass
    try:
        plt.close(fig)
    except:
        pass
    return None