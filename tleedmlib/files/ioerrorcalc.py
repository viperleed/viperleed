# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 17:22:20 2021

@author: Florian Kraushofer

Functions for reading and writing files relevant to the error calculation
"""

import numpy as np
import logging
import re
from scipy import interpolate
from viperleed.tleedmlib.base import range_to_str

try:
    import matplotlib
    matplotlib.rcParams.update({'figure.max_open_warning': 0})
    matplotlib.use('Agg')  # !!! check with Michele if this causes conflicts
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    # import matplotlib.ticker as plticker
    matplotlib.rcParams["mathtext.default"] = "regular"
except Exception:
    plotting = False
else:
    plotting = True

logger = logging.getLogger("tleedm.files.ioerrorcalc")
logger.setLevel(logging.INFO)


def write_errors_csv(errors, filename="Errors.csv", sep=";"):
    """
    Writes errors from the error calculation into a CSV file

    Parameters
    ----------
    errors : R_Error
        Data structure from sections.errorcalc containing the information
    filename : str, optional
        Name of the csv file to write. The default is "Errors.csv".
    sep : str, optional
        The separator to use. The default is ";".

    Returns
    -------
    None.

    """
    # contents of the columns; start with only titles:
    columns = {"at": ["Atoms"],
               "mode": ["Mode"],
               "dir": ["Disp. along (1st atom)"],
               "disp": ["Displacement [A]"],
               "rfac": ["R"]}
    for err in errors:
        ats = range_to_str([at.oriN for at in err.atoms])
        if isinstance(err.displacements[0], np.ndarray):
            dirvec = ((err.displacements[-1] - err.displacements[0])
                      * np.array([1, 1, -1]))
            dirvec = dirvec / np.linalg.norm(dirvec)
            direction = ("[" + ", ".join([str(round(f, 4))
                                          for f in dirvec]) + "]")
        for i in range(0, len(err.rfacs)):
            columns["at"].append(ats)
            columns["mode"].append(err.mode)
            if err.main_element:
                columns["dir"].append(err.main_element)
                disp = err.displacements[i]
            elif isinstance(err.displacements[i], np.ndarray):
                columns["dir"].append(direction)
                disp = np.dot(dirvec * np.array([1, 1, -1]),
                              err.displacements[i])
            else:
                columns["dir"].append("")
                disp = err.displacements[i]
            columns["disp"].append("{:.4f}".format(disp))
            columns["rfac"].append("{:.4f}".format(err.rfacs[i]))

    if all(s == "" for s in columns["dir"][1:]):
        del columns["dir"]

    widths = {}
    for k in columns:
        widths[k] = max(len(s) for s in columns[k]) + 1

    output = ""
    for i in range(len(next(iter(columns.values())))):
        for k in columns:
            output += columns[k][i].rjust(widths[k]) + sep
        output = output[:-1] + "\n"

    try:
        with open(filename, "w") as wf:
            wf.write(output)
    except Exception:
        logger.warning("Failed to write "+filename)
    return


def write_errors_pdf(errors, filename="Errors.pdf", var=None):
    global plotting
    if not plotting:
        logger.debug("Necessary modules for plotting not found. Skipping "
                     "error plotting.")
        return

    fig_order = (3, 2)
    figs_per_page = fig_order[0] * fig_order[1]
    figsize = (5.8, 8.3)
    figs = []

    titles = {"geo": "Geometry",
              "vib": "Vibrational amplitudes",
              "occ": "Site occupation"}

    rmin = min(r for err in errors for r in err.rfacs)
    for mode in ("geo", "vib", "occ"):
        mode_errors = [err for err in errors if err.mode == mode]
        if not mode_errors:
            continue
        err_x = {}
        err_y = {}
        err_x_to_mark = {}
        err_legend = {}
        err_disp = {}

        for err in mode_errors:
            ats = "Atom"
            if len(ats) > 1:
                ats += "s"
            ats += " " + range_to_str([at.oriN for at in err.atoms])
            if mode == "geo":
                dirvec = ((err.displacements[-1] - err.displacements[0])
                          * np.array([1, 1, -1]))
                dirvec = dirvec / np.linalg.norm(dirvec)
                direction = ("[" + ", ".join([str(round(f, 2))
                                              for f in dirvec]) + "]")
                disp = [np.dot(dirvec * np.array([1, 1, -1]), v) * 100
                        for v in err.displacements]
                err_legend[err] = ats + " along " + direction
            else:
                disp = err.displacements[:] * 100   # in pm
                err_legend[err] = ats
            err_disp[err] = disp
        xvals = [d for k in err_disp for d in err_disp[k]]
        xrange = [min(xvals) - abs(max(xvals) - min(xvals)) * 0.05,
                  max(xvals) + abs(max(xvals) - min(xvals)) * 0.05]
        for err in mode_errors:
            disp = err_disp[err]
            tck = interpolate.splrep(disp, err.rfacs)
            x = np.arange(min(xvals), max(xvals)+0.01,
                          (xrange[1] - xrange[0])*1e-4)
            y = interpolate.splev(x, tck)
            ind_to_mark = [np.argmin(abs(x - v)) for v in disp]
            err_x_to_mark[err] = ind_to_mark
            err_x[err] = x
            err_y[err] = y
        # plot combined figure
        rmax = max(r for err in mode_errors for r in err.rfacs)
        fig = plt.figure(figsize=(5.8, 5.8))
        ax = fig.add_subplot(1, 1, 1)
        if mode != "occ":
            ax.set_xlabel('Deviation from bestfit value (pm)')
        else:
            ax.set_xlabel('Site occupation (%)')
        ax.set_ylabel('Pendry R-factor')
        ax.set_title(titles[mode])
        if var and rmin + var < rmax + (rmax-rmin)*0.1:
            ax.plot(xrange, [rmin + var]*2, color="slategray")
            ax.text((xrange[0] + xrange[1])/2, rmin + var + (rmax-rmin)*0.01,
                    "$R_P + var(R_P)$", ha="center", va="bottom")
        for err in mode_errors:
            ax.plot(err_x[err], err_y[err], '-o', label=err_legend[err],
                    markevery=err_x_to_mark[err])
        ax.set_xlim(*xrange)
        ax.set_ylim(rmin - ((rmax-rmin)*0.1), rmax + ((rmax-rmin)*0.1))
        # ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.legend(fontsize="small")
        fig.tight_layout()
        figs.append(fig)
        # now plot individual figures
        figcount = 0
        fig, axs = plt.subplots(fig_order[0], fig_order[1],
                                figsize=figsize, squeeze=True)
        axs = axs.flatten()
        while len(mode_errors) > 0:
            if figcount >= figs_per_page:
                fig.tight_layout(rect=(0, 0, 1, 0.965))
                fig.suptitle(titles[mode])
                figs.append(fig)
                figcount = 0
                fig, axs = plt.subplots(fig_order[0], fig_order[1],
                                        figsize=figsize, squeeze=True)
                axs = axs.flatten()
            err = mode_errors.pop(0)
            rmax = max(r for r in err.rfacs)
            xvals = err_disp[err]
            xrange = [min(xvals) - abs(max(xvals) - min(xvals)) * 0.05,
                      max(xvals) + abs(max(xvals) - min(xvals)) * 0.05]
            if var and rmin + var < rmax + (rmax-rmin)*0.1:
                axs[figcount].plot(xrange, [rmin + var]*2, color="slategray",
                                   linewidth=1)
                axs[figcount].text((xrange[0] + xrange[1])/2,
                                   rmin + var + (rmax-rmin)*0.01,
                                   "$R_P + var(R_P)$", fontsize=6,
                                   ha="center", va="bottom")
            axs[figcount].plot(err_x[err], err_y[err], '-o',
                               markevery=err_x_to_mark[err],
                               linewidth=1, ms=2,
                               label=re.sub("along", "\nalong",
                                            err_legend[err]))
            axs[figcount].set_xlim(*xrange)
            ylim = (rmin - ((rmax-rmin)*0.1), rmax + ((rmax-rmin)*0.1))
            axs[figcount].set_ylim(ylim)
            # determine yticks
            n_ticks = 0
            dec = 2   # tick decimals
            while n_ticks == 0:
                tickbounds = (np.ceil(ylim[0]*(10**dec))/(10**dec),
                              np.floor(ylim[1]*(10**dec))/(10**dec))
                n_ticks = (int((10**dec)*(tickbounds[1] - tickbounds[0]
                                          + 1e-6)) + 1)
                if n_ticks > 8:
                    n_ticks = n_ticks // (n_ticks // 4)
                if n_ticks < 2:
                    dec += 1
                    n_ticks = 0
            tickstep = (tickbounds[1] - tickbounds[0]) / (n_ticks - 1)
            tickstep = np.floor((tickstep+1e-6) * (10**dec)) / (10**dec)
            yticks = (list(np.arange(tickbounds[0], tickbounds[1]+1e-6,
                                     tickstep))[:n_ticks])
            axs[figcount].set_yticks(yticks)
            axs[figcount].set_yticklabels([f"{v:.{dec}f}" for v in yticks])
            axs[figcount].xaxis.set_major_locator(plt.MaxNLocator(5))
            axs[figcount].tick_params(labelsize=6)
            if mode != "occ":
                axs[figcount].set_xlabel('Deviation from bestfit value (pm)',
                                         fontsize=8)
            else:
                axs[figcount].set_xlabel('Site occupation (%)', fontsize=8)
            axs[figcount].set_ylabel('Pendry R-factor', fontsize=8)
            axs[figcount].legend(fontsize="x-small", frameon=False)
            # axs[figcount].set_title(err_legend[err], fontsize=8)
            figcount += 1
        for ax in axs[figcount:]:
            ax.axis('off')
        fig.tight_layout(rect=(0, 0, 1, 0.965))
        fig.suptitle(titles[mode])
        figs.append(fig)
    try:
        pdf = PdfPages(filename)
        for fig in figs:
            pdf.savefig(fig)
    except PermissionError:
        logger.warning("Failed to write to " + filename
                       + ": Permission denied.")
    except KeyboardInterrupt:
        raise
    except Exception:
        logger.warning("Failed to write to "+filename, exc_info=True)
    finally:
        try:
            pdf.close()
        except Exception:
            pass
    for fig in figs:
        try:
            plt.close(fig)
        except Exception:
            pass
