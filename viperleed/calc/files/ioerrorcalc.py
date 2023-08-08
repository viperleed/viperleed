# -*- coding: utf-8 -*-
"""Functions for reading and writing files relevant to the error calculation.
"""
import numpy as np
import logging
import re
from scipy import interpolate
from zipfile import ZipFile, ZIP_DEFLATED

from viperleed.calc.lib.base import range_to_str, max_diff

__authors__ = ["Florian Kraushofer (@fkraushofer)",
               "Alexander M. Imre (@amimre)"]
__created__ = "2021-03-18"

try:
    import matplotlib
    matplotlib.rcParams.update({'figure.max_open_warning': 0})
    matplotlib.use('Agg')  # !!! check with Michele if this causes conflicts
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    # import matplotlib.ticker as plticker
    plt.style.use('viperleed.tleedm')
except Exception:
    _CAN_PLOT = False
else:
    _CAN_PLOT = True

logger = logging.getLogger("tleedm.files.ioerrorcalc")
logger.setLevel(logging.INFO)

def extract_var_r(errors):
    var_r_info = {
        "geo": None,
        "vib": None,
        "occ": None,
    }
    for mode in var_r_info.keys():
        mode_errors = [err for err in errors
                       if (err.mode == mode and err.r_type==1)]
        if mode_errors:
            var_r_info[mode] = mode_errors[0].var_r
    return var_r_info


def write_errors_summary_csv(summary_content, summary_path,
                             summary_fname="Errors_summary.csv"):
    try:
        with open(summary_path/summary_fname, "w", encoding="utf-8") as wf:
            wf.write(summary_content)
    except Exception as err:
        logger.error("Failed to write error calculation summary "
                     f"{summary_fname}:\n{err}")


def write_errors_archive(individual_files,
                         archive_path,
                         compression_level=2,
                         archive_fname="Errors.zip"):
    try:
        with ZipFile(archive_path/archive_fname, 'w',
                     compression=ZIP_DEFLATED,
                     compresslevel=compression_level) as err_archive:
            for fname, content in individual_files.items():
                err_archive.writestr(fname, content)
    except Exception as err:
        logger.error("Failed to write error calculation archive "
                     f"{archive_fname}:\n{err}")


def generate_errors_csv(errors, sep=","):

    summary_columns = {
        "at": ["Atoms"],
        "mode": ["Mode"],
        "dir": ["Direction"],
        "r_min" : ["R_min"],
        "var_r" : ["var(R)"],
        "p_min": ["p_min"],
        "d_p_minus": ["-Δp"],
        "d_p_plus": ["+Δp"],
    }

    # individual files
    indiv_files = {}

    for param_err in errors:

        atoms_str = range_to_str([at.oriN for at in param_err.atoms])
        summary_columns["at"].append(atoms_str)
        summary_columns["dir"].append(param_err.disp_label)
        summary_columns["mode"].append(param_err.mode)
        summary_columns["r_min"].append(param_err.get_r_min)
        summary_columns["var_r"].append(param_err.var_r)
        summary_columns["p_min"].append(param_err.get_p_min)

        # statistical error estimates
        d_p_minus, d_p_plus = param_err.get_error_estimates
        summary_columns["d_p_minus"].append(d_p_minus)
        summary_columns["d_p_plus"].append(d_p_plus)

        # create_filename for individual files
        indiv_fname = f"Errors_{param_err.mode}_atoms#{atoms_str}"
        if param_err.mode == "geo":
            indiv_fname += f"_{param_err.disp_label}"
        indiv_fname += ".csv"

        if param_err.mode == "geo":
            param_columns = geo_errors_csv_content(param_err)
        elif param_err.mode == "vib":
            param_columns = vib_errors_csv_content(param_err)
        elif param_err.mode == "occ":
            param_columns = occ_errors_csv_content(param_err)
        else:
            raise ValueError(f'Unknown mode "{param_err.mode}"')
        file_content = get_string_from_columns(param_columns, sep)
        indiv_files[indiv_fname] = file_content

    summary_csv_content = get_string_from_columns(summary_columns, sep)

    return summary_csv_content, indiv_files


def geo_errors_csv_content(error):
    """Generate columns dict for geometrical errors containing the
    contents of a file to be written into Errors.zip.

    Parameters
    ----------
    error : R_Error
        Error object for geometrical errors.

    Returns
    -------
    dict
        columns dict containing displacements and R-factors.
        
    Raises
    ------
    ValueError
        If error.mode is not "geo".
    """
    
    if error.mode != "geo":
        raise ValueError(f'Cannot format errors of type "{error.mode}"')
    columns = {
        "disp" : [f"Displacement ({error.disp_label}) [Å]"],
        "rfac" : ["R"]
    }
    rfacs = error.rfacs
    for line in range(0, len(error.rfacs)):
        columns["disp"].append(error.lin_disp[line])
        columns["rfac"].append(error.rfacs[line])
    return columns


def vib_errors_csv_content(error):
    """Generate columns dict for vibrational errors containing the
    contents of a file to be written into Errors.zip.

    Parameters
    ----------
    error : R_Error
        Error object for vibrational errors.

    Returns
    -------
    dict
        columns dict containing displacements and R-factors.
        
    Raises
    ------
    ValueError
        If error.mode is not "vib".
    """
    if error.mode != "vib":
        raise ValueError(f'Cannot format errors of type "{error.mode}"')
    columns = {
        "disp" : ["Vib. Amp. change [Å]"],
        "rfac" : ["R"]
    }
    rfacs = error.rfacs
    for line in range(0, len(error.rfacs)):
        columns["disp"].append(error.lin_disp[line])
        columns["rfac"].append(error.rfacs[line])
    return columns


def occ_errors_csv_content(error):
    """Generate columns dict for occupational errors containing the
    contents of a file to be written into Errors.zip.

    Parameters
    ----------
    error : R_Error
        Error object for occupational errors.

    Returns
    -------
    dict
        columns dict containing displacements and R-factors.
        
    Raises
    ------
    ValueError
        If error.mode is not "occ".
    """
    if error.mode != "occ":
        raise ValueError(f'Cannot format errors of type "{error.mode}"')
    columns = {}
    for elem in error.elem_occ.keys():
        columns[elem] = [f"Occupation {elem} [%]",]
    columns["rfac"] = ["R",]

    rfacs = error.rfacs
    for line in range(0, len(error.rfacs)):
        for elem, el_occ in error.elem_occ.items():
            columns[elem].append(el_occ[line]*100)  # in %
        columns["rfac"].append(error.rfacs[line])
    return columns


def get_string_from_columns(columns, sep=","):
    """Formats a columns dictionary into a string conforming to the CSV
    format.

    Parameters
    ----------
    columns : dict
        dict holding the contents to be written in the CSV. Keys are not
        used, values must be a list of entries for each column. Entries 
        can be str, int, float or None.
    sep : str
        CSV separator character to be used. Default is ",".

    Returns
    -------
    str
        str containing the formatted contents of the CSV file.
    """
    widths = {}
    for col in columns:
        widths[col] = max(len(
            format_col_content(col_content)
            ) for col_content in columns[col]) + 1

    content = ""
    for i in range(len(next(iter(columns.values())))):
        for col in columns:
            col_content = format_col_content(columns[col][i])
            content += col_content.rjust(widths[col]) + sep
        content = content[:-1] + "\n"
    return content


def format_col_content(content):
    """Formats a value into a string suitable for writing to errors CSV
    files.

    Parameters
    ----------
    content : str, int, float or None
        _description_

    Returns
    -------
    str
        Formatted string. If float, 4 decimal places are used; if
        content is None, return "N/A".

    Raises
    ------
    ValueError
        If content is neither, str, int, float or None.
    """
    if isinstance(content, str):
        return content
    elif isinstance(content, int):
        return str(content)
    elif isinstance(content, float):
        return f"{content:.4f}"
    elif content is None:
        return "N/A"
    else:
        raise ValueError("Cannot format value for writing to Errors CSV"
                         f" file: {content}")


def make_errors_figs(errors):
    """Creates figures for Errors.pdf.

    Parameters
    ----------
    errors : list of R_Error
        contains the R-factors to be plotted
    filename : str, optional
        Path of file to be written, by default "Errors.pdf"
    """
    global _CAN_PLOT
    if not _CAN_PLOT:
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

    for mode in ("geo", "vib", "occ"):
        mode_errors = [err for err in errors if err.mode == mode]
        if not mode_errors:
            continue

        err_x = {}
        err_y = {}
        err_x_inters = {}  # x-values of intersections with var(R)
        err_x_to_mark = {}
        err_legend = {}
        err_disp = {}

        for err in mode_errors:
            ats = "Atom"
            if len(ats) > 1:
                ats += "s"
            ats += " " + range_to_str([at.oriN for at in err.atoms])
            if mode == "geo":
                direction = err.disp_label
                disp = err.lin_disp  # in A
                err_legend[err] = ats + " along " + direction
            elif mode == "vib":
                disp = err.displacements[:]  # in A
                err_legend[err] = ats
            elif mode == "occ":
                disp = err.displacements[:] * 100   # in % for occ
                err_legend[err] = ats
            else:
                raise ValueError(f'Unknown mode "{mode}"')
            err_disp[err] = disp
        xvals = [d for k in err_disp for d in err_disp[k]]
        xrange = [min(xvals) - abs(max(xvals) - min(xvals)) * 0.05,
                  max(xvals) + abs(max(xvals) - min(xvals)) * 0.05]
        for err in mode_errors:
            rmin = err.get_r_min
            var = err.var_r
            disp = err_disp[err]
            interp_f = interpolate.interp1d(disp, err.rfacs, bounds_error=False)
            x = np.arange(min(xvals), max(xvals),
                          (xrange[1] - xrange[0])*1e-4)
            y = interp_f(x)
            indmark = [np.argmin(abs(x - v)) for v in disp]
            err_x_to_mark[err] = indmark
            err_x[err] = x[indmark[0]:indmark[-1]+1]
            err_y[err] = y[indmark[0]:indmark[-1]+1]
            err_x_to_mark[err] = [v - indmark[0] for v in indmark]
            err_x_inters[err] = []
            if var and any(err_y[err]) > rmin + var:
                rv = rmin + var
                err_x_inters[err] = [
                    x for i, x in enumerate(err_x[err])
                    if 0 < i and (np.sign(err_y[err][i-1]-rv)
                                  != np.sign(err_y[err][i]-rv))]
        # plot combined figure
        rmax = max(r for err in mode_errors for r in err.rfacs)
        fig = plt.figure(figsize=(5.8, 5.8))
        ax = fig.add_subplot(1, 1, 1)
        if mode != "occ":
            ax.set_xlabel('Deviation from bestfit value (Å)')
        else:
            ax.set_xlabel('Site occupation (%)')
        ax.set_ylabel('Pendry R-factor')
        ax.set_title(titles[mode])
        if var and rmin + var < rmax + (rmax-rmin)*0.1:
            ax.plot(xrange, [rmin + var]*2, color="slategray")
            inters = sorted([x for err in mode_errors
                             for x in err_x_inters[err]]
                            + [xrange[0], xrange[1]])
            (ind, diff) = max_diff(inters)
            text_x = inters[ind] - diff/2
            text_y = rmin + var + (rmax-rmin)*0.015
            va = "bottom"
            ind_at_text = {err: np.argmin([abs(x - text_x)
                                           for x in err_x[err]])
                           for err in mode_errors}
            if sum([err_y[err][ind_at_text[err]] > rmin + var
                    for err in mode_errors]) > len(mode_errors) / 2:
                text_y = rmin + var - (rmax-rmin)*0.015
                va = "top"
            ax.text(text_x, text_y, "$R_P + var(R_P)$", ha="center", va=va,
                    bbox=dict(facecolor='white', edgecolor='none',
                              alpha=0.6, pad=0.5))
        for err in mode_errors:
            ax.plot(err_x[err], err_y[err], '-o', label=err_legend[err],
                    markevery=err_x_to_mark[err])
        ax.set_xlim(*xrange)
        ax.set_ylim(rmin - ((rmax-rmin)*0.1), rmax + ((rmax-rmin)*0.1))
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
                inters = sorted(err_x_inters[err] + [xrange[0], xrange[1]])
                (ind, diff) = max_diff(inters)
                text_x = inters[ind] - diff/2
                text_y = rmin + var + (rmax-rmin)*0.015
                va = "bottom"
                ind_at_text = np.argmin([abs(x - text_x)
                                         for x in err_x[err]])
                if err_y[err][ind_at_text] > rmin + var:
                    text_y = rmin + var - (rmax-rmin)*0.015
                    va = "top"
                axs[figcount].text(text_x, text_y, "$R_P + var(R_P)$",
                                   fontsize=6, ha="center", va=va,
                                   bbox=dict(facecolor='white',
                                             edgecolor='none',
                                             alpha=0.6, pad=0.5))
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
                axs[figcount].set_xlabel('Deviation from bestfit value (Å)',
                                         fontsize=8)
            else:
                axs[figcount].set_xlabel('Site occupation (%)', fontsize=8)

            # add uncertainties in plot
            r_min = min(err.rfacs)
            p_best = err.lin_disp[err.rfacs.index(r_min)]
            error_estimates = err.get_error_estimates
            if error_estimates[0]:
                l_bound = p_best-error_estimates[0]
                draw_error(axs[figcount], l_bound, err, r_interval=(rmax-rmin))
            if error_estimates[1]:
                u_bound = p_best+error_estimates[1]
                draw_error(axs[figcount], u_bound, err, r_interval=(rmax-rmin))
            axs[figcount].set_ylabel('Pendry R-factor', fontsize=8)
            axs[figcount].legend(fontsize="x-small", frameon=False)
            # axs[figcount].set_title(err_legend[err], fontsize=8)
            figcount += 1
        for ax in axs[figcount:]:
            ax.axis('off')
        fig.tight_layout(rect=(0, 0, 1, 0.965))
        fig.suptitle(titles[mode])
        figs.append(fig)
    logger.log(1, f'Number of error figures: {len(figs)}')
    return figs


def draw_error(axis, bound, error, r_interval):
    """Adds annotation for statistical error estimates to individual
    error plots.

    Parameters
    ----------
    axis : matplotlib axis
        axis of the error plot to be annotated.
    bound : float
        Error bound == p_min \pm \sigma(p).
    error : R_Error
        Error object used for error estimation.
    r_interval : float
        R-factor range shown in the plot (rmax - rmin). Used for
        graphical scaling.
    """
    r_min = min(error.rfacs)
    p_best = error.lin_disp[error.rfacs.index(r_min)]
    var_r =error.var_r

    # put vertical line at minimum R-factor
    #axis.vlines(x=p_best, color="black", ymin=0, ymax=2, lw=0.5)
    # put dashed line at bound
    axis.vlines(x = bound, ymin = r_min - r_interval*0.024, ymax = r_min+var_r, ls=":", lw=0.5,
            color='slategrey')

    # make arrow "notch"
    axis.vlines(x=p_best,
                ymin= r_min - r_interval*0.032,
                ymax= r_min - r_interval*0.007,
                lw=0.75, color="slategrey")
    # arrow from min R to intersection x
    axis.annotate('', xy=(p_best,r_min - r_interval*0.007),
                  xytext=(bound,r_min - r_interval*0.007),
                  arrowprops=dict(arrowstyle='<-', lw=0.75,
                                  color="slategray", shrinkA=0, shrinkB=0))
    # label error under arrow
    axis.annotate(
        text=f"{format_col_content(bound-p_best)}",
        xy=((p_best+bound)/2,r_min - r_interval*0.022),
        ha='center', va='top', fontsize=4.5)

def write_errors_pdf(figs, filename="Errors.pdf"):
    """Writes a list of figures to a pdf file."""
    if not figs:
        raise ValueError("No figures to write.")
    try:
        pdf = PdfPages(filename)
        for fig in figs:
            pdf.savefig(fig)
    except PermissionError:
        logger.warning("Failed to write to file " + filename
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
