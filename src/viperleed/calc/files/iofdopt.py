"""Functions for writing output from full-dynamic optimization."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2021-10-25'
__license__ = 'GPLv3+'

import copy
import logging

import numpy as np
from numpy.polynomial import Polynomial

from viperleed.calc.files.iorfactor import read_rfactor_columns
from viperleed.calc.files.ivplot import plot_iv
from viperleed.calc.lib.matplotlib_utils import CAN_PLOT
from viperleed.calc.lib.matplotlib_utils import close_figures
from viperleed.calc.lib.matplotlib_utils import log_without_matplotlib
from viperleed.calc.lib.matplotlib_utils import prepare_matplotlib_for_calc

if CAN_PLOT:
    prepare_matplotlib_for_calc()
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    from matplotlib import cm


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # Mute matplotlib debug messages                 # TODO: perhaps nicer to use at_level only in the relevant spots? See also iorfactor and ivplot


def write_fd_opt_csv(points, which, filename="FD_Optimization.csv", sep=","):
    """
    Writes results from the full dynamic optimization into a CSV file

    Parameters
    ----------
    points : numpy.array
        Two columns, containing points and corresponding R-factors
    which : str
        The parameter that was being optimized, for column title
    filename : str, optional
        Name of the csv file to write. The default is "FD_Optimization.csv".
    sep : str, optional
        The separator to use. The default is ",".

    Returns
    -------
    None.

    """

    sorted_points = points[points[:, 0].argsort()]
    title = which
    if which in ["a", "b", "c", "ab", "abc"]:
        title += " scaling"
    width = max(len(title), 12)

    output = title.rjust(width) + sep + "R".rjust(width) + "\n"
    for xy in sorted_points:
        output += "{:.5f}".format(xy[0]).rjust(width) + sep
        output += "{:.5f}".format(xy[1]).rjust(width) + "\n"

    try:
        with open(filename, "w") as wf:
            wf.write(output)
    except Exception as e:
        logger.warning("Failed to write "+filename + ": " + str(e))
    return


@log_without_matplotlib(logger, msg='Skipping error plotting.')
def write_fd_opt_pdf(points, which, filename="FD_Optimization.pdf",
                     parabola=None):
    """
    Plots results from the full dynamic optimization into a pdf file

    Parameters
    ----------
    points : numpy.array
        Two columns, containing points and corresponding R-factors
    which : str
        The parameter that was being optimized, for x-axis label
    filename : str, optional
        File name to write to. The default is "FD_Optimization.pdf".
    parabola : Polynomial, optional
        Pre-existing fit parabola. If not passed, data will be fit again here.

    Returns
    -------
    None.

    """
    title = which
    if which in ["a", "b", "c", "ab", "abc"]:
        title += " scaling"

    fig = plt.figure(figsize=(5.8, 4.1))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel(title)
    ax.set_ylabel('Pendry R-factor')

    if len(points) > 2:
        if parabola is None:
            parabola = Polynomial.fit(points[:, 0], points[:, 1], 2)
        coefs = parabola.convert(domain=[-1, 1]).coef
        p_min = -0.5*coefs[1] / coefs[2]
        p_r = parabola(p_min)
        xx, yy = parabola.linspace()
        ax.plot(xx, yy, lw=2)
        xlims = (min(points[:, 0]), max(points[:, 0]))   # data range x
        ylims = (min(points[:, 1]), max(points[:, 1]))   # data range y
        namePos = ((xlims[1] + xlims[0])*0.5,
                   ylims[1] - 0.1*(ylims[1] - ylims[0]))
        if coefs[2] > 0:
            ax.annotate("Minimum at {:.4f}\nR = {:.4f}".format(p_min, p_r),
                        namePos, fontsize=10, ha="center")
    ax.plot(points[:, 0], points[:, 1], 'o', c='darkslategray')
    fig.tight_layout()

    # write
    try:
        pdf = PdfPages(filename)
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
    close_figures(plt, fig)


def write_fd_opt_beams_pdf(rp, points, which, tmpdirs, best_rfactors,
                           filename="FD_Optimization_beams.pdf"):
    """
    Plots all I(V) curves that were calculated during the optimization.

    Parameters
    ----------
    rp : Rparams
        Run parameters.
    points : numpy.array
        Two columns, containing points and corresponding R-factors.
    which : str
        The parameter that was being optimized, for legend.
    tmpdirs : list of str
        Paths of directories from which I(V) data should be read.
    best_rfactors : list of float
        The R-factors from the best calculation. Will be written into plots.
    filename : str, optional
        File name to write to. The default is "FD_Optimization_beams.pdf".

    Returns
    -------
    None.

    """

    if not CAN_PLOT:
        logger.debug("Necessary modules for plotting not found. Skipping "
                     "error plotting.")
        return

    new_order = points[:, 0].argsort()
    tmpdirs = list(np.array(tmpdirs)[new_order])
    points = points[new_order]
    best_point_ind = np.argmin(points[:, 1])

    theodata = []
    exp_to_use = None

    for i, path in enumerate(tmpdirs):
        try:
            theospec, expspec = read_rfactor_columns(path)
        except Exception:
            logger.warning("Error reading calculated spectra for point {:.4f}"
                           .format(points[i, 0]))
            continue
        theodata.append(theospec)
        if i == best_point_ind:
            exp_to_use = expspec
    if exp_to_use is None:
        logger.warning("Failed to read in experimental spectra from best "
                       "full-dynamic run. Output of collected I(V) spectra "
                       "will be skipped.")
        return
    legends = ["{} = {:.4f}".format(which, points[i, 0])
               for i in range(len(points))] + ["Experimental"]
    annotations = []
    if best_rfactors:
        annotations = ["R = {:.4f}".format(r) for r in best_rfactors]
    if rp.PLOT_IV["overbar"]:
        labelstyle = "overbar"
    else:
        labelstyle = "minus"
    labelwidth = max([beam.getLabel(style=labelstyle)[1]
                      for beam in rp.expbeams])
    labels = [b.getLabel(lwidth=labelwidth, style=labelstyle)[0]
              for b in rp.expbeams]
    formatting = copy.deepcopy(rp.PLOT_IV)  # use to set colors
    formatting['colors'] = (
        list(cm.get_cmap('viridis', len(points)).colors)
        + [np.array([0, 0, 0, 1])])
    formatting['curve_line_widths'] = [0.5] * len(points) + [1.]
    formatting['curve_line_widths'][best_point_ind] = 1.
    try:
        plot_iv(theodata + [exp_to_use], filename,
                labels=labels, annotations=annotations, legends=legends,
                formatting=formatting)
    except Exception:
        logger.warning("Error plotting collected I(V) curves.",
                       exc_info=True)
    return
