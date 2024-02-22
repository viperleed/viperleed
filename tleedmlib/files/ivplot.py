# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 18:47:39 2021

@author: fkrau
"""

import logging

import numpy as np

try:
    import matplotlib
except Exception:
    _CAN_PLOT = False
else:
    _CAN_PLOT = True
    matplotlib.rcParams.update({'figure.max_open_warning': 0})
    matplotlib.use('Agg')  # !!! check with Michele if this causes conflicts
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    import matplotlib.ticker as plticker
    plt.style.use('viperleed.tleedm')

from viperleed.tleedmlib.classes.beam import Beam

logger = logging.getLogger("tleedm.files.ivplot")


def plot_iv(data, filename, labels=[], annotations=[],
            legends=[], formatting=None):
    '''
    Creates a single PDF file containing plots of I(V) curves

    Parameters
    ----------
    data : list of datasets, or single dataset
        dataset : list of Beam or of np.ndarray objects
            contains I(V) data per beam
        All datasets must contain the same beams in the same order, but beams
        can be empty
    filename : str or None
        Name of the file (with or without extension) to which the plots
        will be saved. If None, do not write to file, but return a list
        of matplotlib figure objects instead.
    labels : kwarg, list of str
        Labels for the individual beams, e.g. (h|k) as strings. Must be same
        length as number of beams. If labels are not passed and at least one
        dataset is of type Beam, labels will be generated from that dataset.
    annotations : kwarg, list of str
        Additional information to print for each beam, for example the
        R-factor. Must be same length as the number of beams.
    legends : kwarg, list of str
        How to label the different datasets, in order. Must be same length as
        the number of datasets.
    formatting : kwarg, dict
        dict containing formatting instructions:
        axes : str
            Which axes to draw.
            all:  draw all axes (left, right, top, bottom) for all panels
            lb:   draw only left and bottom axes
            b:    draw only bottom axes
            none: draw no axes except at bottom of last panel in each column
        colors : tuple (str, str)
            Define alternative colors for drawing I(V) curves.
            default: None
        legend : str
            Which legends to print.
            all:   print legend for each panel
            first: print legend only for the first panel on each page
            tr:    print legend only for the top-right panel on each page
            none:  do not print any legend.
        perpage : int or (int, int)
            Define how many figures to plot on each page. Either tuple
            (columns, rows), or single integer (preferred). For single int,
            layout will be adapted automatically. Numbers that are not nicely
            divisible may be rounded up, resulting in some whitespace.
            default: 2 (one column, two rows).

    Returns
    -------
    None

    '''
    global _CAN_PLOT
    if not _CAN_PLOT:
        logger.debug("Necessary modules for plotting not found. Skipping "
                     "R-factor plotting.")
        return

    # check data
    if type(data) not in (list, tuple):
        raise TypeError("Expected data as a list or tuple, found "
                        + str(type(data)))
    if type(data[0]) not in (list, tuple):
        data = [data]       # assume single set of beams
    n_beams = len(data[0])
    for dataset in data:
        if type(dataset) not in (list, tuple):
            raise TypeError(
                "Expected data as a list or tuple, found "
                + str(type(dataset[0])))
        if type(dataset[0]) not in (Beam, np.ndarray):
            raise TypeError(
                "Expected data as a list of tleedmlib Beams or numpy arrays, "
                "found " + str(type(dataset[0])))
        elif len(dataset) != n_beams:
            raise ValueError("Different beam sets or not of equivalent "
                             "length.")

    if labels and len(labels) != n_beams:
        raise ValueError("Number of labels does not match number of beams.")
    if annotations and len(annotations) != n_beams:
        raise ValueError("Number of annotations does not match number of "
                         "beams.")

    # set formatting parameters
    figs_per_page = 2
    plotcolors = []
    linewidths = [1.5] * n_beams
    print_legend = 'all'
    print_axes = 'all'
    if formatting is not None:
        if 'axes' in formatting:
            print_axes = formatting['axes']
        if 'colors' in formatting:
            plotcolors = formatting['colors']
        if 'legend' in formatting:
            print_legend = formatting['legend']
        if 'perpage' in formatting:
            figs_per_page = formatting['perpage']
        if 'linewidths' in formatting:
            linewidths = formatting['linewidths']

    # read data
    readlabels = False
    if not labels:
        readlabels = True   # fill labels from first dataset
    xy_per_beam_per_dataset = []
    for dataset in data:
        if type(dataset[0]) == np.ndarray:
            xy_per_beam = [np.copy(xy) for xy in dataset]
        else:    # list of beams
            xy_per_beam = []
            if readlabels:
                if formatting and formatting["overbar"]:
                    labelstyle = "overbar"
                else:
                    labelstyle = "minus"
                labelwidth = max([beam.getLabel(style=labelstyle)[1]
                                  for beam in dataset])
                labels = [beam.getLabel(lwidth=labelwidth, style=labelstyle)[0]
                          for beam in dataset]
                readlabels = False
            for b in dataset:
                xy = np.array([[e, b.intens[e]]
                               for e in sorted(b.intens.keys())])
                xy_per_beam.append(xy)
        xy_per_beam_per_dataset.append(xy_per_beam)

    # find min and max values of x and y for plotting all curves
    # on the same horizontal scale and leaving a little y space for the legend
    all_xy = [xy for xy_per_beam in xy_per_beam_per_dataset
              for xy in xy_per_beam if len(xy) != 0]
    xmin = min(min(xy[:, 0]) for xy in all_xy)
    xmax = max(max(xy[:, 0]) for xy in all_xy)

    ymin = min(min(xy[:, 1]) for xy in all_xy)
    ymax = max(max(xy[:, 1]) for xy in all_xy)
    dy = ymax - ymin

    # Figure out how to arrange the figures
    if type(figs_per_page) == int:
        xfigs = int(np.round(np.sqrt(figs_per_page/2)))
        yfigs = int(np.ceil(figs_per_page / xfigs))
        if xfigs*yfigs != figs_per_page:
            logger.debug("R-factor plots: Figures per page corrected from {} "
                         "to {}".format(figs_per_page, xfigs*yfigs))
            figs_per_page = xfigs*yfigs
    else:
        xfigs, yfigs = figs_per_page
        figs_per_page = xfigs*yfigs
    figsize = (7., 3.5 * yfigs / xfigs)

    # set ticks spacing to 50 eV and round the x limits to a multiple of it
    tick = 50
    oritick = plticker.MultipleLocator(base=tick)
    majortick = oritick
    minortick = None
    xlims = (np.floor(xmin/tick)*tick,
             np.ceil(xmax/tick)*tick)
    dx = xlims[1] - xlims[0]

    # positions and scales, depending on number of figures per page
    figs = []
    fontscale = 1 / np.sqrt(xfigs)
    linewidth = 1.5 * fontscale
    ylims = (ymin - 0.02*dy, ymax + 0.22*dy/fontscale)
    namePos = (xlims[0] + 0.02*dx, ylims[1] - 0.1*dy/fontscale)  # eg labels
    annotationPos = (namePos[0], namePos[1]-0.085*dy/fontscale)  # eg R-factor

    if dx / (tick * fontscale) > 16:  # too many labelled ticks
        minortick = majortick
        newbase = int(np.ceil(dx / (tick * fontscale * 16))) * tick
        majortick = plticker.MultipleLocator(base=newbase)
    ticklen = 3*fontscale

    # axes helper variables
    axes_visible = {'left': True, 'right': True, 'bottom': True, 'top': True}
    if print_axes != 'all':
        axes_visible['top'] = False
        axes_visible['right'] = False
        if 'l' not in print_axes:
            axes_visible['left'] = False
        if print_axes == 'none':
            axes_visible['bottom'] = False

    # the following will spam the logger with debug messages; disable.
    loglevel = logger.level
    logger.setLevel(logging.INFO)
    try:
        fig_exists = False
        fig_index_on_page = 0
        for ct in range(n_beams):   # iterate through beams
            if all([len(xy_per_beam_per_dataset[i][ct]) == 0
                    for i in range(len(data))]):
                continue   # no data for this beam in any dataset, skip
            if (fig_index_on_page >= figs_per_page) or (not fig_exists):
                # need a new figure
                fig_exists = True # at least one fig exists
                fig_index_on_page = 0
                fig, axs = plt.subplots(yfigs, xfigs, figsize=figsize,
                                        squeeze=False)
                axs = axs.flatten(order='C')  # flatten row-style
                fig.subplots_adjust(left=(0.03 / (xfigs * fontscale)),
                                    right=(1 - 0.03 / (xfigs * fontscale)),
                                    bottom=(0.14 / (yfigs * fontscale)),
                                    top=(1 - 0.02 / (yfigs * fontscale)),
                                    wspace=0.04 / fontscale,
                                    hspace=0.02 / fontscale)
                figs.append(fig)
                [ax.set_xlim(*xlims) for ax in axs]
                [ax.set_ylim(*ylims) for ax in axs]
                [ax.get_yaxis().set_visible(False) for ax in axs]
                [sp.set_linewidth(0.7*linewidth) for ax in axs
                 for sp in ax.spines.values()]
                [ax.xaxis.set_major_locator(majortick) for ax in axs]
                if minortick is not None:
                    [ax.xaxis.set_minor_locator(minortick) for ax in axs]
                for i, ax in enumerate(axs):
                    for k in axes_visible:
                        ax.spines[k].set_visible(axes_visible[k])
                    if (((i//xfigs) + 1 == figs_per_page//xfigs)
                            or n_beams - (ct + i) <= xfigs):
                        ax.set_xlabel("Energy (eV)", fontsize=10*fontscale,
                                      labelpad=4*fontscale)
                        ax.tick_params(
                            which='both', bottom=True,
                            top=axes_visible['top'], axis='x',
                            direction='in', labelsize=9*fontscale,
                            pad=5*fontscale, width=0.7*linewidth,
                            length=ticklen)
                        ax.spines['bottom'].set_visible(True)
                    else:
                        ax.tick_params(
                            which='both', bottom=axes_visible['bottom'],
                            top=axes_visible['top'], axis='x', direction='in',
                            labelbottom=False, width=0.7*linewidth,
                            length=ticklen)
                    if minortick is not None:
                        ax.tick_params(which='minor', length=ticklen*0.5)
            if plotcolors:
                if not all([matplotlib.colors.is_color_like(s)
                            for s in plotcolors]):
                    plotcolors = []
                    logger.warning("plot_iv: Specified colors not "
                                   "recognized, reverting to default colors")
            for i in range(len(data)):
                if legends:
                    label = legends[i]
                else:
                    label = 'Beamset {}'.format(i+1)
                xy = xy_per_beam_per_dataset[i][ct]
                if i < len(linewidths):
                    lw = linewidths[i] * fontscale
                else:
                    lw = linewidth
                if i < len(plotcolors):
                    axs[fig_index_on_page].plot(xy[:, 0], xy[:, 1], label=label,
                                  linewidth=lw,
                                  color=plotcolors[i])
                else:
                    axs[fig_index_on_page].plot(xy[:, 0], xy[:, 1], label=label,
                                  linewidth=lw)
            if labels:
                axs[fig_index_on_page].annotate(labels[ct], namePos, fontsize=10*fontscale)
            if annotations:
                axs[fig_index_on_page].annotate(annotations[ct], annotationPos,
                                  fontsize=10*fontscale)
            if ((print_legend == 'all'
                    or (print_legend == 'first' and fig_index_on_page == 0)
                    or (print_legend == 'tr'
                        and (fig_index_on_page//xfigs == 0 and ((fig_index_on_page+1) % xfigs == 0
                                                  or ct + 1 == n_beams))))
                    and len(data) > 1):
                legendscale = 1.
                if len(data) > 2:
                    legendscale = 1/np.sqrt(len(data)-1)
                legend = axs[fig_index_on_page].legend(fontsize=9*fontscale*legendscale,
                                         loc="upper right", frameon=False,
                                         ncol=(len(data) // 3 + 1))
                legend.get_frame().set_linewidth(linewidth)
            fig_index_on_page += 1

        # finally, in case the last figure is empty (i.e. the number of beams
        # is odd) turn off the last axes (but leave the blank space).
        if n_beams % figs_per_page != 0:
            for a in axs[-(figs_per_page - n_beams % figs_per_page):]:
                a.axis('off')
    except Exception:
        logger.error("plot_iv: Error while compiling figures.", exc_info=True)
        logger.setLevel(loglevel)

    # if a filename is given write to PDF, else return list of figs
    if filename is not None:
        try:
            pdf = PdfPages(filename)
        except PermissionError:
            logger.error(f"writeRfactorPdf: Cannot open file {filename}. "
                         "Aborting.")
            return

        try:
            for fig in figs:
                pdf.savefig(fig)
                plt.close(fig)
        except Exception:
            logger.error("plot_iv: Error while writing {}: ".format(filename),
                        exc_info=True)
        finally:
            pdf.close()
            return
    else:
        return figs

