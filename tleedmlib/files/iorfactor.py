# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 12:17:56 2020

@author: Florian Kraushofer

Functions for reading and writing files relevant to the rfactor calculation
"""

import logging
import re
import numpy as np
from viperleed import fortranformat as ff
import os

from viperleed.tleedmlib.files.beams import writeAUXEXPBEAMS
from viperleed.tleedmlib.leedbase import getBeamCorrespondence, getYfunc

try:
    import matplotlib
    matplotlib.rcParams.update({'figure.max_open_warning': 0})
    matplotlib.use('Agg')  # !!! check with Michele if this causes conflicts
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    import matplotlib.ticker as plticker
except Exception:
    plotting = False
else:
    plotting = True

logger = logging.getLogger("tleedm.files.iorfactor")


def readROUT(filename="ROUT"):
    """
    Reads the ROUT file.

    Parameters
    ----------
    filename : string, optional
        Pass if you want to read from a file other than 'ROUT'

    Returns
    -------
    tuple (r, r_int, r_frac), float v0rshift, list rfaclist
        r, r_int, r_frac: tuple of floats:
            average R-factor for all beams, integer-order beams and
            fractional order beams.
        v0rshift: inner potential shift corresponding to r, r_int, r_frac
        rfaclist: list of floats, r-factors per beam in order of experimental
            beams

    """
    try:
        with open(filename, 'r') as rf:
            lines = rf.readlines()
    except Exception:
        logger.error("Could not open " + filename + " file")
        raise
    line = ""
    i = 0
    while "AVERAGE R-FACTOR =" not in line and i+1 < len(lines):
        i += 1
        line = lines[-i]
    if line == "":
        return 0, 0, []
    rfac = 0
    rfac_int = -1
    rfac_frac = -1
    v0rshift = 0
    rgx = re.compile(r'.+SHIFT\s*(?P<shift>[-0-9.]+).+=\s*(?P<rfac>[-0-9.]+)')
    m = rgx.match(line)
    if m:
        try:
            rfac = float(m.group("rfac"))
        except ValueError:
            logger.error("Could not read R-factor from "+filename)
        try:
            v0rshift = float(m.group("shift"))
        except ValueError:
            logger.error("Could not read inner potential shift from "
                         + filename)
    else:
        return (0, 0, 0), 0, []
    # now read the R-factors per beam at v0rshift
    rfaclist = []
    for line in [li for li in lines if len(li) >= 70]:
        values = line[18:69].split()
        # limits to 999 beams; use [17:69] if more are required
        try:
            index = int(values[0])
            v0r = float(values[2])
            rav = float(values[-1])
        except (ValueError, IndexError):
            pass    # ignore line
        else:
            if v0r == v0rshift and index > 0:
                if index != len(rfaclist)+1:
                    logger.warning("Unexpected index mismatch in readROUT. "
                                   "Reading R-factors per beam will fail.")
                else:
                    rfaclist.append(rav)
            elif v0r == v0rshift and index == -1:
                if line.startswith("AV.-INT"):
                    rfac_int = rav
                elif line.startswith("AV.-FRAC"):
                    rfac_frac = rav
    return (rfac, rfac_int, rfac_frac), v0rshift, rfaclist


def readROUTSHORT(filename="ROUTSHORT"):
    """
    Reads the ROUTSHORT file. This is very minimalist and just contains one
    average R-factor per line.

    Parameters
    ----------
    filename : string, optional
        Pass if you want to read from a file other than 'ROUTSHORT'

    Returns
    -------
    rfaclist : list of floats
        r-factors from the ROUTSHORT file

    """
    rfaclist = []
    try:
        with open(filename, 'r') as rf:
            lines = rf.readlines()
    except Exception:
        logger.error("Could not open " + filename + " file")
        raise
    for line in lines:
        line = line.strip()
        if not line:
            continue
        try:
            rfaclist.append(float(line))
        except ValueError:
            logger.warning("Unexpected value in " + filename + " file, could "
                           "not transform to float: "+line)
    return rfaclist


def writeWEXPEL(sl, rp, theobeams, filename="WEXPEL", for_error=False):
    """
    Writes input file WEXPEL for R-factor calculation.

    Parameters
    ----------
    sl : Slab
        The Slab object containing atom information.
    rp : Rparams
        The run parameters.
    theobeams : list of Beam
        The theoretical beams, containing I(V) data.
    filename : str, optional
        Name of the file that will be written. The default is "WEXPEL".

    Returns
    -------
    None.

    """
    theoEnergies = []
    for b in theobeams:
        theoEnergies.extend([k for k in b.intens if k not in theoEnergies])
    theoEnergies.sort()
    expEnergies = []
    for b in rp.expbeams:
        expEnergies.extend([k for k in b.intens if k not in expEnergies])
    expEnergies.sort()
    minen = max(min(expEnergies), rp.THEO_ENERGIES[0])
    maxen = min(max(expEnergies), rp.THEO_ENERGIES[1])
    if not for_error:
        real_iv_shift = rp.IV_SHIFT_RANGE[:2]
    else:
        real_iv_shift = [rp.best_v0r] * 2
    # extend energy range if they are close together
    if abs(min(expEnergies) - rp.THEO_ENERGIES[0]) < abs(real_iv_shift[0]):
        minen = (max(min(expEnergies), rp.THEO_ENERGIES[0])
                 - real_iv_shift[0])
    if abs(max(expEnergies) - rp.THEO_ENERGIES[1]) < abs(real_iv_shift[1]):
        maxen = (min(max(expEnergies), rp.THEO_ENERGIES[1])
                 + real_iv_shift[1]) + 0.01
    step = min(expEnergies[1]-expEnergies[0], theoEnergies[1]-theoEnergies[0])
    if rp.IV_SHIFT_RANGE[2] > 0:
        vincr = rp.IV_SHIFT_RANGE[2]
        # step = min(step, vincr)
    else:
        vincr = step
    # find correspondence experimental to theoretical beams:
    beamcorr = getBeamCorrespondence(sl, rp)
    # integer & fractional beams
    iorf = []
    for (i, beam) in enumerate(rp.expbeams):
        if beam.hk[0] % 1.0 != 0.0 or beam.hk[1] % 1.0 != 0.0:
            iorf.append(2)
        else:
            iorf.append(1)
    iorf.extend([0]*(len(rp.ivbeams)-len(rp.expbeams)))

    f72 = ff.FortranRecordWriter("F7.2")
    i3x25 = ff.FortranRecordWriter("25I3")
    i3 = ff.FortranRecordWriter("I3")
    output = " &NL1\n"
    output += (" EMIN=" + f72.write([minen]).rjust(9) + ",\n")
    output += (" EMAX=" + f72.write([maxen]).rjust(9) + ",\n")
    output += (" EINCR=" + f72.write([step]).rjust(8) + ",\n")
    output += " LIMFIL=      1,\n"  # number of consecutive input files
    output += " IPR=         0,\n"  # output formatting
    output += (" VI=" + f72.write([rp.V0_IMAG]).rjust(11) + ",\n")
    output += " V0RR=      0.0,\n"
    output += (" V01=" + f72.write([real_iv_shift[0]]).rjust(10) + ",\n")
    output += (" V02=" + f72.write([real_iv_shift[1]]).rjust(10) + ",\n")
    output += (" VINCR=" + f72.write([vincr]).rjust(8) + ",\n")
    output += " ISMOTH=" + i3.write([rp.R_FACTOR_SMOOTH]).rjust(7) + ",\n"
    output += " EOT=         0,\n"
    output += " PLOT=        1,\n"
    output += " GAP=         0,\n"
    output += " &END\n"
    output += (i3x25.write([n+1 for n in beamcorr]) + "\n")
    if len(beamcorr) % 25 == 0:
        output += "\n"
    for i in range(0, 2):  # redundant since indices are already taken care of
        output += i3x25.write([n+1 for n in range(0, len(rp.expbeams))]) + "\n"
        if len(rp.expbeams) % 25 == 0:
            output += "\n"
    output += i3x25.write(iorf) + "\n"
    if len(iorf) % 25 == 0:
        output += "\n"
    output += "&NL2\n"
    output += " NSSK=    0,\n"
    if rp.R_FACTOR_TYPE == 1:
        output += " WR=      0.,0.,1.,\n"  # Pendry
    elif rp.R_FACTOR_TYPE == 2:
        output += " WR=      1.,0.,0.,\n"  # R2
    else:
        output += " WR=      0.,1.,0.,\n"  # Zanazzi-Jona
    output += """ &END
 &NL3
 NORM=           1,
 INTMAX=    999.99,
 PLSIZE=   1.0,1.0,
 XTICS=         50,
 &END
 """
    auxexpbeams = writeAUXEXPBEAMS(rp.expbeams, header=rp.systemName,
                                   write=True, numbers=False)
    output += auxexpbeams
    output += "\n"
    # information about gaps in the experimental spectra would go here
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except Exception:
        logger.error("Failed to write "+filename)
        raise
    logger.debug("Wrote to R-factor input file "+filename+" successfully")
    return


def writeRfactPARAM(rp, theobeams, for_error=False, only_vary=None):
    """
    Generates the PARAM file for the rfactor calculation.

    Parameters
    ----------
    rp : Rparams
        The run parameters.
    theobeams : list of Beam
        The theoretical beams, containing I(V) data.

    Returns
    -------
    None.

    """
    theoEnergies = []
    for b in theobeams:
        theoEnergies.extend([k for k in b.intens if k not in theoEnergies])
    theoEnergies.sort()
    expEnergies = []
    for b in rp.expbeams:
        expEnergies.extend([k for k in b.intens if k not in expEnergies])
    expEnergies.sort()
    minen = min(min(expEnergies), min(theoEnergies))
    maxen = max(max(expEnergies), max(theoEnergies))
    if rp.IV_SHIFT_RANGE[2] > 0:
        step = rp.IV_SHIFT_RANGE[2]
    else:
        step = min(expEnergies[1]-expEnergies[0], rp.THEO_ENERGIES[2])
    ngrid = int(np.ceil(((maxen-minen)/step)*1.1))
    n_var = 1
    if for_error:
        if not only_vary:
            logger.warning("Rfactor PARAM for error: Parameters under "
                           "variation not passed.")
            only_vary = [sp for sp in rp.searchpars
                         if sp.atom in rp.search_atlist]
        n_var = max([sp.steps for sp in only_vary])
    output = """
C  MNBED  : number of beams in experimental spectra before averaging
C  MNBTD  : number of beams in theoretical spectra before averaging

      PARAMETER (MNBED = {}, MNBTD = {})""".format(len(rp.expbeams),
                                                   len(theobeams))
    output += """

C  MNET   : number of energies in theoretical beam at time of reading in
C  MNGP  : greater equal number of grid points in energy working grid (ie after
C           interpolation)
C  MNS    : number of geometries including those, that are skipped

      PARAMETER (MNET = {}, MNGP = {})
      PARAMETER (MNS = {})""".format(len(theoEnergies), ngrid, n_var)
    output += """

C  MNGAP  : number of gaps in the experimental spectra (NOTE: if there are no
C           gaps in the spectra set MNGAP to 1 to avoid zero-sized arrays)

      PARAMETER (MNGAP = 1)
"""
    # write PARAM
    try:
        with open("PARAM", "w") as wf:
            wf.write(output)
    except Exception:
        logger.error("Failed at writing PARAM file for R-factor calculation.")
        raise
    return


def writeRfactorPdf(beams, colsDir='', outName='Rfactor_plots.pdf',
                    analysisFile='', v0i=0., formatting=None):
    '''
    Creates a single PDF file containing the plots of R-factors

    Parameters
    ----------
    beams : list of (name, R) tuples
           name: str
                 formatted fractional index of the beam
           R: float
              R-factor value
    colsDir : kwarg, str
        path to folder containing the files theo.column and exp.column
        generated by the R-factor routine of TensErLEED.
        default: current path
    outName : kwarg, str
        name of the file (with or without extension) to which the plots
        will be saved.
        default: 'Rfactor_plots.pdf'
    analysisFile : kwarg, string
        if not empty, a more extensive R-factor analysis pdf with
        calculated Y-functions and absolute errors will be written to the
        given file name.
    v0i : kwarg, float
        imaginary part of the inner potential for calculating Y-functions.
        Should always be passed if analysisFile is passed.
    formatting : kwarg, dict
        dict containing formatting instructions:
        axes : str
            Which axes to draw.
            all:  draw all axes (left, right, top, bottom) for all panels
            lb:   draw only left and bottom axes
            b:    draw only bottom axes
            none: draw no axes except at bottom of last panel in each column
        colors : tuple (str, str)
            Define alternative colors for experimental and theoretical beams.
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
    global plotting
    if not plotting:
        logger.debug("Necessary modules for plotting not found. Skipping "
                     "R-factor plotting.")
        return

    # set formatting parameters
    figs_per_page = 2
    plotcolors = None
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

    fnames = ['theo.column', 'exp.column']

    if not hasattr(beams, '__len__'):
        logger.error("writeRfactorPdf: First argument should be list, not "
                     + str(type(beams)))
        return

    xxyy = []
    for fname in fnames:
        try:
            f = open(os.path.join(colsDir, fname), 'r')
        except FileNotFoundError:
            logger.error("writeRfactorPdf: File {} not found. Aborting."
                         .format(fname))
            return
        except PermissionError:
            logger.error("writeRfactorPdf: Cannot open file {}. Aborting."
                         .format(fname))
            return

        cols = [[float(col) for col in line.split()] for line in f]

        if(np.shape(cols)[1] != 2*len(beams)):
            logger.error("writeRfactorPdf: Number of beams in file {} does "
                         "not match the input. Aborting.".format(fname))
            return

        cols = np.array(cols)
        xy = np.split(cols, len(beams), axis=1)
        # xy is now a list of 2D arrays.
        # Each array has the form [[en1, intens1], ...]
        #
        # for each beam, get rid of the points that have (en, intens) = (0, 0)
        # so that they don't screw up the plots later
        xy = [coords[~np.all(coords < 1e-3, axis=1)] for coords in xy]
        if xy:
            xxyy.append(xy)

    xyTheo = xxyy[0]
    xyExp = xxyy[1]

    # find min and max values of x and y for plotting all curves
    # on the same horizontal scale and leaving a little y space for the legend
    xmin = min(min(xy[:, 0]) for xy in [*xyTheo, *xyExp] if len(xy) != 0)
    xmax = max(max(xy[:, 0]) for xy in [*xyTheo, *xyExp] if len(xy) != 0)

    ymin = min(min(xy[:, 1]) for xy in [*xyTheo, *xyExp] if len(xy) != 0)
    ymax = max(max(xy[:, 1]) for xy in [*xyTheo, *xyExp] if len(xy) != 0)
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

    # Set up stuff needed for the plots

    try:
        pdf = PdfPages(outName)
    except PermissionError:
        logger.error("writeRfactorPdf: Cannot open file {}. Aborting."
                     .format(outName))
        return

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
    namePos = (xlims[0] + 0.02*dx, ylims[1] - 0.1*dy/fontscale)
    rPos = (namePos[0], namePos[1]-0.085*dy/fontscale)

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
        for ct, (name, rfact, theo, exp) in enumerate(zip(*zip(*beams),
                                                          xyTheo, xyExp)):
            if len(exp) == 0:
                continue
            if ct % figs_per_page == 0:
                # need a new figure
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
                            or len(beams) - (ct + i) <= xfigs):
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
            idx = ct % figs_per_page
            if plotcolors is not None:
                if not all([matplotlib.colors.is_color_like(s)
                            for s in plotcolors]):
                    plotcolors = None
                    logger.warning("writeRfactorPdf: Specified colors not "
                                   "recognized, reverting to default colors")
            if plotcolors is None:
                axs[idx].plot(theo[:, 0], theo[:, 1], label='Theoretical',
                              linewidth=linewidth)
                axs[idx].plot(exp[:, 0], exp[:, 1], label='Experimental',
                              linewidth=linewidth)
            else:
                axs[idx].plot(theo[:, 0], theo[:, 1], label='Theoretical',
                              color=plotcolors[0], linewidth=linewidth)
                axs[idx].plot(exp[:, 0], exp[:, 1], label='Experimental',
                              color=plotcolors[1], linewidth=linewidth)
            axs[idx].annotate(name, namePos, fontsize=10*fontscale)
            axs[idx].annotate("R = {:.4f}".format(rfact), rPos,
                              fontsize=10*fontscale)
            if (print_legend == 'all'
                    or (print_legend == 'first' and idx == 0)
                    or (print_legend == 'tr'
                        and (idx//xfigs == 0 and ((idx+1) % xfigs == 0
                                                  or ct + 1 == len(beams))))):
                legend = axs[idx].legend(fontsize=9*fontscale,
                                         loc="upper right", frameon=False)
                legend.get_frame().set_linewidth(linewidth)

        # finally, in case the last figure is empty (i.e. the number of beams
        # is odd) turn off the last axes (but leave the blank space).
        if len(beams) % figs_per_page != 0:
            for a in axs[-(figs_per_page - len(beams) % figs_per_page):]:
                a.axis('off')

        for fig in figs:
            pdf.savefig(fig)
            plt.close(fig)
    except Exception:
        logger.error("writeRfactorPdf: Error while writing rfactor pdf: ",
                     exc_info=True)
    finally:
        pdf.close()
        logger.setLevel(loglevel)

    if not analysisFile:
        return

    # write R-factor analysis
    try:
        pdf = PdfPages(analysisFile)
    except PermissionError:
        logger.error("writeRfactorPdf: Cannot open file {}. Aborting."
                     .format(analysisFile))
        return

    figsize = (5.8, 8.3)
    figs = []

    ylims = (ymin - 0.02*dy, ymax + 0.22*dy)
    namePos = (xlims[0] + 0.45*dx, ylims[1] - 0.1*dy)
    rPos = (namePos[0], namePos[1]-0.085*dy)
    # the following will spam the logger with debug messages; disable.
    loglevel = logger.level
    logger.setLevel(logging.INFO)
    try:
        for i, (name, rfact, theo, exp) in enumerate(zip(*zip(*beams),
                                                         xyTheo, xyExp)):
            if len(exp) == 0:
                continue
            fig, axs = plt.subplots(3, figsize=figsize,
                                    squeeze=True)
            fig.subplots_adjust(left=0.06, right=0.94,
                                bottom=0.07, top=0.98,
                                wspace=0, hspace=0.08)
            figs.append(fig)
            [ax.set_xlim(*xlims) for ax in axs]
            axs[0].set_ylim(*ylims)
            [ax.get_yaxis().set_ticks([]) for ax in axs]
            [ax.tick_params(bottom=True,
                            top=True,
                            axis='x', direction='in') for ax in axs]
            [ax.xaxis.set_major_locator(oritick) for ax in axs]
            axs[0].set_ylabel("Intensity (arb. units)")
            axs[1].set_ylabel("Y")
            axs[2].set_ylabel("\u2211(\u0394Y)\u00b2")    # sum delta Y ^2
            axs[2].set_xlabel("Energy (eV)")

            ytheo = getYfunc(theo, v0i)
            yexp = getYfunc(exp, v0i)
            dy = np.array([(ytheo[j, 0], yexp[j, 1] - ytheo[j, 1])
                           for j in range(0, min(len(ytheo), len(yexp)))])
            dysq = np.copy(dy)
            dysq[:, 1] = dysq[:, 1]**2
            idysq = np.array([dysq[0]])
            for j in range(1, len(dysq)):
                idysq = np.append(idysq, [[dysq[j, 0],
                                           idysq[j-1, 1]+dysq[j, 1]]], axis=0)

            axs[1].plot(xlims, [0., 0.], color='grey', alpha=0.2)
            if plotcolors is not None:
                if not all([matplotlib.colors.is_color_like(s)
                            for s in plotcolors]):
                    plotcolors = None
                    logger.warning("writeRfactorPdf: Specified colors not "
                                   "recognized, reverting to default colors")
            if plotcolors is None:
                axs[0].plot(theo[:, 0], theo[:, 1], label='Theoretical')
                axs[0].plot(exp[:, 0], exp[:, 1], label='Experimental')
                axs[1].plot(ytheo[:, 0], ytheo[:, 1], label='Theoretical')
                axs[1].plot(yexp[:, 0], yexp[:, 1], label="Experimental")
            else:
                axs[0].plot(theo[:, 0], theo[:, 1], label='Theoretical',
                            color=plotcolors[0])
                axs[0].plot(exp[:, 0], exp[:, 1], label='Experimental',
                            color=plotcolors[1])
                axs[1].plot(ytheo[:, 0], ytheo[:, 1], label='Theoretical',
                            color=plotcolors[0], linewidth=0.75)
                axs[1].plot(yexp[:, 0], yexp[:, 1], label="Experimental",
                            color=plotcolors[0], linewidth=0.75)
            axs[1].plot(dy[:, 0], dy[:, 1], label="\u0394Y", color="black",
                        linewidth=0.5)
            axs[1].fill_between(dy[:, 0], dy[:, 1], 0., facecolor='grey',
                                alpha=0.5)
            axs[2].plot(idysq[:, 0], idysq[:, 1], color="black",
                        drawstyle="steps-mid")

            axs[0].annotate(name, namePos, fontsize=10)
            axs[0].annotate("R = {:.4f}".format(rfact), rPos, fontsize=10)
            axs[0].legend()
            axs[1].legend()

        for fig in figs:
            pdf.savefig(fig)
            plt.close(fig)
    except Exception:
        logger.error("writeRfactorPdf: Error while writing analysis pdf: ",
                     exc_info=True)
    finally:
        pdf.close()
        logger.setLevel(loglevel)
    return
