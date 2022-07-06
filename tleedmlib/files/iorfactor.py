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
from viperleed.tleedmlib.files.ivplot import plot_iv

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
    # finds line at end of ROUT file that states best v0r and corresponding R-factor
    while "AVERAGE R-FACTOR =" not in line and i+1 < len(lines):
        i += 1
        line = lines[-i]
    if line == "":
        return 0, 0, []
    rfac = 0
    rfac_int = -1
    rfac_frac = -1
    v0rshift = 0
    # writes best V0r to v0rshift and best R-factor to rfac
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
    # now read the R-factors per beam only at best V0r
    rfaclist = []
    for line in [li for li in lines if len(li) >= 70]:
        line = line.strip()
        if line.endswith("<---"):
            line = line[:-4]

        # rfactor.f reserves 5A4 i.e. 20 characters for the beam name
        values = line[19:].split()
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
    step = min(0.5, expEnergies[1]-expEnergies[0],
               theoEnergies[1]-theoEnergies[0])
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
    if rp.TL_VERSION < 1.7:
        beam_formatter = ff.FortranRecordWriter("25I3")
    else:
        beam_formatter = ff.FortranRecordWriter("25I4")
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
    output += (beam_formatter.write([n+1 for n in beamcorr]) + "\n")
    if len(beamcorr) % 25 == 0:
        output += "\n"
    for i in range(0, 2):  # redundant since indices are already taken care of
        output += beam_formatter.write(
            [n+1 for n in range(0, len(rp.expbeams))]) + "\n"
        if len(rp.expbeams) % 25 == 0:
            output += "\n"
    output += beam_formatter.write(iorf) + "\n"
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


def read_rfactor_columns(cols_dir=''):
    """
    Reads data from the theo.column and exp.column files in a given directory.

    Parameters
    ----------
    cols_dir : str, optional
        The directory to read from. The default is ''.

    Returns
    -------
    list [theo_beams, exp_beams]
        Both the theoretical and the experimental beams are formatted as lists
        of 2D numpy arrays, each of which has the form [[en1, intens1], ...]
    """

    fnames = ['theo.column', 'exp.column']
    xxyy = []
    for fname in fnames:
        try:
            f = open(os.path.join(cols_dir, fname), 'r')
        except FileNotFoundError:
            logger.error("read_rfactor_columns: File {} not found. Aborting."
                         .format(fname))
            return []
        except PermissionError:
            logger.error("read_rfactor_columns: Cannot open file {}. Aborting."
                         .format(fname))
            return []

        cols = np.array([[float(col) for col in line.split()] for line in f])
        xy = np.split(cols, np.shape(cols)[1]/2, axis=1)
        # xy is now a list of 2D arrays.
        # Each array has the form [[en1, intens1], ...]
        #
        # for each beam, get rid of the points that have (en, intens) = (0, 0)
        # so that they don't screw up the plots later
        xy = [coords[~np.all(coords < 1e-3, axis=1)] for coords in xy]
        if xy:
            xxyy.append(xy)
        else:
            logger.warning("File " + fname + " contains no usable data.")
    # xxyy now contains first the theoretical, then the experimental beams
    return xxyy


def writeRfactorPdf(beams, colsDir='', outName='Rfactor_plots.pdf',
                    analysisFile='', v0i=0., formatting=None):
    '''
    Creates a single PDF file containing the plots of R-factors, using plot_iv.
    If analysisFile is defined, a second 'analysis' PDF will be generated.

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

    xyTheo, xyExp = read_rfactor_columns(cols_dir=colsDir)
    labels, rfacs = zip(*beams)
    rfac_str = ["R = {:.4f}".format(r) for r in rfacs]
    plot_iv([xyTheo, xyExp], outName, legends=['Theoretical', 'Experimental'],
            labels=labels, annotations=rfac_str, formatting=formatting)

    if not analysisFile:
        return

    figs, figsize, namePos, oritick, plotcolors, rPos, xlims, ylims = prepare_analysis_plot(formatting, xyExp, xyTheo)

    try:
        pdf = PdfPages(analysisFile)
    except PermissionError:
        logger.error("writeRfactorPdf: Cannot open file {}. Aborting."
                     .format(analysisFile))
        return
    # the following will spam the logger with debug messages; disable.
    loglevel = logger.level
    logger.setLevel(logging.INFO)
    try:
        for i, (name, rfact, theo, exp) in enumerate(zip(*zip(*beams),
                                                         xyTheo, xyExp)):

            if len(exp) == 0:
                continue
            ytheo = getYfunc(theo, v0i)
            yexp = getYfunc(exp, v0i)

            plot_analysis(exp, figs, figsize, name, namePos, oritick, plotcolors, rPos, rfact, theo, xlims, yexp, ylims,
                          ytheo)

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


def prepare_analysis_plot(formatting, xyExp, xyTheo):
    # write R-factor analysis
    # find min and max values of x and y for plotting all curves
    # on the same horizontal scale and leaving a little y space for the legend
    xmin = min(min(xy[:, 0]) for xy in [*xyTheo, *xyExp] if len(xy) != 0)
    xmax = max(max(xy[:, 0]) for xy in [*xyTheo, *xyExp] if len(xy) != 0)
    ymin = min(min(xy[:, 1]) for xy in [*xyTheo, *xyExp] if len(xy) != 0)
    ymax = max(max(xy[:, 1]) for xy in [*xyTheo, *xyExp] if len(xy) != 0)
    dy = ymax - ymin
    # set ticks spacing to 50 eV and round the x limits to a multiple of it
    tick = 50
    oritick = plticker.MultipleLocator(base=tick)
    xlims = (np.floor(xmin / tick) * tick,
             np.ceil(xmax / tick) * tick)
    dx = xlims[1] - xlims[0]
    # set formatting parameters
    plotcolors = []
    if formatting is not None:
        if 'colors' in formatting:
            plotcolors = formatting['colors']
    figsize = (5.8, 8.3)
    figs = []
    ylims = (ymin - 0.02 * dy, ymax + 0.22 * dy)
    namePos = (xlims[0] + 0.45 * dx, ylims[1] - 0.1 * dy)
    rPos = (namePos[0], namePos[1] - 0.085 * dy)
    return figs, figsize, namePos, oritick, plotcolors, rPos, xlims, ylims


def plot_analysis(exp, figs, figsize, name, namePos, oritick, plotcolors, rPos, rfact, theo, xlims, yexp, ylims, ytheo):
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
    axs[2].set_ylabel(r'$\sum\Delta Y^2 / \left(Y_1^2 + Y_2^2\right)$')
    axs[2].set_xlabel("Energy (eV)")
    y_theo_sq = np.copy(ytheo)
    y_theo_sq[:, 1] = y_theo_sq[:, 1] ** 2
    y_exp_sq = np.copy(yexp)
    y_exp_sq[:, 1] = y_exp_sq[:, 1] ** 2
    dy = np.array([(ytheo[j, 0], yexp[j, 1] - ytheo[j, 1])
                   for j in range(0, min(len(ytheo), len(yexp)))])
    dysq = np.copy(dy)
    dysq[:, 1] = dysq[:, 1] ** 2
    norm_y_squares = np.array(
        [[dysq[j, 0], (dysq[j, 1]
                       / (y_theo_sq[j, 1] + y_exp_sq[j, 1]))]
         for j in range(len(dysq))])
    sum_norm_y_squares = np.array([norm_y_squares[0]])
    for j in range(1, len(norm_y_squares)):
        sum_norm_y_squares = (
            np.append(sum_norm_y_squares,
                      [[norm_y_squares[j, 0],
                        (sum_norm_y_squares[j - 1, 1]
                         + norm_y_squares[j, 1])]],
                      axis=0))
    axs[1].plot(xlims, [0., 0.], color='grey', alpha=0.2)
    if not plotcolors:
        if not all([matplotlib.colors.is_color_like(s)
                    for s in plotcolors]):
            plotcolors = []
            logger.warning("writeRfactorPdf: Specified colors not "
                           "recognized, reverting to default colors")
    if not plotcolors:
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
    axs[2].plot(sum_norm_y_squares[:, 0], sum_norm_y_squares[:, 1],
                color="black", drawstyle="steps-mid")
    axs[0].annotate(name, namePos, fontsize=10)
    axs[0].annotate("R = {:.4f}".format(rfact), rPos, fontsize=10)
    axs[0].legend()
    axs[1].legend()


def writeRfactorPdf_new(n_beams, labels, rfactor_beams,
                        energies, id_start,
                        n_E_beams,
                        int_1, int_2, y_1, y_2 ,
                        outName='Rfactor_plots.pdf',
                        analysisFile='', v0i = 0., formatting=None):

    # after applying the V0r shift outside, the id_start and n_E_beams should be same for experiment and theory
    global plotting
    if not plotting:
        logger.debug("Necessary modules for plotting not found. Skipping "
                     "R-factor plotting.")
        return


    # get data
    exp_xy = []
    theo_xy = []
    for i in range(n_beams):
        xy = np.empty([n_E_beams[i], 2])
        xy[:, 0] = energies[id_start[i] -1: id_start[i] + n_E_beams[i] -1]
        xy[:, 1] = int_1[id_start[i] -1: id_start[i] + n_E_beams[i] -1, i]
        # normalize to max of beam:
        xy[:, 1] /= np.nanmax(xy[:, 1])
        exp_xy.append(xy)


        xy = np.empty([n_E_beams[i], 2]) # want this at same range as exp only!
        xy[:, 0] = energies[id_start[i] -1: id_start[i] + n_E_beams[i] -1]
        xy[:, 1] = int_2[id_start[i] -1: id_start[i] + n_E_beams[i] -1, i]
        # normalize to max of beam:
        xy[:, 1] /= np.nanmax(xy[:, 1])
        theo_xy.append(xy)

    data = [theo_xy, exp_xy]

    rfac_str = ["R = {:.4f}".format(r) for r in rfactor_beams]
    plot_iv(data, outName, legends=['Theoretical', 'Experimental'],
            labels=labels, annotations=rfac_str, formatting=formatting)

    if not analysisFile:
        return

    figs, figsize, namePos, oritick, plotcolors, rPos, xlims, ylims = \
        prepare_analysis_plot(formatting, exp_xy, theo_xy)


    try:
        pdf = PdfPages(analysisFile)
    except PermissionError:
        logger.error("writeRfactorPdf: Cannot open file {}. Aborting."
                     .format(analysisFile))
        return
    # the following will spam the logger with debug messages; disable.
    loglevel = logger.level
    logger.setLevel(logging.INFO)

    # proper minus character
    labels = [label.replace("-", "−") for label in labels]

    try:
        for i in range(n_beams):
            exp = exp_xy[i]
            theo = theo_xy[i]
            beam_energies = energies[id_start[i] -1: id_start[i] + n_E_beams[i] -1]
            y_exp = np.empty([n_E_beams[i], 2])
            y_theo = np.empty([n_E_beams[i], 2])
            y_exp[:, 0] = beam_energies
            y_exp[:, 1] = y_1[id_start[i] -1: id_start[i] + n_E_beams[i] -1, i]
            y_theo[:, 0] = beam_energies
            y_theo[:, 1] = y_2[id_start[i] -1: id_start[i] + n_E_beams[i] -1, i]

            plot_analysis(exp, figs, figsize, labels[i], namePos, oritick, plotcolors, rPos, rfactor_beams[i],
                          theo, xlims, y_exp, ylims, y_theo)

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


def beamlist_to_array(beams):
    # turn list of Beam objects into an array of intensities

    n_beams = len(beams)
    energies = sorted({e for b in beams for e in b.intens})
    in_grid = np.array(energies)
    n_E = in_grid.shape[0]
    
    # fill with NaNs as default value
    beam_arr = np.full([n_E, n_beams], fill_value=np.NaN)

    id_start = np.int32(np.zeros([n_beams]))
    n_E_beams = np.int32(np.zeros([n_beams]))

    for i, b in enumerate(beams):
        # write beams into colums of beam_arr
        id_start[i] = np.where(np.isclose(energies, min(b.intens.keys())))[0][0]
        n_E_beams[i] = len(b.intens)
        beam_arr[id_start[i]: id_start[i] + n_E_beams[i], i] = list(b.intens.values())


    return in_grid, id_start, n_E_beams, beam_arr

