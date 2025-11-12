"""Module iorfactor of viperleed.calc.files.

Functions for reading and writing files relevant to the
R-factor calculation.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-08-19'
__license__ = 'GPLv3+'

import logging
import os
import re

import fortranformat as ff
import numpy as np

from viperleed.calc.classes.rparams.special.energy_range import EnergyRange
from viperleed.calc.files.beams import writeAUXEXPBEAMS
from viperleed.calc.files.ivplot import plot_iv
from viperleed.calc.lib import leedbase
from viperleed.calc.lib.log_utils import at_level
from viperleed.calc.lib.matplotlib_utils import CAN_PLOT
from viperleed.calc.lib.matplotlib_utils import log_without_matplotlib
from viperleed.calc.lib.matplotlib_utils import prepare_matplotlib_for_calc
from viperleed.calc.lib.matplotlib_utils import skip_without_matplotlib
from viperleed.calc.lib.version import Version

if CAN_PLOT:
    prepare_matplotlib_for_calc()
    from matplotlib import pyplot as plt
    from matplotlib import ticker as plticker
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import is_color_like

logger = logging.getLogger(__name__)


# How many extra points to take at the boundaries to prevent wasting
# data for interpolation of theoretical beams? The FORTRAN code likes
# to have 2 intervals left and right. Notice that this value is
# specific of the FORTRAN code, which uses 3rd order interpolation.
# It is made into a const so we can change it in only one spot should
# I (@michele-riva) have gotten it wrong. It SHOULD NOT BE USED for
# other purposes.
_N_EXPAND_THEO = 2


class RfactorError(Exception):
    """Base exception for R-factor calculations."""


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
        return (0, 0, 0), 0, []
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


def check_theobeams_energies(rpars, theobeams):
    """Complain if the energies in theobeams are inconsistent with rpars."""
    theo_grid = sorted_energies_from_beams(theobeams)
    theo_energies = EnergyRange.from_sorted_grid(theo_grid)
    if not theo_energies.contains(rpars.THEO_ENERGIES, ignore_step=True):
        raise ValueError(
            f'theobeams has theo_energies={theo_energies}, which does not '
            'contain all the energies in the THEO_ENERGIES parameter '
            f'({rpars.THEO_ENERGIES}). Did you load the wrong calculation?'
            )


def prepare_rfactor_energy_ranges(rpars, theobeams=None, for_error=False,
                                  n_expand=_N_EXPAND_THEO):
    f"""Return EnergyRange objects for experiment, theory and iv_shifts.

    Parameters
    ----------
    rpars : Rparams
        The object holding information about the current PARAMETERS.
        Attributes used: THEO_ENERGIES, expbeams, IV_SHIFT_RANGE
    theobeams : Sequence or None, optional
        Elements are Beam objects. They are the beams loaded from
        the calculation against which experimental data should be
        compared. If given, it is checked for consistency with the
        rpars.THEO_ENERGIES attribute. theobeams are considered to
        be consistent if rpars.THEO_ENERGIES is a subset of the
        energies in theobeams with the same step size. No check is
        performed if not given or None. Default is None.
    for_error : bool, optional
        Whether the R-factor calculation is meant to be used for
        error curves. Default is False.
    n_expand : int, optional
        How many steps should the theory range be expanded at both
        ends not to disregard available values that may be useful
        for interpolation. This is on top of the expansion that is
        considered due to shifting. Default is {_N_EXPAND_THEO}.

    Returns
    -------
    exp_range : EnergyRange
        The range of energies of the experimental beams.
    theo_range : TheoEnergies
        The relevant range of energies. This is derived from the
        rpars.THEO_ENERGIES value. It is appropriately expanded to
        account for relative shifting of the curves. Additionally,
        it is also expanded by n_expand steps to the left and right.
    iv_shift : IVShiftRange
        The range of shifts to be applied. If for_error, this is
        fixed to rpars.best_v0r value. Otherwise, it is the rpars
        attribute IV_SHIFT_RANGE. Notice that, in the latter case,
        the attribute is modified in place if it does not yet have
        a step defined. Its step is set to interp_step. Notice also
        that, as a consequence of this change, rpars.IV_SHIFT_RANGE
        may be expanded so that its bounds are integer multiples of
        the step.
    interp_step : float
        The energy step to be used for interpolation. Notice that
        this is ALWAYS taken from rpars.IV_SHIFT_RANGE.

    Raises
    ------
    ValueError
        If theobeams was given, but it is inconsistent with rpars.
    RfactorError
        If rpars.THEO_ENERGIES has any undefined attribute, if there
        are no beams in rpars.expbeams, or if rpars.IV_SHIFT_RANGE
        has any of its bounds undefined.
    RfactorError
        When for_error and either (i) no rpars.best_v0r is present,
        or (ii) the rpars.best_v0r value is inconsistent with the
        interpolation step, or (iii) rpars.IV_SHIFT_RANGE has no
        step.
    """
    if theobeams is not None:
        check_theobeams_energies(rpars, theobeams)  # May ValueError
    if not rpars.THEO_ENERGIES.defined:
        raise RfactorError('Cannot run an R-factor calculation before '
                           f'rpars.THEO_ENERGIES={rpars.THEO_ENERGIES} '
                           'has stop, start, and step defined')
    if not rpars.expbeams:
        raise RfactorError('No experimental beams loaded')

    # Prepare experimental range
    exp_grid = sorted_energies_from_beams(rpars.expbeams)
    exp_range = EnergyRange.from_sorted_grid(exp_grid)

    # Now handle iv_shift. Notice that _prepare_iv_shift_range always
    # sets the step to the one of IV_SHIFT_RANGE to ensure consistency
    iv_shift = _prepare_iv_shift_range(rpars, exp_range, for_error)

    small_step = min(exp_range.step, rpars.THEO_ENERGIES.step)
    if not iv_shift.is_fixed and iv_shift.step > small_step:
        # Warn only if we're not fixed. If we are, it may just be
        # that the user gave us this specific step value because
        # it is exactly an integer multiple of the fixed value
        # they desired
        logger.warning(f'Using interpolation step {iv_shift.step} eV, which '
                       'is coarser than the smallest step between experiment '
                       f'and theory ({small_step} eV). Consider using '
                       f'{small_step} or not defining a step for the '
                       'IV_SHIFT_RANGE parameter')

    # Finally, the relevant range of theory energies
    theo_range = _find_common_theory_experiment_range(rpars.THEO_ENERGIES,
                                                      exp_range, iv_shift,
                                                      n_expand)
    return exp_range, theo_range, iv_shift, iv_shift.step


def _prepare_iv_shift_range(rpars, experiment, for_error):
    """Return an iv_shift appropriate for an R-factor calculation."""
    if for_error and rpars.best_v0r is None:
        raise RfactorError(
            'Cannot run an error calculation if rpars.best_v0r is '
            'undefined. Did you forget to run a normal R-factor first?'
            )
    if not rpars.IV_SHIFT_RANGE.has_bounds:
        raise RfactorError('Cannot run an R-factor calculation if '
                           f'rpars.IV_SHIFT_RANGE={rpars.IV_SHIFT_RANGE} '
                           'misses one of its bounds')

    # Set up the step. This also ensures that the bounds are consistent
    if not for_error:
        # Notice that here we purposely do not take as step
        # min(step, IV_SHIFT_RANGE.step). This would otherwise
        # require us to potentially force a modification of the
        # step of IV_SHIFT_RANGE and, likely, its bounds (because
        # they are integer multiples). This seems like too much
        # fiddling with the user input. See also discussion at
        # viperleed/commit/d4626116f11fb0bf9bef6c228413047a7207d441
        step = min(rpars.THEO_ENERGIES.step, experiment.step)
        try:
            rpars.IV_SHIFT_RANGE.set_undefined_step(step)
        except RuntimeError:  # Was fixed and we're trying to unfix it
            start = rpars.IV_SHIFT_RANGE.min
            raise RfactorError(
                f'Cannot fix IV_SHIFT_RANGE to {start}. The automatic step '
                f'from experiment and theory ({step} eV) is inappropriate: '
                f'{start} must be an integer  multiple of step={step}. '
                'Provide an explicit step by setting IV_SHIFT_RANGE'
                ) from None
        return rpars.IV_SHIFT_RANGE

    # In the for_error case, always make sure to
    # use the step stored in rpars.IV_SHIFT_RANGE
    if not rpars.IV_SHIFT_RANGE.has_step:
        raise RfactorError(
            'Cannot run an error calculation if rpars.IV_SHIFT_RANGE has '
            'no step defined. Did you forget to run a normal R-factor first?'
            )
    iv_shift = rpars.IV_SHIFT_RANGE.fixed(rpars.best_v0r)
    try:
        iv_shift.set_undefined_step(rpars.IV_SHIFT_RANGE.step)
    except RuntimeError as exc:
        # We have tried to unfix a fixed range: best_v0r
        # is for the wrong interpolation step
        raise RfactorError(f'rpars.best_v0r={rpars.best_v0r} is inconsistent '
                           f'with the interpolation step ({iv_shift.step} eV).'
                           ' Did you change the step of rpars.IV_SHIFT_RANGE '
                           'since the last R-factor calculation?') from exc
    return iv_shift


def _find_common_theory_experiment_range(theory, experiment,
                                         iv_shift, n_expand):
    """Return the range of theory energies in common with experiment."""
    # This is a little bit tricky. Here is what we want:
    #             |-----------------------|    theory
    #                     |-------------|      experiment
    #                     .             .
    #  |------------------vvvvv|        .      theory, shift max left
    #                |----vvvvvvvvvvvvvvv----| theory, shift max right
    # The 'v's mark regions of the shifted theory in common
    # with experiments, on the 'shifted' scales. (These are
    # the extremes. We want all the continuous shifts in
    # between too.) However, we want the union of all the
    # 'v' parts ON THE ORIGINAL ENERGY SCALE of the theory,
    # i.e.,
    #             |-----------------------|    theory
    #             |------------------vvvvv|    left, on same scale
    #             |----vvvvvvvvvvvvvvv----|    right, on same scale
    #             |----vvvvvvvvvvvvvvvvvvv|    Union of all shifts
    # Rather than doing this mess, we can equivalently EXPAND the
    # experimental range: on the left by the max right shift, on
    # the right by the max left shift, like so:
    #                     |-------------|            experiment
    #                  |RR|-------------|LLLLLLLLLL| expanded exp.
    #             |-----------------------|          theory
    #             |----vvvvvvvvvvvvvvvvvvv|          common range
    expanded_exp = experiment.copy()
    expanded_exp.start -= iv_shift.stop  # Expand when stop > 0
    expanded_exp.stop -= iv_shift.start  # Expand when start < 0

    # Expand by n_expand left & right, finally intersect
    theory_larger = theory.expanded_by(n_expand)
    try:
        return theory_larger.intersected(expanded_exp)
    except ValueError:  # No intersection
        raise RfactorError('Cannot run R-factor calculation: No common '
                           'energies between experiment and theory') from None


def sorted_energies_from_beams(beams):
    """Return a list of sorted energies from a list of Beam objects."""
    return sorted({e for beam in beams for e in beam.intens})


def writeWEXPEL(sl, rp, theobeams, filename="WEXPEL", for_error=False):
    """Write input file WEXPEL for R-factor calculation.

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
    for_error : bool, optional
        Whether the R-factor calculation is used to determine error
        curves. The default is False.

    Returns
    -------
    None.

    Raises
    ------
    RfactorError
        If selected R-factor type is not supported by TensErLEED
    """
    (_, theo_range,
     iv_shift, vincr) = prepare_rfactor_energy_ranges(rp, theobeams,
                                                      for_error,
                                                      n_expand=_N_EXPAND_THEO)

    # find correspondence experimental to theoretical beams:
    beamcorr = leedbase.getBeamCorrespondence(sl, rp)
    # integer & fractional beams
    iorf = []
    for (i, beam) in enumerate(rp.expbeams):
        if beam.hk[0] % 1.0 != 0.0 or beam.hk[1] % 1.0 != 0.0:
            iorf.append(2)
        else:
            iorf.append(1)
    iorf.extend([0]*(len(rp.ivbeams)-len(rp.expbeams)))

    f72 = ff.FortranRecordWriter('F7.2')
    if rp.TL_VERSION < Version('1.7.0'):
        beam_formatter = ff.FortranRecordWriter('25I3')
    else:
        beam_formatter = ff.FortranRecordWriter('25I4')
    i3 = ff.FortranRecordWriter('I3')

    # EMIN/EMAX are just two bounds to select which theory beams are
    #     read in. All energies outside the closed [EMIN, EMAX] range
    #     are not read in (and do not fill up any array). Notice that
    #     we take the maximum a bit larger, just to make sure that
    #     integer divisions in the FORTRAN code don't discard the
    #     upper limit.
    # EINCR is the interpolation step.
    # LIMFIL is number of consecutive input files. Hardcode to 1.
    # IPR is output formatting. Hardcode to 0.
    # **IMPORTANT**: VINCR should be an integer multiple of EINCR.
    # If the two are not multiples, the V0R shifts are wrong. Here,
    # for simplicity, we take them to be the exact same.
    output = f'''\
 &NL1
 EMIN={f72.write([theo_range.min]):>9},
 EMAX={f72.write([theo_range.max + 0.1 * vincr]):>9},
 EINCR={f72.write([vincr]):>8},
 LIMFIL=      1,
 IPR=         0,
 VI={f72.write([rp.V0_IMAG]):>11},
 V0RR=      0.0,
 V01={f72.write([iv_shift.start]):>10},
 V02={f72.write([iv_shift.stop]):>10},
 VINCR={f72.write([vincr]):>8},
 ISMOTH={i3.write([rp.R_FACTOR_SMOOTH]):>7},
 EOT=         0,
 PLOT=        1,
 GAP=         0,
 &END
'''
    output += beam_formatter.write([n+1 for n in beamcorr]) + '\n'
    if len(beamcorr) % 25 == 0:
        output += "\n"
    for i in range(0, 2):  # redundant since indices are already taken care of
        output += beam_formatter.write(
            [n+1 for n in range(len(rp.expbeams))]) + '\n'
        if len(rp.expbeams) % 25 == 0:
            output += '\n'
    output += beam_formatter.write(iorf) + '\n'
    if len(iorf) % 25 == 0:
        output += '\n'
    output += '&NL2\n'
    output += ' NSSK=    0,\n'
    if str(rp.R_FACTOR_TYPE) == 'pendry':
        output += ' WR=      0.,0.,1.,\n'
    elif str(rp.R_FACTOR_TYPE) == 'r2':
        output += ' WR=      1.,0.,0.,\n'
    elif str(rp.R_FACTOR_TYPE) == 'zj':
        output += ' WR=      0.,1.,0.,\n'
    else:
        msg = (
            f'R factor type {rp.R_FACTOR_TYPE} not supported by '
            'TensErLEED backend.'
        )
        logger.error(msg)
        raise RfactorError(msg)
    output += '''\
 &END
 &NL3
 NORM=           1,
 INTMAX=    999.99,
 PLSIZE=   1.0,1.0,
 XTICS=         50,
 &END
 '''
    auxexpbeams = writeAUXEXPBEAMS(rp.expbeams, header=rp.systemName,
                                   write=True, numbers=False)
    output += auxexpbeams + '\n'
    # information about gaps in the experimental spectra would go here
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except Exception:
        logger.error(f'Failed to write {filename}')
        raise
    logger.debug(f'Wrote to R-factor input file {filename} successfully')


def largest_nr_grid_points(rpars, theobeams, for_error,
                           n_expand=_N_EXPAND_THEO):
    """Return the largest possible number of grid points."""
    # The number of grid point is the quantity used to set up the
    # dimensions of the arrays that will contain the IV curves.
    # It should be at least as large as:
    # - the number of experiment energies, on the original grid
    # - the number of theory energies, on the original grid
    # - the number of experiment energies, on the interpolated grid
    # - the number of theory energies, on the interpolated grid
    # - in principle, there is also some relation to the number
    #   of energies in GAPS, if there's gaps in the data. Since
    #   we don't allow this, we do not care. Should we ever change
    #   our mind, we'd have to look more closely at the code.

    # Notice that, in principle, we would not need to keep extra room
    # for shifting the theory relative to experiments. The FORTRAN code
    # translates the 'shifted' indices to indices in the interpolated
    # version of the theory beams. However, it seems cleaner to use the
    # very same logic as in WEXPEL. This way we only need to maintain
    # one piece of code.
    (experiment, theory,
     _, interp_step) = prepare_rfactor_energy_ranges(rpars, theobeams,
                                                     for_error, n_expand)
    interp_exp = EnergyRange(experiment.start, experiment.stop, interp_step)
    section = 'refcalc' if not for_error else 'superpos'

    # Make sure to use the whole stored range of energies if available
    # (which corresponds also to what is in fd.out) and not just the
    # overlapping one!
    if rpars.theobeams[section]:
        interp_theo = EnergyRange(
            min([min(b.energies) for b in rpars.theobeams[section]]),
            max([max(b.energies) for b in rpars.theobeams[section]]),
            interp_step)
    else:
        interp_theo = EnergyRange(theory.start, theory.stop, interp_step)

    n_max = max((experiment.n_energies,
                 theory.n_energies,
                 interp_exp.n_energies,
                 interp_theo.n_energies))
    return round(np.ceil(n_max * 1.1))  # 10% headroom, just in case


def writeRfactPARAM(rp, theobeams, for_error=False, only_vary=None):
    """Generate the PARAM file for the R-factor calculation.

    Parameters
    ----------
    rp : Rparams
        The run parameters.
    theobeams : list of Beam
        The theoretical beams, containing I(V) data.
    for_error : bool, optional
        Whether this R-factor calculation is used to produce error
        curves. Default is False.
    only_vary : list or None, optional
        Which parameter should be varied. Items are SearchPar objects.
        This argument is used only if for_error. In that case, it is
        used solely to compute how many distinct parameters values
        are expected. If not given or None, all the parameters are
        used. Default is False.

    Raises
    ------
    RfactorError
        If this function is called before writeWEXPEL.
    """
    if not rp.IV_SHIFT_RANGE.has_step:
        raise RfactorError('Cannot writeRfactPARAM without interpolation '
                           'step. Did you forget to call writeWEXPEL first?')
    ngrid = largest_nr_grid_points(rp, theobeams, for_error)

    n_var = 1
    if for_error:
        if not only_vary:
            logger.warning('Rfactor PARAM for error: Parameters under '
                           'variation not passed.')
            only_vary = [sp for sp in rp.searchpars
                         if sp.atom in rp.search_atlist]
        n_var = max([sp.steps for sp in only_vary])
    output = f'''
C  MNBED  : number of beams in experimental spectra before averaging
C  MNBTD  : number of beams in theoretical spectra before averaging

      PARAMETER (MNBED = {len(rp.expbeams)}, MNBTD = {len(theobeams)})

C  MNET   : number of energies in theoretical beam at time of reading in
C  MNGP  : greater equal number of grid points in energy working grid (ie after
C           interpolation)
C  MNS    : number of geometries including those, that are skipped

      PARAMETER (MNET = {rp.THEO_ENERGIES.n_energies}, MNGP = {ngrid})
      PARAMETER (MNS = {n_var})

C  MNGAP  : number of gaps in the experimental spectra (NOTE: if there are no
C           gaps in the spectra set MNGAP to 1 to avoid zero-sized arrays)

      PARAMETER (MNGAP = 1)
'''
    # write PARAM
    try:
        with open('PARAM', 'w') as wf:
            wf.write(output)
    except Exception:
        logger.error('Failed at writing PARAM file for R-factor calculation.')
        raise


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
            return [], []
        except PermissionError:
            logger.error("read_rfactor_columns: Cannot open file {}. Aborting."
                         .format(fname))
            return [], []

        cols = np.array([[float(col) for col in line.split()] for line in f])
        f.close()
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


@log_without_matplotlib(logger, msg='Skipping R-factor plotting.')
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
    xyTheo, xyExp = read_rfactor_columns(cols_dir=colsDir)
    labels, rfacs = zip(*beams)
    rfac_str = ["R = {:.4f}".format(r) for r in rfacs]
    plot_iv([xyTheo, xyExp], outName, legends=['Theoretical', 'Experimental'],
            labels=labels, annotations=rfac_str, formatting=formatting)

    if not analysisFile:
        return

    figs, figsize, namePos, oritick, plotcolors, rPos, xlims, ylims = prepare_analysis_plot(formatting, xyExp, xyTheo)

    try:
        # Pylint can't tell that we will not execute this,
        # as per decorator, if we fail to import matplotlib
        # pylint: disable-next=possibly-used-before-assignment
        pdf = PdfPages(analysisFile)
    except PermissionError:
        logger.error("writeRfactorPdf: Cannot open file {}. Aborting."
                     .format(analysisFile))
        return

    # The following will spam the logger with debug messages; disable.
    with at_level(logger, logging.INFO):
        try:
            for i, (name, rfact, theo, exp) in enumerate(zip(*zip(*beams),
                                                             xyTheo, xyExp)):
                if len(exp) == 0:
                    continue
                ytheo = leedbase.getYfunc(theo, v0i)
                yexp = leedbase.getYfunc(exp, v0i)
                plot_analysis(exp, figs, figsize, name, namePos, oritick,
                              plotcolors, rPos, rfact, theo, xlims, yexp,
                              ylims, ytheo, v0i)
            for fig in figs:
                pdf.savefig(fig)
                # Pylint can't tell that we will not execute this,
                # as per decorator, if we fail to import matplotlib
                # pylint: disable-next=possibly-used-before-assignment
                plt.close(fig)
        except Exception:
            logger.error("writeRfactorPdf: Error while writing analysis pdf: ",
                         exc_info=True)
        finally:
            pdf.close()


@skip_without_matplotlib
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
    # Pylint can't tell that we will not execute this,
    # as per decorator, if we fail to import matplotlib
    # pylint: disable-next=possibly-used-before-assignment
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


def plot_analysis(exp, figs, figsize, name, namePos, oritick, plotcolors, rPos, rfact, theo, xlims, yexp, ylims, ytheo, v0i):
    fig, axs = plt.subplots(3, figsize=figsize,
                            squeeze=True)
    fig.subplots_adjust(left=0.06, right=0.94,
                        bottom=0.07, top=0.98,
                        wspace=0, hspace=0.08)
    figs.append(fig)
    [ax.set_xlim(*xlims) for ax in axs]
    axs[0].set_ylim(*ylims)
    if v0i is not None: # if is None, don't scale
        # Y function plot limits:
        # Y osciallates between +/- 1/(2*abs(v0i))
        y_range_uni_direc = 1/(2*abs(v0i)) * 1.1 # 10% spacing on either side
        axs[1].set_ylim(-1*y_range_uni_direc, y_range_uni_direc)
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

    eps = 1e-10 # needed to avoid division by 0. Only important for plotting.
    norm_y_squares = np.array(
        [[dysq[j, 0], (dysq[j, 1]
                       / (y_theo_sq[j, 1] + y_exp_sq[j, 1] + eps))]
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
        if not all(is_color_like(s) for s in plotcolors):
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
                    color=plotcolors[1], linewidth=0.75)
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


@log_without_matplotlib(logger, msg='Skipping R-factor plotting.')
def writeRfactorPdf_new(n_beams, labels, rfactor_beams,
                        energies, id_start,
                        n_E_beams,
                        int_1, int_2, y_1, y_2 ,
                        outName='Rfactor_plots.pdf',
                        analysisFile='', v0i = 0., formatting=None):
    # after applying the V0r shift outside, the id_start and n_E_beams should be same for experiment and theory
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
        # Pylint can't tell that we will not execute this,
        # as per decorator, if we fail to import matplotlib
        # pylint: disable-next=possibly-used-before-assignment
        pdf = PdfPages(analysisFile)
    except PermissionError:
        logger.error("writeRfactorPdf: Cannot open file {}. Aborting."
                     .format(analysisFile))
        return
    # The following will spam the logger with debug messages; disable.
    with at_level(logger, logging.INFO):
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
                              theo, xlims, y_exp, ylims, y_theo, v0i)
            for fig in figs:
                pdf.savefig(fig)
                # Pylint can't tell that we will not execute this,
                # as per decorator, if we fail to import matplotlib
                # pylint: disable-next=possibly-used-before-assignment
                plt.close(fig)
        except Exception:
            logger.error("writeRfactorPdf: Error while writing analysis pdf: ",
                         exc_info=True)
        finally:
             pdf.close()


def beamlist_to_array(beams):
    # turn list of Beam objects into an array of intensities

    n_beams = len(beams)
    energies = sorted_energies_from_beams(beams)
    in_grid = np.array(energies)
    n_E = in_grid.shape[0]

    # fill with NaNs as default value
    beam_arr = np.full([n_E, n_beams], fill_value=np.nan)

    id_start = np.zeros(n_beams, dtype=np.int32)
    n_E_beams = np.zeros(n_beams, dtype=np.int32)

    for i, b in enumerate(beams):
        # write beams into colums of beam_arr
        id_start[i] = np.where(np.isclose(energies, min(b.intens.keys())))[0][0]
        n_E_beams[i] = len(b.intens)
        beam_arr[id_start[i]: id_start[i] + n_E_beams[i], i] = list(b.intens.values())


    return in_grid, id_start, n_E_beams, beam_arr

