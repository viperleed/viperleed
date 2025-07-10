"""Module iodeltas of viperleed.calc.files.

Defines functions for reading and writing files relevant to the
delta-amplitudes calculation.
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
from pathlib import Path
import shutil

import fortranformat as ff
import numpy as np

from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.files import beams
from viperleed.calc.lib.version import Version

logger = logging.getLogger(__name__)


def collect_static_input_files(slab, rpars):
    """Collect the contents of input files common to all delta calculations.

    Parameters
    ----------
    slab : Slab
        The slab for which delta calculations are to be performed.
    rpars : Rparams
        The parameters corresponding to `slab`.

    Returns
    -------
    delta_basic : str
        The contents of the input file in common to all deltas for
        `slab` and `rpars`.
    auxbeams : str
        The contents of the AUXBEAMS file.
    phaseshifts
        The contents of the PHASESHIFTS file.

    Raises
    ------
    Exception
        If AUXBEAMS is not present and writing it anew fails.
    OSError
        If reading AUXBEAMS of PHASESHIFTS fails.
    """
    delta_basic = generateDeltaBasic(slab, rpars)
    auxbeams = _read_file_with_newline(_fetch_auxbeams(rpars))
    phaseshifts = _read_file_with_newline('PHASESHIFTS')
    return delta_basic, auxbeams, phaseshifts


def checkDelta(filename, at, el, rp):
    """Checks whether a given delta file corresponds to the requested
    displacements of a given atom. Returns True or False."""
    eps = 1e-4
    fgeo = []                               # found geo disp
    fvib = []                               # found vib disp
    at.mergeDisp(el)
    if el == "vac":
        dgeo = [np.array([0.0, 0.0, 0.0])]  # requested geo disp
        dvib = [0.0]                        # requested vib disp
    else:
        if el in at.disp_geo:
            dgeo = at.disp_geo[el]
        else:
            dgeo = at.disp_geo["all"]
        if el in at.disp_vib:
            dvib = at.disp_vib[el]
        else:
            dvib = at.disp_vib["all"]
    try:
        with open(filename, "r") as rf:
            lines = rf.readlines()
    except FileNotFoundError:
        logger.error("Error reading file "+filename)
        raise
    if rp.TL_VERSION < Version('1.7.0'):  # formatting changed to I5 in version 1.7
        intlen = 3
    else:
        intlen = 5
    try:
        nbeams = int(lines[1][0:intlen])        # no. of beams
        nvar = int(lines[1][2*intlen:3*intlen]) # no. of variations (geo*vib)
    except Exception:
        logger.error("Error parsing file " + filename)
        raise
    if nbeams != len(rp.ivbeams):
        return False
    if nvar != len(dgeo)*len(dvib):
        return False
    beams = []      # read beams from delta file
    atline = 2
    for i in range(atline, len(lines)):  # read and check beams, then break
        try:        # read hk of beams formatted as F10.5
            fl = [float(s) for s in [lines[i][j:j+10]
                                     for j in range(0, len(lines[i]) - 10,
                                                    10)]]
        except Exception:
            logger.error("Error parsing file "+filename)
            raise
        for j in range(0, int(len(fl)/2)):
            beams.append((fl[2*j], fl[(2*j) + 1]))
        if len(beams) == nbeams:
            atline = i+1
            break
    for (i, hk) in enumerate(beams):   # check beams
        if not rp.ivbeams[i].isEqual_hk(hk, eps=eps):
            return False
    atline += 1   # skip the line after beams
    # geo displacements start here
    rf74x10 = ff.FortranRecordReader("10F7.4")
    parselist = []
    repeats = False
    endgeo = False
    entrycount = 0
    for i in range(atline, len(lines)):
        try:
            fl = [f for f in rf74x10.read(lines[i]) if f is not None]
        except Exception:
            logger.error("Error parsing file "+filename)
            raise
        if len(fl) < 10:
            atline = i+1
            endgeo = True  # short line -> end of geo block
        parselist = parselist + fl
        while len(parselist) >= 3:
            v = parselist[:3]
            new = np.array([v[1], v[2], v[0]])
            parselist = parselist[3:]
            entrycount += 1
            if not repeats:
                append = True
                if fgeo and np.linalg.norm(new - fgeo[0]) < eps:
                    repeats = True
                    append = False
                if append:
                    fgeo.append(new)
                if len(fgeo) == nvar:
                    atline = i+1
                    endgeo = True
                    break
            elif not endgeo:
                if (entrycount-1) % len(fgeo) == 0:  # should repeat here
                    if np.linalg.norm(new - fgeo[0]) > eps:
                        atline = i
                        endgeo = True
                        break
        if endgeo:
            break
    ngeo = len(fgeo)
    # check geometry
    if ngeo != len(dgeo):
        return False
    for (i, gd) in enumerate(fgeo):
        if np.linalg.norm(gd - dgeo[i]) > eps:
            return False
    # vib displacement starts here
    nvib = nvar / ngeo
    if int(nvib) - nvib > 1e-4:
        logger.error("Error reading file "+filename+": number of geometry "
                     "variations found does not match header.")
        return False
    nvib = int(nvib)
    for i in range(atline, len(lines)):
        try:
            fl = [f for f in rf74x10.read(lines[i]) if f is not None]
        except Exception:
            logger.error("Error parsing file "+filename)
            raise
        parselist = parselist + fl
        while len(parselist) >= ngeo:
            fvib.append(parselist[0])
            if any([f != 0. for f in parselist[1:ngeo]]):
                logger.warning("File "+filename+": Found unexpected entries "
                               "in list of vibration displacements.")
            parselist = parselist[ngeo:]
        if len(fvib) >= nvib:
            break
    # check vibrations:
    if el.lower() == "vac":
        voff = 0.
    else:
        voff = at.site.vibamp[el]
    for (i, f) in enumerate(fvib):
        bv = round(round(dvib[i] + voff, 4) / 0.529177, 4)
        # in bohr radii; rounding twice to account for 1. writing to
        #  delta-input, 2. reading from DELTA file. Precision taken from
        #  TensErLEED GLOBAL
        if abs(f - bv) >= 1e-4:
            return False
    return True


def generateDeltaInput(atom, targetel, sl, rp, deltaBasic, auxbeams,
                       phaseshifts):
    """
    Generates a PARAM file and delta input for one element of one atom.

    Parameters
    ----------
    atom : Atom
        Atom object for which input should be generated.
    targetel : str
        The element of that atom for which input should be generated. This may
        not be the main atom element.
    sl : Slab
        The Slab object containing atom information.
    rp : Rparams
        The run parameters object.
    deltaBasic : str, optional
        Part of delta input that is the same for all atoms. Use
        generateDeltaBasic for creating this.
    auxbeams : str
        The contents of the AUXBEAMS file. Should end with a newline.
    phaseshifts : str
        The contents of the PHASESHIFTS file. Should end with a newline

    Returns
    -------
    (str, str, str).
        The delta input, a shortened version of that input for logging, and
        the contents of the required PARAM file.
    """
    MNLMB = [19, 126, 498, 1463, 3549, 7534, 14484, 25821, 43351, 69322,
             106470, 158067, 227969, 320664, 441320, 595833, 790876, 1033942]
    try:
        beamlist, _, _ = beams.writeAUXBEAMS(ivbeams=rp.ivbeams,
                                             beamlist=rp.beamlist,
                                             write=False)
    except Exception:
        logger.error("writeDeltaInput: Exception while getting data from "
                     "writeAUXBEAMS")
        raise

    # merge offsets with displacement lists
    atom.mergeDisp(targetel)

    # generate delta.in
    din = ("""   1                         FORMOUT - 1: formatted output
-------------------------------------------------------------------
--- chemical nature of displaced atom                           ---
-------------------------------------------------------------------
""")
    if targetel.lower() == "vac":
        iel = 0
    else:
        # find number of target element
        i = 0
        for el in sl.elements:
            # this reproduces the order of blocks contained in PHASESHIFTS:
            if el in rp.ELEMENT_MIX:
                chemelList = rp.ELEMENT_MIX[el]
            else:
                chemelList = [el]
            siteList = [s for s in sl.sitelist if s.el == el]
            for cel in chemelList:
                for s in siteList:
                    i += 1
                    if s.isEquivalent(atom.site) and cel == targetel:
                        iel = i
    i4 = ff.FortranRecordWriter("I4")
    f74x3 = ff.FortranRecordWriter('3F7.4')
    ol = i4.write([iel])
    din += ol.ljust(29) + "IEL  - element in PHASESHIFT list"
    if rp.TL_VERSION < Version('1.7.0'):
        din += """
-------------------------------------------------------------------
--- undisplaced position of atomic site in question             ---
-------------------------------------------------------------------
"""
        ol = f74x3.write([0.0, 0.0, 0.0])
        din += ol.ljust(29) + "CUNDISP - displacement offset"
    din += """
-------------------------------------------------------------------
--- displaced positions of atomic site in question              ---
-------------------------------------------------------------------
"""
    if targetel == "vac":
        geolist = [(0., 0., 0.)]
    elif targetel in atom.disp_geo:
        geolist = atom.disp_geo[targetel]
    else:
        geolist = atom.disp_geo["all"]
    geosteps = len(geolist)
    ol = i4.write([geosteps])
    din += ol.ljust(29) + "NCSTEP - number of displaced positions\n"
    for disp in geolist:
        ol = f74x3.write([disp[2], disp[0], disp[1]])
        din += ol.ljust(29)+"CDISP(z,x,y) - z pointing towards bulk\n"
        # TODO: should we allow this if e.g. HALTING = 1 and just warn instead?
        if any([abs(d) > 1. for d in disp]):
            logger.error(
                "Displacements for delta amplitudes have to be smaller than "
                "one Angstrom! Larger displacements are not reasonable "
                "within the tensor LEED approximation.\n"
                "Found displacement {} for {}.".format(disp, atom))
            raise ValueError("Excessive displacements (>1A) detected.")
    din += (
        """-------------------------------------------------------------------
--- vibrational displacements of atomic site in question        ---
-------------------------------------------------------------------
""")
    if targetel == "vac":
        viblist = [0.]
    elif targetel in atom.disp_vib:
        viblist = atom.disp_vib[targetel]
    else:
        viblist = atom.disp_vib["all"]
    vibsteps = len(viblist)
    ol = i4.write([vibsteps])
    din += ol.ljust(29) + "NDEB - number of vib. amplitudes to be considered\n"
    f74 = ff.FortranRecordWriter("F7.4")
    if targetel.lower() != "vac":
        # "default" vibamp + offset, not just offset
        vibamps = [v + atom.site.vibamp[targetel] for v in viblist]
        if any([v <= 0 for v in vibamps]):
            logger.warning(
                "Vibration amplitudes for {} contain values <= 0 "
                "(smallest: {:.4f}). Shifting displacement list to avoid "
                "non-positive numbers.".format(atom, min(vibamps)))
            corr = min([v for v in vibamps if v > 0]) - min(vibamps)
            vibamps = [v + corr for v in vibamps]
            for i in range(len(viblist)):
                # can't be done by list comprehension because it should modify
                #   the list that 'viblist' is pointing to, not make a copy
                viblist[i] += corr
    else:
        vibamps = [0.]
    for vibamp in vibamps:
        ol = f74.write([vibamp])
        din += ol.ljust(29)+"DRPER_A\n"

    din_main = deltaBasic + auxbeams + phaseshifts + din
    din_short = deltaBasic + "[AUXBEAMS]\n" + "[PHASESHIFTS]\n" + din

    # write PARAM
    param = ("""C  Parameter statements for delta amplitude calculation, v1.2
C  parameters must be consistent with preceding reference calculation!

C  MLMAX: maximum angular momentum to be considered in calculation
C  MNLMB: number of Clebsh-Gordon coefficients needed in tmatrix() subroutine -
C         set according to current LMAX
"""
             "C         MLMAX:  1  2   3    4    5    6    7     8     9     "
             "10    11     12     13     14     15     16     17     18\n"
             "C         MNLMB: 19 126 498 1463 3549 7534 14484 25821 43351 "
             "69322 106470 158067 227969 320664 441320 595833 790876 1033942\n" """
C  MNPSI: number of phase shift values tabulated in phase shift file
C  MNEL : number of elements for which phase shifts are tabulated
C  MNT0 : number of beams for which delta amplitude calculation is required
C  MNATOMS: currently must be set to 1. In principle number of different atomic
C      positions in a superlattice wrt the reference periodicity when computing
C      TLEED beams for a superlattice not present in the reference structure
C  MNDEB: number of thermal variation steps to be performed (outer var. loop)
C  MNCSTEP: number of geometric variation steps to be performed """
             + "(inner var. loop)\n\n")
    param += "      PARAMETER( MLMAX = {} )\n".format(rp.LMAX.max)
    param += "      PARAMETER( MNLMB = {} )\n".format(MNLMB[rp.LMAX.max-1])
    param += ("      PARAMETER( MNPSI = {}, MNEL = {} )\n"
              .format(len(rp.phaseshifts), (len(rp.phaseshifts[0][1]))))
    param += "      PARAMETER( MNT0 = {} )\n".format(len(beamlist))
    param += "      PARAMETER( MNATOMS = 1 )\n"
    param += "      PARAMETER( MNDEB = {} )\n".format(vibsteps)
    param += "      PARAMETER( MNCSTEP = {} )\n".format(geosteps)
    return din_main, din_short, param


def generateDeltaBasic(sl, rp):
    """Generates the part of the input for delta-amplitudes that is the same
    for all atoms, and returns it as a string."""
    if rp.TL_VERSION < Version('1.7.0'):
        formatter = {'energies': ff.FortranRecordWriter('3F7.2'),
                     'uc': ff.FortranRecordWriter('2F7.4'),
                     'angles': ff.FortranRecordWriter('2F7.4'),
                     }
        lj = 24  # ljust spacing
    else:
        formatter = {'energies': ff.FortranRecordWriter('2F9.2'),
                     'uc': ff.FortranRecordWriter('2F9.4'),
                     'angles': ff.FortranRecordWriter('2F9.4'),
                     'int': ff.FortranRecordWriter('I4'),
                     }
        lj = 30  # ljust spacing
    output = ""
    output += rp.systemName+" "+rp.timestamp+"\n"
    output += (formatter['energies'].write(
        [rp.THEO_ENERGIES.start, rp.THEO_ENERGIES.stop+0.01]).ljust(lj)
        + 'EI,EF\n')
    ucsurf = sl.ab_cell.T
    if sl.bulkslab is None:
        sl.make_bulk_slab(rp)
    ucbulk = sl.bulkslab.ab_cell.T
    output += formatter['uc'].write(ucbulk[0]).ljust(lj) + 'ARA1\n'
    output += formatter['uc'].write(ucbulk[1]).ljust(lj) + 'ARA2\n'
    output += formatter['uc'].write(ucsurf[0]).ljust(lj) + 'ARB1\n'
    output += formatter['uc'].write(ucsurf[1]).ljust(lj) + 'ARB2\n'
    output += (formatter['angles'].write([rp.THETA, rp.PHI]).ljust(lj)
               + 'THETA PHI\n')
    if rp.TL_VERSION >= Version('1.7.0'):
        # TODO: if phaseshifts are calculated differently, change format here
        output += (formatter['int'].write([1]).ljust(lj)
                   + 'PSFORMAT  1: Rundgren_v1.6; 2: Rundgren_v1.7\n')
    return output


def write_delta_input_file(compile_tasks, run_tasks):
    """Write a collection of the inputs for all delta calculations.

    The delta-input file is meant for users' debug purposes (or
    manual execution of a delta calculation). It collated the
    PARAM files (i.e., array dimensions) for the delta-amplitude
    executables that are compiled, as well as a short version of
    the input piped to these executables when called to produce
    single delta files.

    Parameters
    ----------
    compile_tasks : Sequence of DeltaCompileTask
        Information about which executables need to be compiled.
    run_tasks : Sequence of DeltaRunTask
        Information about which delta calculations should be
        executed.

    Returns
    -------
    None.
    """
    fpath = Path('delta-input')
    dinput = '''\
# ABOUT THIS FILE:
# Input for the delta-calculations is collected here. The blocks of data are
# new 'PARAM' files, which are used to recompile the fortran code, and input
# for generation of specific DELTA files. Lines starting with '#' are comments
# on the function of the next block of data.
# In the DELTA file blocks, [AUXBEAMS] and [PHASESHIFTS] denote where the
# entire contents of the AUXBEAMS and PHASESHIFTS files should be inserted.
'''
    for compile_task in compile_tasks:
        dinput += f'''
#### NEW 'PARAM' FILE: ####

{compile_task.param}
'''
        for run_task in run_tasks:
            if run_task.comptask is not compile_task:
                continue
            dinput += f'''
#### INPUT for new DELTA file {run_task.deltaname}: ####

{run_task.din_short}
'''
    try:
        fpath.write_text(dinput, encoding='utf-8')
    except OSError:
        logger.warning(f'Failed to write file {fpath.name!r}. This will '
                       'not affect TensErLEED execution, proceeding...')


def _fetch_auxbeams(rpars):
    """Collect an existing AUXBEAMS file, or write a new one."""
    auxbeams = Path('AUXBEAMS')
    if not auxbeams.is_file():  # Try fetching it from SUPP
        try:
            shutil.copy2(DEFAULT_SUPP/auxbeams, auxbeams)
        except FileNotFoundError:
            pass  # Will write a new one below
        except OSError:
            logger.warning(f'Failed to copy {auxbeams} from {DEFAULT_SUPP} '
                           'folder. Generating new file...')
    if not auxbeams.is_file():
        try:
            beams.writeAUXBEAMS(ivbeams=rpars.ivbeams, beamlist=rpars.beamlist)
        except Exception:                                                       # TODO: better exception
            logger.error('Exception during writeAUXBEAMS: ')
            raise
    return auxbeams


def _read_file_with_newline(file_path):
    """Return the contents of `file_path`, with a terminating newline."""
    file_path = Path(file_path)
    try:
        contents = file_path.read_text(encoding='utf-8')
    except OSError:
        logger.error(f'Could not read {file_path.name} for delta input')
        raise
    if not contents.endswith('\n'):
        contents += '\n'
    return contents
