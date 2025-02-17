"""Module iorefcalc of viperleed.calc.files.

Defines functions for reading and writing files relevant to the
reference calculation.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-19'
__license__ = 'GPLv3+'

import copy
import logging
import os

import fortranformat as ff
import numpy as np

from viperleed.calc.classes.beam import Beam
from viperleed.calc.files.beams import writeAUXBEAMS
from viperleed.calc.lib import leedbase
from viperleed.calc.lib.fortran_utils import wrap_fortran_line
from viperleed.calc.lib.string_utils import rsplit_once
from viperleed.calc.lib.version import Version

logger = logging.getLogger(__name__)


def combine_tensors(oripath=".", targetpath=".", buffer=0):
    """Combines Tensor files in 'oripath' into common files in 'targetpath'.
    The 'buffer' argument specified how much of a Tensor file can be stored in
    memory at a given time (in bytes). Set 0 to read the entire file at once
    (faster, but may be a problem if memory is limited)."""
    tensorfiles = [f for f in os.listdir(oripath) if f.startswith("T_")
                   and os.path.isfile(os.path.join(oripath, f))
                   and len(f.split("_")) == 3]
    tnum = {}
    for f in tensorfiles:
        try:
            num = int(f.split("_")[1])
            float(f.split("_")[2][:-2])
        except ValueError:
            continue   # not a real Tensor file
        if num not in tnum:
            tnum[num] = [f]
        else:
            tnum[num].append(f)
    if len(tnum) == 0:
        logger.error("Found no Tensor files to combine.")
        raise RuntimeError("No Tensor files found")
    for num in tnum:
        tlist = sorted(tnum[num], key=lambda x: float(x.split("_")[2][:-2]))
        name = os.path.join(targetpath, rsplit_once(tlist[0], "_")[0])
        with open(name, "w") as wf:
            # having the 'if' this far outside duplicates code, but speeds
            #  up the loop
            if buffer == 0:
                for tf in tlist:
                    with open(os.path.join(oripath, tf), "r") as rf:
                        wf.write(rf.read())
                    try:
                        os.remove(os.path.join(oripath, tf))
                    except Exception:
                        logger.warning("Failed to delete file " + tf)
            else:
                for tf in tlist:
                    with open(os.path.join(oripath, tf), "r") as rf:
                        while True:
                            data = rf.read(buffer)
                            if data:
                                wf.write(data)
                            else:
                                break
                    try:
                        os.remove(os.path.join(oripath, tf))
                    except Exception:
                        logger.warning("Failed to delete file " + tf)
    return


def combine_fdout(oripath=".", targetpath="."):
    """Combines fd.out and amp.out files in oripath into a common file at
    targetpath."""
    for filetype in ("fd", "amp"):
        outlines = []
        startstr = filetype + "_"
        filelist = [f for f in os.listdir(oripath)
                    if os.path.isfile(os.path.join(oripath, f))
                    and f.startswith(startstr) and f.endswith("eV.out")]
        i = 0
        while i < len(filelist):
            try:
                float(filelist[i].split(startstr)[1].split("eV.out")[0])
                i += 1
            except ValueError:
                filelist.pop(i)
        if len(filelist) == 0:
            if filetype == "amp":
                return   # return without error, we don't need amp files
            logger.error("Found no fd.out files to combine.")
            raise RuntimeError("No fd.out files found")
        filelist.sort(key=lambda x: float(x.split(startstr)[1]
                                          .split("eV.out")[0]))
        nbeams = 0
        for f in filelist:
            with open(os.path.join(oripath, f), "r") as rf:
                lines = rf.readlines()
            lines = [s for s in lines if ".  CORRECT TERMINATION" not in s
                     and len(s.strip()) >= 1]
            if len(outlines) == 0:
                try:
                    nbeams = int(lines[1].strip().split()[0])
                except Exception:
                    logger.error("Failed to combine {}.out".format(filetype),
                                 exc_info=True)
                    raise
                outlines += lines
            else:
                outlines += lines[nbeams+2:]
            if (nbeams % 5) == 4:
                outlines += "\n"  # expects an empty line (fortran formatting)
            try:
                os.remove(os.path.join(oripath, f))
            except Exception:
                logger.warning("Failed to delete file " + f)
        try:
            with open(os.path.join(targetpath, (filetype + ".out")),
                      "w") as wf:
                wf.write("".join(outlines))
        except Exception:
            logger.error("Failed to write combined {}.out".format(filetype))
            raise
    return


def readFdOut(readfile="fd.out", for_error=False, ampfile="amp.out"):
    """Reads the fd.out file produced by the refcalc and returns a list of
    Beam objects. If 'ampfile' is set, will attempt to read complex amplitudes
    from the given file as well."""

    def parse_data(lines, which, filename, theobeams, nbeams):
        out_str = ""
        blocks = []
        for line in lines:
            if ".  CORRECT TERMINATION" in line:
                # sometimes superpos output is not in right line. sync issue?
                continue
            out_str += line + "\n"
            llist = line.split()
            if len(llist) == 0:
                continue  # skip empty lines
            try:
                float(llist[0])
            except ValueError:
                break  # end of data
            blocks.extend(llist)

        # Each block contains exactly nbeams+2 entries (double for amp)
        # reshape blocks so that each line corresponds to one block
        if which == "fd":
            blocks = np.reshape(blocks, (-1, nbeams+2))
        else:
            blocks = np.reshape(blocks, (-1, 2*nbeams+2))

        # and parse the data
        warned = False
        for block in blocks:
            if block[1] != "0.0001":
                if not for_error and not warned:
                    warned = True
                    logger.warning("File " + filename + "contains unexpected "
                                   "data for more than one structural "
                                   "variation. Only the first block was read.")
                continue
            en = float(block[0])
            values = [float(s) for s in block[2:]]
            for (j, beam) in enumerate(theobeams):
                if which == "fd":
                    beam.intens[en] = values[j]
                else:
                    beam.complex_amplitude[en] = complex(values[2*j],
                                                         values[2*j + 1])
        return out_str

    try:
        with open(readfile, 'r') as rf:
            filelines = [line[:-1] for line in rf.readlines()]
            rf.seek(0)
    except FileNotFoundError:
        logger.error("Error in readFdOut: file "+str(readfile)+" not found.")
        raise
    if filelines[0].startswith(" SOME ERROR OCCURED WHILE READING"):
        logger.error("File "+readfile+" reports error: "+filelines[0])
        raise Exception("File "+readfile+" reports error.")
    fdout = ""
    theobeams = []
    i = 1   # number lines as in text editor - be careful about indexing!
    nbeams = 1e10   # some large number, just to enter the while loop
    while i < nbeams+3:
        fdout += filelines[i-1] + "\n"
        llist = filelines[i-1].split()
        if i == 1:
            pass  # header - skip
        elif i == 2:
            nbeams = int(llist[0])
        else:
            theobeams.append(Beam((float(llist[1]), float(llist[2]))))
        i += 1

    # re-label the beams to get the correct number of characters and formatting
    mw = max([beam.getLabel()[1] for beam in theobeams])
    for beam in theobeams:
        beam.label = beam.getLabel(lwidth=mw)[0]

    fdout += parse_data(filelines[i-1:], "fd", readfile, theobeams ,nbeams)

    if ampfile and os.path.isfile(ampfile):
        with open(ampfile, 'r') as rf:
            amplines = [line[:-1] for line in rf.readlines()]
        # check header
        if any(amplines[i] != filelines[i] for i in range(1, nbeams + 2)):
            logger.warning("Failed to read " + ampfile + ": Header does not "
                           "match " + readfile)
            return theobeams, fdout
        # now read the rest
        parse_data(amplines[nbeams+2:], "amp", ampfile, theobeams, nbeams)
    if not theobeams:
        raise RuntimeError(f"No beams found in {readfile}")
    return theobeams, fdout


def writePARAM(sl, rp, lmax=-1):
    """Creats the contents of the PARAM file for the reference calculation.
    If no LMAX is passed, will use maximum LMAX from rp. Returns str."""
    if lmax == -1:
        lmax = rp.LMAX.max
    try:
        beamlist, beamblocks, beamN = writeAUXBEAMS(
            ivbeams=rp.ivbeams, beamlist=rp.beamlist, write=False)
    except Exception:
        logger.error("generatePARAM: Exception while getting data from "
                     "writeAUXBEAMS")
        raise

    # define Clebsh-Gordon coefficient tables:
    mnlmo = [1, 70, 264, 759, 1820, 3836, 7344, 13053, 21868, 34914, 53560,
             79443, 114492, 160952, 221408, 298809, 396492, 518206]
    mnlm = [1, 76, 284, 809, 1925, 4032, 7680, 13593, 22693, 36124, 55276,
            81809, 117677, 165152, 226848, 305745, 405213, 529036]

    # start generating output
    output = ('C  Dimension statements for Tensor LEED reference calculation, '
              '\nC  version v1.2\n\n')
    output += 'C  1. lattice symmetry\n\n'
    nl1, nl2 = leedbase.get_superlattice_repetitions(rp.SUPERLATTICE)
    ideg = 2  # any 2D point grid is at least 2fold symmetric
    # if sl.planegroup in ['p2','pmm','pmg','pgg','cmm','rcmm']:
    #     ideg = 2
    if sl.planegroup in ['p3', 'p3m1', 'p31m']:
        ideg = 3
    elif sl.planegroup in ['p4', 'p4m', 'p4g']:
        ideg = 4
    elif sl.planegroup in ['p6', 'p6m']:
        ideg = 3
        # should be 6, but according to TensErLEED fortran comments,
        #   3 works better
    output += ('      PARAMETER (MIDEG='+str(ideg)+',MNL1='+str(nl1)
               + ',MNL2='+str(nl2)+')\n')
    output += '      PARAMETER (MNL = MNL1*MNL2)\n'
    output += '\nC  2. General calculational quantities\n\n'
    output += '      PARAMETER (MKNBS = '+str(beamblocks)+')\n'
    output += '      PARAMETER (MKNT =  '+str(beamN)+')\n'
    output += ('      PARAMETER (MNPUN = '+str(len(beamlist))+', MNT0 = '
               + str(len(beamlist))+')\n')
    output += ('      PARAMETER (MNPSI = '+str(len(rp.phaseshifts))+', MNEL = '
               + str(len(rp.phaseshifts[0][1]))+')\n')
    output += '      PARAMETER (MLMAX = '+str(lmax)+')\n'
    output += ('      PARAMETER (MNLMO = '+str(mnlmo[lmax-1])+', MNLM = '
               + str(mnlm[lmax-1])+')\n')
    output += '\nC  3. Parameters for (3D) geometry within (2D) unit mesh\n\n'
    output += '      PARAMETER (MNSITE  = '+str(len(sl.sitelist))+')\n'
    output += '      PARAMETER (MNLTYPE = '+str(sl.n_layers)+')\n'
    mnbrav = 0
    mnsub = 0
    mnstack = 0
    if sl.bulkslab is None:
        sl.make_bulk_slab(rp)
    for layer in sl.non_bulk_layers:
        mnstack += 1
        if layer.n_atoms == 1:
            mnbrav += 1
        if layer.n_atoms > mnsub:
            mnsub = layer.n_atoms
    for i, layer in enumerate(sl.bulk_layers):
        if sl.bulkslab.layers[i].n_atoms == 1:
            mnbrav += 1
        if sl.bulkslab.layers[i].n_atoms > mnsub:
            mnsub = layer.n_atoms
    output += '      PARAMETER (MNBRAV  = '+str(mnbrav)+')\n'
    output += '      PARAMETER (MNSUB   = '+str(mnsub)+')\n'
    output += '      PARAMETER (MNSTACK = '+str(mnstack)+')\n'
    output += ('\nC  4. some derived quantities that must be treated '
               'explicitly (dummy dimensions for\n')
    output += 'C     special cases necessary\n\n'
    output += '      PARAMETER (MLMAX1=MLMAX+1)\n'
    output += '      PARAMETER (MLMMAX = MLMAX1*MLMAX1)\n\n'
    output += ('      PARAMETER (MNBRAV2 = '+('MNBRAV' if mnbrav > 0
                                              else '1')+')\n\n')
    output += ('      PARAMETER (MNCOMP= '+('MNLTYPE-MNBRAV' if mnsub > 1
                                            else '1 ')+')\n')
    output += ('      PARAMETER (MLMT  = '+('MNSUB*MLMMAX' if mnsub > 1
                                            else '1 ')+')\n')
    output += ('      PARAMETER (MNSUB2= '+('MNSUB * (MNSUB-1)/2' if mnsub > 1
                                            else '1 ')+')\n')
    output += ('      PARAMETER (MLMG  = '+('MNSUB2*MLMMAX*2' if mnsub > 1
                                            else '1 ')+')\n')
    output += ('      PARAMETER (MLMN  = '+('MNSUB * MLMMAX' if mnsub > 1
                                            else '1 ')+')\n')
    output += ('      PARAMETER (MLM2N = '+('2*MLMN' if mnsub > 1
                                            else '1 ')+')\n')
    output += ('      PARAMETER (MLMNI = '+('MNSUB*MLMMAX' if mnsub > 1
                                            else '1 ')+')\n')
    return output


def collectFIN(version=0.):
    """Combines AUXLATGEO, BEAMLIST, AUXNONSTRUCT, PHASESHIFTS, AUXBEAMS
    and AUXGEO into one string (input for refcalc), which it returns. Pass
    beamlist to avoid reading it again."""
    if version < Version('1.7.3'):
        filenames = ["AUXLATGEO", "BEAMLIST", "AUXNONSTRUCT", "PHASESHIFTS",
                     "AUXBEAMS", "AUXGEO"]
    else:
        filenames = ["AUXLATGEO", "BEAMLIST", "AUXNONSTRUCT", "PHASESHIFTS",
                     "AUXGEO"]
    fin = ""
    for fn in filenames:
        with open(fn, "r") as rf:
            fin += rf.read()
        if fin[-1] != "\n":
            fin += "\n"
    return fin


def writeAUXLATGEO(sl, rp):
    """Writes AUXLATGEO, which is part of the input FIN for the refcalc."""
    if rp.TL_VERSION < Version('1.7.0'):
        formatter = {'energies': ff.FortranRecordWriter('3F7.2'),
                     'uc': ff.FortranRecordWriter('2F7.4'),
                     'x_ase_wf': ff.FortranRecordWriter('2F7.4'),
                     }
        lj = 24  # ljust spacing
    else:
        formatter = {'energies': ff.FortranRecordWriter('3F9.2'),
                     'uc': ff.FortranRecordWriter('2F9.4'),
                     'x_ase_wf': ff.FortranRecordWriter('3F9.2'),
                     }
        lj = 30  # ljust spacing
    output = ''
    output += rp.systemName+' '+rp.timestamp+'\n'
    ens = [rp.THEO_ENERGIES.start,
           rp.THEO_ENERGIES.stop + 0.01,
           rp.THEO_ENERGIES.step]
    output += formatter['energies'].write(ens).ljust(lj) + 'EI,EF,DE\n'
    ucsurf = sl.ab_cell.T
    if sl.bulkslab is None:
        sl.make_bulk_slab(rp)
    ucbulk = sl.bulkslab.ab_cell.T
    output += formatter['uc'].write(ucbulk[0]).ljust(lj) + 'ARA1\n'
    output += formatter['uc'].write(ucbulk[1]).ljust(lj) + 'ARA2\n'
    if rp.TL_VERSION < Version('1.7.0'):
        output += ' 0.0    0.0             SS1\n'
        output += ' 0.0    0.0             SS2\n'
        output += ' 0.0    0.0             SS3\n'
        output += ' 0.0    0.0             SS4\n'
    output += formatter['uc'].write(ucsurf[0]).ljust(lj) + 'ARB1\n'
    output += formatter['uc'].write(ucsurf[1]).ljust(lj) + 'ARB2\n'
    if rp.TL_VERSION < Version('1.7.0'):
        output += ' 0.0    0.0             SO1\n'
        output += ' 0.0    0.0             SO2\n'
        output += ' 0.0    0.0             SO3\n'
        output += (formatter['x_ase_wf'].write([0.5, rp.V0_Z_ONSET]).ljust(lj)
                   + 'FR ASE\n')
    else:
        # Different parameters here in version 1.7! Previously in muftin
        output += (formatter['x_ase_wf'].write([rp.V0_IMAG, rp.V0_Z_ONSET,
                                                rp.FILAMENT_WF]).ljust(lj)
                   + 'V0i,ASE,WORKFN\n')
    try:
        with open('AUXLATGEO', 'w') as wf:
            wf.write(output)
    except Exception:
        logger.error("Failed to write AUXLATGEO file")
        raise
    logger.debug("Wrote to AUXLATGEO successfully.")
    return


def writeAUXNONSTRUCT(sl, rp):
    """Writes AUXNONSTRUCT, which is part of the input FIN for the refcalc."""
    try:
        beamnums, _, _ = writeAUXBEAMS(ivbeams=rp.ivbeams,
                                       beamlist=rp.beamlist,
                                       write=False)
    except Exception:
        logger.error("generatePARAM: Exception while getting data from "
                     "writeAUXBEAMS")
        raise
    if rp.TL_VERSION < Version('1.7.0'):
        formatter = {'tst': ff.FortranRecordWriter('F7.4'),
                     'incidence': ff.FortranRecordWriter('(F7.2, F6.1)'),
                     'beamnums': ff.FortranRecordWriter('15I4'),
                     'eps': ff.FortranRecordWriter('F7.4'),
                     'ints': ff.FortranRecordWriter('I3'),
                     }
    else:
        formatter = {'tst': ff.FortranRecordWriter('F11.8'),
                     'incidence': ff.FortranRecordWriter('2F9.4'),
                     'beamnums': ff.FortranRecordWriter('15I5'),
                     'eps': ff.FortranRecordWriter('F7.4'),
                     'ints': ff.FortranRecordWriter('I3'),
                     }
    output = ''

    output += (formatter['tst'].write([rp.ATTENUATION_EPS]).ljust(18)
               + '>>>>> ! <<<<<              TST\n')
    output += formatter['beamnums'].write(beamnums)+'\n'
    output += (formatter['incidence'].write([rp.THETA, rp.PHI]).ljust(45)
               + 'THETA FI\n')
    output += formatter['eps'].write([rp.BULKDOUBLING_EPS]).ljust(45) + 'EPS\n'
    output += (formatter['ints'].write([rp.BULKDOUBLING_MAX]).ljust(45)
               + 'LITER\n')
    output += formatter['ints'].write([rp.LMAX.max]).ljust(45) + 'LMAX\n'
    if rp.TL_VERSION >= Version('1.7.0'):
        # TODO: if phaseshifts are calculated differently, change format here
        output += (formatter['ints'].write([1]).ljust(45)
                   + 'PSFORMAT  1: Rundgren_v1.6; 2: Rundgren_v1.7\n')
    if rp.TL_VERSION >= Version('1.7.3'):
        output += (formatter['ints'].write([1]).ljust(45)
                   + 'IFORM - formatted input and output\n')
    try:
        with open('AUXNONSTRUCT', 'w') as wf:
            wf.write(output)
    except Exception:
        logger.error("Failed to write AUXNONSTRUCT file")
        raise
    logger.debug("Wrote to AUXNONSTRUCT successfully.")
    return


def writeAUXGEO(sl, rp):
    """Writes AUXGEO, which is part of the input FIN for the refcalc."""
    if rp.TL_VERSION < Version('1.7.0'):
        formatter = {'vibrocc': ff.FortranRecordWriter('2F7.4'),
                     'geo': ff.FortranRecordWriter('3F7.4'),
                     }
        lj = 26
    else:
        formatter = {'vibrocc': ff.FortranRecordWriter('2F9.4'),
                     'geo': ff.FortranRecordWriter('3F9.4'),
                     }
        lj = 32
    slab_c = sl.c_vector.copy()
    if rp.LAYER_STACK_VERTICAL:
        sl = copy.deepcopy(sl)
        sl.project_c_to_z()
        sl.update_layer_coordinates()
    output = ''
    output += ('---------------------------------------------------------'
               '----------\n')
    output += ('--- define chem. and vib. properties for different atomic '
               'sites ---\n')
    output += ('---------------------------------------------------------'
               '----------\n')
    i3 = ff.FortranRecordWriter('I3')
    ol = i3.write([len(sl.sitelist)]).ljust(lj)
    output += ol + 'NSITE: number of different site types\n'
    for i, site in enumerate(sl.sitelist):
        output += '-   site type '+str(i+1)+' ---\n'

        for el in sl.elements:
            # this reproduces the order of blocks contained in PHASESHIFTS:
            if el in rp.ELEMENT_MIX:
                chemelList = rp.ELEMENT_MIX[el]
            else:
                chemelList = [el]
            siteList = [s for s in sl.sitelist if s.el == el]
            for cel in chemelList:
                for s in siteList:
                    if s.isEquivalent(site):
                        occ, vib = site.occ[cel], site.vibamp[cel]
                        comment = ('Occ & VibAmp for '+cel+' in '+site.label
                                   + ' site')
                    else:
                        occ, vib = 0., 0.
                        comment = ''
                    try:
                        ol = formatter['vibrocc'].write([occ, vib]).ljust(lj)
                    except Exception:
                        logger.error(
                            "Exception while trying to write occupation / "
                            "vibration amplitude for site " + site.label,
                            exc_info=True)
                    output += ol + comment + '\n'

    output += ('-----------------------------------------------------'
               '--------------\n')
    output += ('--- define different layer types                     '
               '           ---\n')
    output += ('-----------------------------------------------------'
               '--------------\n')
    ol = i3.write([sl.n_layers]).ljust(lj)
    output += ol + 'NLTYPE: number of different layer types\n'
    blayers = sl.bulk_layers
    layerOffsets = [np.zeros(3) for _ in range(sl.n_layers + 1)]
    if sl.bulkslab is None:
        sl.make_bulk_slab(rp)
    for i, layer in enumerate(sl.layers):
        output += '-   layer type '+str(i+1)+' ---\n'
        if layer.is_bulk:
            output += ('  2'.ljust(lj) + 'LAY = 2: layer type no. '
                       + str(i+1) + ' has bulk lateral periodicity\n')
        else:
            output += ('  1'.ljust(lj) + 'LAY = 1: layer type no. '
                       + str(i+1) + ' has overlayer lateral periodicity\n')
        if layer.is_bulk:
            bl = sl.bulkslab.layers[blayers.index(layer)]
            bulknums = {at.num for at in bl}
            bulkUnique = [at for at in layer if at.num in bulknums]
            natoms = len(bulkUnique)
            # sanity check: ratio of unit cell areas (given simply by
            #  SUPERLATTICE) should match ratio of written vs skipped atoms:
            arearatio = 1 / abs(np.linalg.det(rp.SUPERLATTICE))
            atomratio = len(bulkUnique) / layer.n_atoms
            if abs(arearatio - atomratio) > 1e-3:
                logger.warning(
                    'Ratio of bulk atoms inside/outside the bulk unit cell '
                    'differs from bulk/slab unit cell size ratio. This means '
                    'that the actual periodicity of the POSCAR does not match '
                    'the periodicity given in the SUPERLATTICE parameter. '
                    'Check SUPERLATTICE parameter and bulk symmetry!')
                rp.setHaltingLevel(2)
        else:
            natoms = layer.n_atoms
        ol = i3.write([natoms]).ljust(lj)
        output += ol+'number of Bravais sublayers in layer '+str(i+1)+'\n'
        if layer.is_bulk:
            writelist = bulkUnique
        else:
            writelist = layer.atlist
        writelist.sort(key=lambda atom: -atom.pos[2])
        for atom in writelist:
            writepos = atom.cartpos - atom.layer.cartori                        # TODO: .cartpos[2]. Issue #174
            ol = i3.write([sl.sitelist.index(atom.site)+1])
            if natoms != 1:
                ol += formatter['geo'].write([writepos[2],
                                              writepos[0], writepos[1]])
            else:
                # Bravais layers need to have coordinate (0., 0., 0.)
                #  -> store actual position for later, it will go into the
                #  interlayer vector
                ol += formatter['geo'].write([0., 0., 0.])
                layerOffsets[layer.num] += writepos
                layerOffsets[layer.num + 1] -= writepos
            ol = ol.ljust(lj)
            output += f'{ol}Atom N={atom.num} ({atom.el})\n'
    output += ('--------------------------------------------------------------'
               '-----\n')
    output += ('--- define bulk stacking sequence                             '
               '  ---\n')
    output += ('--------------------------------------------------------------'
               '-----\n')
    output += ('  0'.ljust(lj) + 'TSLAB = 0: compute bulk using layer '
               'doubling\n')

    # determine ASA
    if type(rp.BULK_REPEAT) == np.ndarray:
        bvectors_ASA = np.copy(rp.BULK_REPEAT) * np.array([1, 1, -1])
        if bvectors_ASA[2] < 0:
            bvectors_ASA = -bvectors_ASA
        bulkc = bvectors_ASA[2]
    else:
        bulkc = np.copy(rp.BULK_REPEAT)
        bvectors_ASA = -slab_c * bulkc/slab_c[2]
    # bulkc is now the repeat length along c. now correct for layer thickness:
    if rp.N_BULK_LAYERS == 2:
        bvectors_ASA[2] = bulkc - (blayers[1].cartbotz - blayers[0].cartori[2])  # TODO: .cartpos[2]. Issue #174
    else:
        bvectors_ASA[2] = bulkc - (blayers[0].cartbotz - blayers[0].cartori[2])  # TODO: .cartpos[2]. Issue #174

    # determine ASBULK - interlayer vector between bulk layers
    if rp.N_BULK_LAYERS == 2:
        # add layerOffsets for Bravais layers:
        bvectors_ASA += layerOffsets[blayers[0].num+1]
        # calculate ASBULK:
        bvectors_ASBULK = blayers[1].cartori - blayers[0].cartori               # TODO: .cartpos[2]. Issue #174
        bvectors_ASBULK[2] = blayers[1].cartori[2] - blayers[0].cartbotz
        bl2num = blayers[1].num
        # add layerOffsets for Bravais layers:
        bvectors_ASBULK -= layerOffsets[blayers[0].num+1]
        # correct ASA xy in case ASBULK already includes some:
        bvectors_ASA[:2] -= bvectors_ASBULK[:2]
    else:
        bl2num = blayers[0].num
        bvectors_ASBULK = bvectors_ASA

    ol = formatter['geo'].write([bvectors_ASA[2],
                                 bvectors_ASA[0], bvectors_ASA[1]]).ljust(lj)
    output += ol + 'ASA interlayer vector between different bulk units\n'
    ol = i3.write([blayers[0].num+1]).ljust(lj)
    output += (ol + 'top layer of bulk unit: layer type '+str(blayers[0].num+1)
               + '\n')
    ol = i3.write([bl2num+1]).ljust(lj)
    output += ol + 'bottom layer of bulk unit: layer type '+str(bl2num+1)+'\n'
    ol = formatter['geo'].write([bvectors_ASBULK[2], bvectors_ASBULK[0],
                                 bvectors_ASBULK[1]]).ljust(lj)
    output += (ol + 'ASBULK between the two bulk unit layers (may differ from '
               'ASA)\n')
    output += ('--------------------------------------------------------------'
               '-----\n')
    output += ('--- define layer stacking sequence and Tensor LEED output     '
               '  ---\n')
    output += ('--------------------------------------------------------------'
               '-----\n')
    nonbulk = sl.n_layers - rp.N_BULK_LAYERS
    if len(rp.TENSOR_OUTPUT) < nonbulk:    # check TENSOR_OUTPUT parameter
        if len(rp.TENSOR_OUTPUT) == 1:
            # interpret one value as concerning all
            rp.TENSOR_OUTPUT = rp.TENSOR_OUTPUT * nonbulk
        elif rp.TENSOR_OUTPUT:
            logger.warning(
                'Parameters TENSOR_OUTPUT is defined, but contains fewer '
                'values than there are non-bulk layers. Missing values '
                'will be set to 1.')
            rp.setHaltingLevel(1)
        for i in range(0, nonbulk-len(rp.TENSOR_OUTPUT)):
            rp.TENSOR_OUTPUT.append(1)
    if len(rp.TENSOR_OUTPUT) > nonbulk:
        logger.warning(
            'Parameters TENSOR_OUTPUT is defined, but contains more values '
            'than there are non-bulk layers. Excess values will be ignored.')
        rp.setHaltingLevel(1)
    ol = i3.write([sl.n_layers - rp.N_BULK_LAYERS]).ljust(lj)
    output += ol + 'NSTACK: number of layers stacked onto bulk\n'
    for layer in reversed(sl.non_bulk_layers):
        n = layer.num + 1
        v = sl.layers[n].cartori - layer.cartori
        v[2] = sl.layers[n].cartori[2] - layer.cartbotz                         # TODO: .cartpos[2]. Issue #174
        v = v + layerOffsets[n]   # add layerOffsets for Bravais layers
        ol = i3.write([n]) + formatter['geo'].write([v[2],
                                                     v[0], v[1]])
        output += (ol.ljust(lj) + 'layer '+str(n)+': layer type '+str(n)
                   + ', interlayer vector below\n')
        ol = i3.write([rp.TENSOR_OUTPUT[layer.num]]).ljust(lj)
        output += (ol + '0/1: Tensor output is required for this layer '
                   '(TENSOR_OUTPUT)\n')
        if rp.TENSOR_OUTPUT[layer.num] == 0:
            continue   # don't write the Tensor file names
        for i, atom in enumerate(layer):
            ol = f'T_{atom.num}'.ljust(lj)
            output += (ol + 'Tensor file name, current layer, sublayer '
                       + str(i+1) + '\n')
    output += ('--------------------------------------------------------------'
               '-----\n')
    output += ('--- end geometrical input                                     '
               '  ---\n')
    output += ('--------------------------------------------------------------'
               '-----\n')

    try:
        with open('AUXGEO', 'w') as wf:
            wf.write(output)
    except Exception:
        logger.error("Failed to write AUXGEO file")
        raise
    logger.debug("Wrote to AUXGEO successfully.")
    return


def writeMuftin(rp):
    """Writes a muftin.f file, which will be compiled for the refcalc."""
    output = """
C  Subroutine muftin contains explicit energy dependence of inner
C  potentials. The functional form should be left to the user entirely,
C  thus the subroutine is included in the input explicitly.

C  All quantities should be set in eV. Conversion to Hartree (atomic units)
C  is handled by ref-calc.f itself. Both the real and imaginary part should
C  be given as positive values (the program then handles the sign correctly).

      subroutine muftin(EEV,VO,VV,VPI,VPIS,VPIO)

C  global variable

C  EEV :  electron energy in the vacuum region
C  VO  :  difference between real part of inner potential in the bulk and in
C         topmost layer. Usually irrelevant -> set to zero.
C  VV  :  real part of the inner potential in the bulk.
C  VPI :  imaginary part of the inner potential.
C  VPIS:  in case different values for the imaginary part of the inner
C  VPIO:  potential for bulk and topmost layer were desired, VPIS would
C         correspond to the bulk value, VPIO would be the respective
C         value for the topmost layer. Usually, the effect is irrelevant,
C         i.e. both are set equal to VPI.

      real EEV,VO,VV,VPI,VPIS,VPIO

C  local variable

C  workfn: Work function of the LEED filament material. Theoretical predictions
C         for the inner potential usually have their energy zero fixed at the
C         Fermi level; the experimental energy scale (-> EEV) nominally
C         does the same, except for the fact that electrons need to overcome
C         the work function of the cathode first. Note that the this value is
C         thus formally determined by a LEED fit, yet don't trust its accuracy.

      real workfn

C  set work function of cathode
C  work function should be positive (added to exp. energy EEV)

      workfn = """
    output += str(round(rp.FILAMENT_WF, 4))+"\n"
    output += """
C  set real part of inner potential

"""
    if type(rp.V0_REAL) == list:
        oline = ("      VV = workfn-max({:.2f}, (({:.2f})+({:.2f})/sqrt("
                 "EEV+workfn+({:.2f}))))".format(*rp.V0_REAL))
    else:
        oline = "      VV = "+rp.V0_REAL
    output += wrap_fortran_line(oline) + "\n"
    output += """
      write(6,*) workfn, EEV
      write(6,*) VV

c  set difference between bulk and overlayer inner potential

      VO = 0.

c  set imaginary part of inner potential - energy independent value used here

"""
    oline = "      VPI = "+str(rp.V0_IMAG)
    output += wrap_fortran_line(oline) + "\n"
    output += """
C  set substrate / overlayer imaginary part of inner potential

"""
    oline = "      VPIS = "+str(rp.V0_IMAG)
    output += wrap_fortran_line(oline) + "\n"
    oline = "      VPIO = "+str(rp.V0_IMAG)
    output += wrap_fortran_line(oline) + "\n"
    output += """
      return
      end"""
    try:
        with open("muftin.f", "w") as wf:
            wf.write(output)
        logger.debug("Wrote to muftin.f successfully.")
    except Exception:
        logger.error("Exception while writing muftin.f file: ",
                     exc_info=True)
        raise
    return
