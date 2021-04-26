# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 11:55:15 2020

@author: Florian Kraushofer

Functions for reading and writing files relevant to the reference calculation
"""

import numpy as np
import logging
import os
from viperleed import fortranformat as ff

import viperleed.tleedmlib as tl
from viperleed.tleedmlib.base import splitMaxRight
from viperleed.tleedmlib.files.beams import writeAUXBEAMS

logger = logging.getLogger("tleedm.files.iorefcalc")


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
        name = os.path.join(targetpath, splitMaxRight(tlist[0], "_")[0])
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
    """Combines fd.out files in oripath into a common file at targetpath."""
    outlines = []
    fdfiles = [f for f in os.listdir(oripath)
               if os.path.isfile(os.path.join(oripath, f))
               and f.startswith("fd_") and f.endswith("eV.out")]
    i = 0
    while i < len(fdfiles):
        try:
            float(fdfiles[i].split("fd_")[1].split("eV.out")[0])
            i += 1
        except ValueError:
            fdfiles.pop(i)
    if len(fdfiles) == 0:
        logger.error("Found no fd.out files to combine.")
        raise RuntimeError("No fd.out files found")
    fdfiles.sort(key=lambda x: float(x.split("fd_")[1].split("eV.out")[0]))
    nbeams = 0
    for f in fdfiles:
        with open(os.path.join(oripath, f), "r") as rf:
            lines = rf.readlines()
        lines = [s for s in lines if ".  CORRECT TERMINATION" not in s
                 and len(s.strip()) >= 1]
        if len(outlines) == 0:
            try:
                nbeams = int(lines[1].strip().split()[0])
            except Exception:
                logger.error("Failed to combine fd.out", exc_info=True)
                raise
            outlines += lines
        else:
            outlines += lines[nbeams+2:]
        try:
            os.remove(os.path.join(oripath, f))
        except Exception:
            logger.warning("Failed to delete file " + f)
    try:
        with open(os.path.join(targetpath, "fd.out"), "w") as wf:
            wf.write("".join(outlines))
    except Exception:
        logger.error("Failed to write combine fd.out")
        raise
    return


def readFdOut(readfile="fd.out", for_error=False):
    """Reads the fd.out file produced by the refcalc and returns a list of
    Beam objects."""
    try:
        with open(readfile, 'r') as rf:
            filelines = [line[:-1] for line in rf.readlines()]
            rf.seek(0)
            fdout = rf.read()
    except FileNotFoundError:
        logger.error("Error in readFdOut: file "+str(readfile)+" not found.")
        raise
    if filelines[0].startswith(" SOME ERROR OCCURED WHILE READING"):
        logger.error("File "+readfile+" reports error: "+filelines[0])
        raise Exception("File "+readfile+" reports error.")
    if ".  CORRECT TERMINATION" in fdout:   # happens for superpos output...
        fdout = fdout.split(".  CORRECT TERMINATION")[0]
    theobeams = []
    i = 1   # number lines as in text editor - be careful about indexing!
    nbeams = 1e10   # some large number, just to enter the while loop
    while i < nbeams+3:
        llist = filelines[i-1].split()
        if i == 1:
            pass  # header - skip
        elif i == 2:
            nbeams = int(llist[0])
        else:
            theobeams.append(tl.Beam((float(llist[1]), float(llist[2]))))
        i += 1

    # re-label the beams to get the correct number of characters and formatting
    mw = max([beam.lwidth for beam in theobeams])
    for beam in theobeams:
        beam.lwidth = mw
        beam.getLabel()

    # From now on, all the lines correspond to blocks of beam intensities.
    # Collect them, skipping empty lines
    blocks = []
    for line in filelines[i-1:]:
        llist = line.split()
        if len(llist) == 0:
            continue  # skip empty lines
        try:
            float(llist[0])
        except ValueError:
            break  # end of data
        blocks.extend(llist)

    # Each block contains exactly nbeams+2 entries
    # reshape blocks so that each line corresponds to one block
    blocks = np.reshape(blocks, (-1, nbeams+2))

    # and parse the data
    warned = False
    for block in blocks:
        if block[1] != "0.0001":
            if not for_error and not warned:
                warned = True
                logger.warning("File " + readfile + "contains unexpected "
                               "data for more than one structural "
                               "variation. Only the first block was read.")
            continue
        en = float(block[0])
        values = [float(s) for s in block[2:]]
        for (j, beam) in enumerate(theobeams):
            beam.intens[en] = values[j]
    return theobeams, fdout


def writePARAM(sl, rp, lmax=-1):
    """Creats the contents of the PARAM file for the reference calculation.
    If no LMAX is passed, will use maximum LMAX from rp. Returns str."""
    if lmax == -1:
        lmax = rp.LMAX[1]
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
    m = rp.SUPERLATTICE.copy()
    if m[1, 1] != 0:      # m[1] not parallel to a_bulk
        if m[0, 1] != 0:  # m[0] not parallel to a_bulk
            # find basis in which m[0] is parallel to a_bulk
            f = tl.base.lcm(abs(int(m[0, 1])), abs(int(m[1, 1])))
            m[0] *= f/m[0, 1]
            m[1] *= f/m[1, 1]
            m[0] -= m[1]*np.sign(m[0, 1])*np.sign(m[1, 1])
        nl1 = abs(int(round(m[0, 0])))
        nl2 = abs(int(round(abs(np.linalg.det(rp.SUPERLATTICE))/nl1)))
    else:               # m[1] already parallel to a_bulk
        nl2 = abs(int(round(m[1, 0])))
        nl1 = abs(int(round(abs(np.linalg.det(rp.SUPERLATTICE))/nl2)))
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
    output += '      PARAMETER (MNLTYPE = '+str(len(sl.layers))+')\n'
    mnbrav = 0
    mnsub = 0
    mnstack = 0
    if sl.bulkslab is None:
        sl.bulkslab = sl.makeBulkSlab(rp)
    for layer in [lay for lay in sl.layers if not lay.isBulk]:
        mnstack += 1
        if len(layer.atlist) == 1:
            mnbrav += 1
        if len(layer.atlist) > mnsub:
            mnsub = len(layer.atlist)
    for i, layer in enumerate([lay for lay in sl.layers if lay.isBulk]):
        if len(sl.bulkslab.layers[i].atlist) == 1:
            mnbrav += 1
        if len(sl.bulkslab.layers[i].atlist) > mnsub:
            mnsub = len(layer.atlist)
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


def collectFIN():
    """Combines AUXLATGEO, BEAMLIST, AUXNONSTRUCT, PHASESHIFTS, AUXBEAMS
    and AUXGEO into one string (input for refcalc), which it returns."""
    filenames = ["AUXLATGEO", "BEAMLIST", "AUXNONSTRUCT", "PHASESHIFTS",
                 "AUXBEAMS", "AUXGEO"]
    fin = ""
    for fn in filenames:
        with open(fn, "r") as rf:
            fin += rf.read()
        if fin[-1] != "\n":
            fin += "\n"
    return fin


def writeAUXLATGEO(sl, rp):
    """Writes AUXLATGEO, which is part of the input FIN for the refcalc."""
    output = ''
    output += rp.systemName+' '+rp.timestamp+'\n'
    f72x3 = ff.FortranRecordWriter('3F7.2')
    ens = [rp.THEO_ENERGIES[0], rp.THEO_ENERGIES[1]+0.01, rp.THEO_ENERGIES[2]]
    ol = f72x3.write(ens).ljust(24)
    output += ol + 'EI,EF,DE\n'
    f74x2 = ff.FortranRecordWriter('2F7.4')
    ucsurf = np.transpose(sl.ucell[:2, :2])
    if sl.bulkslab is None:
        sl.bulkslab = sl.makeBulkSlab(rp)
    ucbulk = sl.bulkslab.ucell[:2, :2].T
    ol = f74x2.write(ucbulk[0]).ljust(24)
    output += ol + 'ARA1\n'
    ol = f74x2.write(ucbulk[1]).ljust(24)
    output += ol + 'ARA2\n'
    output += ' 0.0    0.0             SS1\n'
    output += ' 0.0    0.0             SS2\n'
    output += ' 0.0    0.0             SS3\n'
    output += ' 0.0    0.0             SS4\n'
    ol = f74x2.write(ucsurf[0]).ljust(24)
    output += ol + 'ARB1\n'
    ol = f74x2.write(ucsurf[1]).ljust(24)
    output += ol + 'ARB2\n'
    output += ' 0.0    0.0             SO1\n'
    output += ' 0.0    0.0             SO2\n'
    output += ' 0.0    0.0             SO3\n'
    ol = f74x2.write([0.5, rp.V0_Z_ONSET])
    ol = ol.ljust(24)
    output += ol + 'FR ASE\n'
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
    output = ''
    f74 = ff.FortranRecordWriter('F7.4')
    ol = f74.write([rp.ATTENUATION_EPS])
    output += ol+'           >>>>> ! <<<<<              TST\n'
    i4x15 = ff.FortranRecordWriter('15I4')
    ol = i4x15.write(beamnums)
    output += ol+'\n'
    f72f61 = ff.FortranRecordWriter('(F7.2, F6.1)')
    ol = f72f61.write([rp.THETA, rp.PHI]).ljust(45)
    output += ol+'THETA FI\n'
    ol = f74.write([rp.BULKDOUBLING_EPS]).ljust(45)
    output += ol+'EPS\n'
    i3 = ff.FortranRecordWriter('I3')
    ol = i3.write([rp.BULKDOUBLING_MAX]).ljust(45)
    output += ol+'LITER\n'
    ol = i3.write([rp.LMAX[1]]).ljust(45)
    output += ol+'LMAX\n'
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
    output = ''
    output += ('---------------------------------------------------------'
               '----------\n')
    output += ('--- define chem. and vib. properties for different atomic '
               'sites ---\n')
    output += ('---------------------------------------------------------'
               '----------\n')
    i3 = ff.FortranRecordWriter('I3')
    ol = i3.write([len(sl.sitelist)])
    ol = ol.ljust(26)
    output += ol + 'NSITE: number of different site types\n'
    f74x2 = ff.FortranRecordWriter('2F7.4')
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
                        ol = f74x2.write([occ, vib])
                    except Exception:
                        logger.error(
                            "Exception while trying to write occupation / "
                            "vibrational amplitude for site " + site.label,
                            exc_info=True)
                    ol = ol.ljust(26)
                    output += ol + comment + '\n'

    output += ('-----------------------------------------------------'
               '--------------\n')
    output += ('--- define different layer types                     '
               '           ---\n')
    output += ('-----------------------------------------------------'
               '--------------\n')
    ol = i3.write([len(sl.layers)])
    ol = ol.ljust(26)
    output += ol + 'NLTYPE: number of different layer types\n'
    f74x3 = ff.FortranRecordWriter('3F7.4')
    blayers = [lay for lay in sl.layers if lay.isBulk]
    nblayers = [lay for lay in sl.layers if not lay.isBulk]
    layerOffsets = [np.zeros(3) for _ in range(len(sl.layers) + 1)]
    if sl.bulkslab is None:
        sl.bulkslab = sl.makeBulkSlab(rp)
    for i, layer in enumerate(sl.layers):
        output += '-   layer type '+str(i+1)+' ---\n'
        if layer.isBulk:
            output += ('  2                       LAY = 2: layer type no. '
                       + str(i+1) + ' has bulk lateral periodicity\n')
        else:
            output += ('  1                       LAY = 1: layer type no. '
                       + str(i+1) + ' has overlayer lateral periodicity\n')
        if layer.isBulk:
            bl = sl.bulkslab.layers[blayers.index(layer)]
            bulknums = [at.oriN for at in bl.atlist]
            bulkUnique = [at for at in layer.atlist if at.oriN in bulknums]
            natoms = len(bulkUnique)
            # sanity check: ratio of unit cell areas (given simply by
            #  SUPERLATTICE) should match ratio of written vs skipped atoms:
            arearatio = 1 / abs(np.linalg.det(rp.SUPERLATTICE))
            atomratio = len(bulkUnique) / len(layer.atlist)
            if abs(arearatio - atomratio) > 1e-3:
                logger.warning(
                    'Ratio of bulk atoms inside/outside the bulk unit cell '
                    'differs from bulk/slab unit cell size ratio. This means '
                    'that the actual periodicity of the POSCAR does not match '
                    'the periodicity given in the SUPERLATTICE parameter. '
                    'Check SUPERLATTICE parameter and bulk symmetry!')
                rp.setHaltingLevel(2)
        else:
            natoms = len(layer.atlist)
        ol = i3.write([natoms])
        ol = ol.ljust(26)
        output += ol+'number of Bravais sublayers in layer '+str(i+1)+'\n'
        if layer.isBulk:
            writelist = bulkUnique
        else:
            writelist = layer.atlist
        writelist.sort(key=lambda atom: -atom.pos[2])
        for atom in writelist:
            writepos = atom.cartpos - atom.layer.cartori
            if rp.LAYER_STACK_VERTICAL:
                writepos += np.append(atom.layer.cartori[:2]
                                      - blayers[0].cartori[:2], 0)
            ol = i3.write([sl.sitelist.index(atom.site)+1])
            if natoms != 1:
                ol += f74x3.write([writepos[2], writepos[0], writepos[1]])
            else:
                # Bravais layers need to have coordinate (0., 0., 0.)
                #  -> store actual position for later, it will go into the
                #  interlayer vector
                ol += f74x3.write([0., 0., 0.])
                layerOffsets[layer.num] += writepos
                layerOffsets[layer.num + 1] -= writepos
            ol = ol.ljust(26)
            output += ol+'Atom N='+str(atom.oriN)+' ('+atom.el+')\n'
    output += ('--------------------------------------------------------------'
               '-----\n')
    output += ('--- define bulk stacking sequence                             '
               '  ---\n')
    output += ('--------------------------------------------------------------'
               '-----\n')
    output += ('  0                       TSLAB = 0: compute bulk using layer '
               'doubling\n')

    if type(rp.BULK_REPEAT) == np.ndarray:
        bulkc = np.copy(rp.BULK_REPEAT) * np.array([1, 1, -1])
        if bulkc[2] < 0:
            bulkc = -bulkc
        if rp.N_BULK_LAYERS == 2:
            zdiff = bulkc[2] - (blayers[1].cartbotz - blayers[0].cartbotz)
            asaZ = bulkc[2] - (blayers[1].cartbotz - blayers[0].cartori[2])
        else:
            zdiff = bulkc[2]
            asaZ = bulkc[2] - (blayers[0].cartbotz - blayers[0].cartori[2])
        bulkc = bulkc * zdiff/bulkc[2]
        bulkc[2] = asaZ
        bvectors_ASA = bulkc
    else:
        if rp.N_BULK_LAYERS == 2:
            zdiff = rp.BULK_REPEAT - (blayers[1].cartbotz
                                      - blayers[0].cartbotz)
            asaZ = rp.BULK_REPEAT - blayers[1].cartbotz + blayers[0].cartori[2]
        else:
            zdiff = rp.BULK_REPEAT
            asaZ = rp.BULK_REPEAT - blayers[0].cartbotz + blayers[0].cartori[2]
        # zdiff = rp.ASAZ + blayers[0].cartbotz - blayers[0].cartori[2]
        bvectors_ASA = [-sl.ucell[0][2] * zdiff/sl.ucell[2][2],
                        -sl.ucell[1][2] * zdiff/sl.ucell[2][2],
                        asaZ]

    # determine ASBULK - interlayer vector between bulk layers
    if rp.N_BULK_LAYERS == 2:
        # add layerOffsets for Bravais layers:
        bvectors_ASA += layerOffsets[blayers[0].num+1]
        # calculate ASBULK:
        bvectors_ASBULK = blayers[1].cartori - blayers[0].cartori
        bvectors_ASBULK[2] = blayers[1].cartori[2] - blayers[0].cartbotz
        bl2num = blayers[1].num
        # add layerOffsets for Bravais layers:
        bvectors_ASBULK -= layerOffsets[blayers[0].num+1]
        # bvectors_ASBULK += (layerOffsets[blayers[1].num+1]
        #                     + layerOffsets[blayers[0].num])
    else:
        bl2num = blayers[0].num
        bvectors_ASBULK = bvectors_ASA

    ol = f74x3.write([bvectors_ASA[2], bvectors_ASA[0], bvectors_ASA[1]])
    ol = ol.ljust(26)
    output += ol + 'ASA interlayer vector between different bulk units\n'
    ol = i3.write([blayers[0].num+1])
    ol = ol.ljust(26)
    output += (ol + 'top layer of bulk unit: layer type '+str(blayers[0].num+1)
               + '\n')
    ol = i3.write([bl2num+1])
    ol = ol.ljust(26)
    output += ol + 'bottom layer of bulk unit: layer type '+str(bl2num+1)+'\n'
    ol = f74x3.write([bvectors_ASBULK[2], bvectors_ASBULK[0],
                      bvectors_ASBULK[1]])
    ol = ol.ljust(26)
    output += (ol + 'ASBULK between the two bulk unit layers (may differ from '
               'ASA)\n')
    output += ('--------------------------------------------------------------'
               '-----\n')
    output += ('--- define layer stacking sequence and Tensor LEED output     '
               '  ---\n')
    output += ('--------------------------------------------------------------'
               '-----\n')
    nonbulk = len(sl.layers)-rp.N_BULK_LAYERS
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
    ol = i3.write([len(sl.layers)-rp.N_BULK_LAYERS])
    ol = ol.ljust(26)
    output += ol + 'NSTACK: number of layers stacked onto bulk\n'
    for layer in list(reversed(nblayers)):
        n = layer.num + 1
        if not rp.LAYER_STACK_VERTICAL:
            v = sl.layers[n].cartori - layer.cartori
        else:
            v = np.zeros(3)
        v[2] = sl.layers[n].cartori[2] - layer.cartbotz
        v = v + layerOffsets[n]   # add layerOffsets for Bravais layers
        ol = i3.write([n]) + f74x3.write([v[2], v[0], v[1]])
        ol = ol.ljust(26)
        output += (ol + 'layer '+str(n)+': layer type '+str(n)+', interlayer '
                   'vector below\n')     # every layer is also a layer type
        ol = i3.write([rp.TENSOR_OUTPUT[layer.num]])
        ol = ol.ljust(26)
        output += (ol + '0/1: Tensor output is required for this layer '
                   '(TENSOR_OUTPUT)\n')
        i = 1
        for atom in layer.atlist:
            # ol = 'T_'+atom.el+str(atom.oriN)
            ol = 'T_'+str(atom.oriN)
            ol = ol.ljust(26)
            output += (ol + 'Tensor file name, current layer, sublayer '+str(i)
                       + '\n')
            i += 1
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


def writeMuftin(sl, rp):
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
    oline = "      VV = "+rp.V0_REAL
    output += tl.base.fortranContLine(oline) + "\n"
    output += """
      write(6,*) workfn, EEV
      write(6,*) VV

c  set difference between bulk and overlayer inner potential

      VO = 0.

c  set imaginary part of inner potential - energy independent value used here

"""
    oline = "      VPI = "+str(rp.V0_IMAG)
    output += tl.base.fortranContLine(oline) + "\n"
    output += """
C  set substrate / overlayer imaginary part of inner potential

"""
    oline = "      VPIS = "+str(rp.V0_IMAG)
    output += tl.base.fortranContLine(oline) + "\n"
    oline = "      VPIO = "+str(rp.V0_IMAG)
    output += tl.base.fortranContLine(oline) + "\n"
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
