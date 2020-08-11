# -*- coding: utf-8 -*-

"""
Created on Aug 11 2020

@author: Florian Kraushofer

Tensor LEED Manager section Initialization
"""

import os
import logging
import copy

import tleedmlib as tl

logger = logging.getLogger("tleedm.initialization")

def initialization(sl, rp):
    """Runs the initialization. Returns 0 on success."""
    # check for experimental beams:
    expbeamsname = ""
    for fn in ["EXPBEAMS.csv", "EXPBEAMS"]:
        if os.path.isfile(fn):
            expbeamsname = fn
            break
    if expbeamsname:
        if len(rp.THEO_ENERGIES) == 0:
            er = []
        else:
            er = rp.THEO_ENERGIES[:2]
        if not rp.fileLoaded["EXPBEAMS"]:
            try:
                rp.expbeams = tl.readOUTBEAMS(fn, enrange=er)
                rp.fileLoaded["EXPBEAMS"] = True
            except:
                logger.error("Error while reading file "+fn, exc_info=True)
    rp.initTheoEnergies()  # may be initialized based on exp. beams
    # check whether _PHASESHIFTS are present & consistent:
    newpsGen, newpsWrite = True, True
                          # new phaseshifts need to be generated/written
    if os.path.isfile(os.path.join(".","_PHASESHIFTS")):
        try:
            (rp.phaseshifts_firstline, rp.phaseshifts,
                 newpsGen, newpsWrite) = tl.readPHASESHIFTS(sl, rp)
        except:
            logger.warning("Found a _PHASESHIFTS file but could not "
                "read it. A new _PHASESHIFTS file will be generated."
                "The exception during read was: ", exc_info=True)
            rp.setHaltingLevel(1)
    if newpsGen:
        try:
            rundgrenpath = os.path.join('.', 'source', 'EEASiSSS.x')
            serneliuspath = os.path.join('.', 'source', 'seSernelius')
            logger.info("Generating phaseshifts data... ")
            (rp.phaseshifts_firstline, 
                        rp.phaseshifts) = tl.runPhaseshiftGen(sl, rp,
                                               psgensource = rundgrenpath, 
                                               excosource = serneliuspath)
            logger.debug("Finished generating phaseshift data")
        except:
            logger.error("Exception while calling phaseshiftgen: ")
            raise
    if newpsWrite:
        try:
            tl.writePHASESHIFTS(rp.phaseshifts_firstline, rp.phaseshifts)
        except:
            logger.error("Exception during writePHASESHIFTS: ")
            raise
    rp.fileLoaded["PHASESHIFTS"] = True
    rp.updateDerivedParams()
    rp.manifest.append("_PHASESHIFTS")

    # if necessary, run findSymmetry:
    if sl.planegroup == "unknown":
        sl.findSymmetry(rp)
        sl.enforceSymmetry(rp)
    
    # generate new POSCAR
    tmpslab = copy.deepcopy(sl)
    tmpslab.sortOriginal()
    try:
        tl.writeCONTCAR(tmpslab, filename='POSCAR', comments='all')
    except:
        logger.error("Exception occurred while writing new POSCAR")
        raise
    rp.manifest.append('POSCAR')
    # generate POSCAR_oricell
    tmpslab.revertUnitCell()
    try:
        tl.writeCONTCAR(tmpslab, filename='POSCAR_oricell', comments='nodir')
    except:
        logger.error("Exception occurred while writing POSCAR_oricell, "
                      "execution will continue...")

    # create bulk slab:
    if sl.bulkslab == tl.DEFAULT:
        sl.bulkslab = sl.makeBulkSlab(rp)
    bsl = sl.bulkslab
    # find minimum in-plane unit cell for bulk:
    logger.info("Checking bulk unit cell...")
    changecell, mincell = bsl.getMinUnitCell(rp)
    if changecell:
        sl.changeBulkCell(rp, mincell)
        bsl = sl.bulkslab
    if not rp.superlattice_defined:
        ws = tl.writeWoodsNotation(rp.SUPERLATTICE) 
                # !!! replace the writeWoodsNotation from baselib with 
                #   the one from guilib
        si = rp.SUPERLATTICE.astype(int)
        if ws:
            logger.info("Found SUPERLATTICE = "+ws)
        else:
            logger.info("Found SUPERLATTICE M = {} {}, {} {}".format(
                                    si[0,0], si[0,1], si[1,0], si[1,1]))
    
    # bulk plane group detection:
    logger.info("Initializing bulk symmetry search...")
    bsl.findSymmetry(rp, bulk=True, output=False)
    bsl.revertUnitCell() # keep origin matched with main slab
    logger.info("Found bulk plane group: "+bsl.foundplanegroup)
    bsl.findBulkSymmetry(rp)
    
    # write POSCAR_bulk
    bsl = copy.deepcopy(sl.bulkslab)
    bsl.sortOriginal()
    try:
        tl.writeCONTCAR(bsl, filename='POSCAR_bulk', comments='bulk')
    except:
        logger.error("Exception occurred while writing POSCAR_bulk")
        raise
    
    # generate beamlist
    logger.info("Generating _BEAMLIST...")
    try:
        bgenpath = os.path.join('.', 'source', 'beamgen3.out')
        tl.runBeamGen(sl,rp,beamgensource = bgenpath)
        # this does NOT read the resulting file!
    except:
        logger.error("Exception occurred while calling beamgen.")
        raise
    rp.manifest.append("_BEAMLIST")
    try:
        rp.beamlist = tl.readBEAMLIST()
        rp.fileLoaded["BEAMLIST"] = True
    except:
        logger.error("Error while reading required file _BEAMLIST")
        raise
    
    tl.writePatternInfo(sl, rp)
    
    # if EXPBEAMS was loaded, it hasn't been check yet - check now
    if rp.fileLoaded["EXPBEAMS"]:
        tl.checkEXPBEAMS(sl, rp)
    # write and sort IVBEAMS
    if not rp.fileLoaded["IVBEAMS"]:
        try:
            rp.ivbeams = tl.writeIVBEAMS(sl, rp)
            rp.ivbeams_sorted = False
            rp.fileLoaded["IVBEAMS"] = True
            rp.manifest.append("IVBEAMS")
        except:
            logger.error("Error while writing IVBEAMS file based on "
                          "EXPBEAMS data.")
            raise
    if not rp.ivbeams_sorted:
        rp.ivbeams = tl.sortIVBEAMS(sl, rp)
        rp.ivbeams_sorted = True
    return 0