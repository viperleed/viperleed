# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 15:32:12 2020
major rework Sep&Oct 2021

@author: Florian Kraushofer & Alexander M. Imre
"""

import os
import numpy as np
import logging
import random
import subprocess
import re
import shutil

import viperleed.tleedmlib as tl
from viperleed import fortranformat as ff

logger = logging.getLogger("tleedm.psgen")


def runPhaseshiftGen_old(sl, rp,
                     psgensource=os.path.join('tensorleed', 'EEASiSSS.x'),
                     excosource=os.path.join('tensorleed', 'seSernelius'),
                     atdenssource=os.path.join('tensorleed',
                                               'atom_density_files')):
    """Creates required input for EEASiSSS.x, then runs it. Reads the output
    files and extracts information for PHASESHIFTS file, then returns that
    information (without writing PHASESHIFTS)."""

    test_new = True
    '''
    if test_new:
        runPhaseshiftGen_new(sl, rp,
                             psgensource=os.path.join('tensorleed', 'EEASiSSS.x'),
                             excosource=os.path.join('tensorleed', 'seSernelius'),
                             atdenssource=os.path.join('tensorleed',
                                                       'atom_density_files'))
    '''
    shortpath = rp.sourcedir
    if len(os.path.relpath(rp.sourcedir)) < len(shortpath):
        shortpath = os.path.relpath(rp.sourcedir)

    if len(shortpath) > 62:
        # too long - need to copy stuff here
        manual_copy = True
        os.makedirs("tensorleed", exist_ok=True)
        shutil.copy2(os.path.join(shortpath, excosource), excosource)
        shortpath = "."
    else:
        manual_copy = False

    psgensource = os.path.join(rp.sourcedir, psgensource)
    excosource = os.path.join(shortpath, excosource)

    lmax = 16   # this could be a variable, for now set fixed...
    nsl, newbulkats = sl.addBulkLayers(rp)
    outvals = {}
    # dict containing lists of values to output: outvals[energy][block][L]

    # The following originally calculated phaseshifts first for only the bulk,
    #   then for the slab, taking phaseshifts from the bulk calculations if the
    #   site is present at all in the bulk. However, this means that different
    #   muffin tin parameters are applied for the different phaseshifts.
    # New solution:
    #   Take phaseshifts from the slab calculation only. Average only over
    #   atoms NOT added as "new bulk" ("newbulkats" list above) to avoid
    #   influence of the (false) bottom "surface".
#    for bulk in [True,False]:
#    wsl = bsl if bulk else nsl  # working slab for this iteration
    # before starting on unit cell, determine whether supercell is needed:
    blocks = []  # tuples of sites (with poscar elements) and elements
    #        (same as POSCAR if no ELEMENT_MIX, ELEMENT_MIX elements if not)
    for site in nsl.sitelist:
        if site.el in rp.ELEMENT_MIX:
            for el in rp.ELEMENT_MIX[site.el]:
                blocks.append((site, el))
        else:
            blocks.append((site, site.el))
    scsize = 1 # Supercell size
    if len(rp.ELEMENT_MIX) > 0:
        minnum = -1
        for (site, el) in [(site, el) for (site, el) in blocks if site.el
                           in rp.ELEMENT_MIX and (site.occ[el] > 0. or
                                                  el in site.mixedEls)]:
            al = [at for at in nsl.atlist if at.site == site]
            atcount = len(al)*site.occ[el]
            if minnum < 0 or (minnum > atcount >= 0):
                minnum = atcount
        # we want at least 2 atoms of each element in each site type:
        if 0 < minnum < 2.0:
            scsize = int(np.ceil(2.0/minnum))
        elif minnum == 0:
            scsize = 100  # large number, will be decreased below
    if scsize > 1:  # some checks to make sure it doesn't get too large
        maxcells = 20  # upper limit on supercell size
        maxats = 500   # upper limit on atoms in supercell
        if scsize > maxcells:
            scsize = maxcells
            # don't warn - this is a large unit cell either way.
        if len(nsl.atlist) * scsize > maxats:
            logger.debug(
                "Phaseshift generation: Given element "
                "concentrations would require a very large supercell. "
                "Element concentrations for low-occupancy elements will be "
                "increased to avoid this. This only concerns the phaseshifts "
                "calculation and should not cause problems.")
            # determine minimum size to have 2 of each element
            minsize = 1
            for site in [s for s in nsl.sitelist if s.el in rp.ELEMENT_MIX]:
                ats = len([at for at in nsl.atlist if at.site == site])
                els = len([el for el in rp.ELEMENT_MIX[site.el]
                           if site.occ[el] > 0.])
                minsize = max(minsize, int(np.ceil(2*els / ats)))
            scsize = max(minsize, int(maxats / len(nsl.atlist)))

    subatlists = {}     # atlist per block tuple
    if scsize > 1:  # construct supercell to get enough atoms
        xsize = int(np.ceil(np.sqrt(scsize)))  # if scsize is not prime, try
        while scsize % xsize != 0:             # making it close to square
            xsize += 1
        ysize = int(scsize / xsize)
        cpatlist = nsl.atlist[:]
        for at in cpatlist:
            for i in range(0, xsize):
                for j in range(0, ysize):
                    if i == j == 0:
                        continue
                    tmpat = at.duplicate()
                    tmpat.pos[0] += i
                    tmpat.pos[1] += j
        nsl.getCartesianCoordinates()
        nsl.ucell = np.dot(np.array([[xsize, 0, 0], [0, ysize, 0], [0, 0, 1]]),
                           nsl.ucell)
        nsl.getFractionalCoordinates()

    for site in nsl.sitelist:
        if site.el in rp.ELEMENT_MIX:
            occdict = {}
            for (k, v) in site.occ.items():
                if v > 0.0 or k in rp.ELEMENT_MIX[site.el]:
                    occdict[k] = v
            # sort by occupancy values
            occdict = dict(sorted(occdict.items(),
                                  key=lambda kv: (kv[1], kv[0])))
            al = [at for at in nsl.atlist if at.site == site]
            totats = len(al)
            for el in occdict:
                subatlists[(site, el)] = []
                reqats = int(np.ceil(totats * site.occ[el]))
                reqats = max(2, reqats)
                while reqats > 0 and len(al) > 0:
                    at = random.choice(al)
                    subatlists[(site, el)].append(at)
                    al.remove(at)
                    reqats -= 1
            if len(al) > 0:  # should never happen
                logger.warning("Error in PHASESHIFTS file "
                               "generation: Not all atoms were distributed!")
        else:
            subatlists[(site, site.el)] = [at for at in nsl.atlist
                                           if at.site == site]
    blocks = [(site, el) for (site, el) in blocks
              if len(subatlists[(site, el)]) > 0]
#    if bulk:
#        bulksites = [site for (site,el) in blocks]
#    else:
#        # site objects are different in the new slab; copy to equivalent
#        oldbulksites = bulksites[:]
#        bulksites = []
#        for oldsite in oldbulksites:
#            for newsite in [site for (site,el) in blocks]:
#                if newsite.isEquivalent(oldsite):
#                    bulksites.append(newsite)
    # start writing output, will be input for EEASiSSS code:
    output = ''
    output += "STRUCTURE:\n"
    output += rp.systemName+" "+rp.timestamp+"\n"
#    if bulk:
#        output += ("'b'  1.00 16      BulkOrSlab('b' or 's'), "
#                  "LatticeConstant(Angstroms), nshell\n")
#    else:
    output += ("'s'  1.00 16      BulkOrSlab('b' or 's'), "
               "LatticeConstant(Angstroms), nshell\n")
    uct = nsl.ucell.transpose()
    for i in range(0, 3):
        ol = ''
        for j in range(0, 3):
            s = str(round(uct[i, j], 4))+' '
            ol += s.ljust(8)
        if i == 0:
            ol += '      CoordinatesOfUnitCell(LCunits)\n'
        else:
            ol += '\n'
        output += ol

    output += (str(len(nsl.atlist))+"  "+str(len(nsl.atlist))
               + "                  #AtomTypes,#OccupiedAtomTypes\n")
    ptl = [el.lower() for el in tl.leedbase.periodic_table]

    chemels = {}
    chemelspaths = {}
    for (site, el) in blocks:
        if el in rp.ELEMENT_RENAME:
            chemel = rp.ELEMENT_RENAME[el]
        elif el.lower() in ptl:
            chemel = el.capitalize()
        else:
            logger.error("Error generating PHASESHIFTS file: Could not "
                         "identify "+el+" as a chemical element. Define "
                         "ELEMENT_RENAME or ELEMENT_MIX parameter.")
            raise
        subpath = os.path.join(atdenssource, chemel, "chgden"+chemel)
        chgdenrelpath = os.path.join(shortpath, subpath)
        if manual_copy:
            os.makedirs(os.path.join(os.path.dirname(subpath)), exist_ok=True)
            shutil.copy2(os.path.join(rp.sourcedir, subpath), subpath)
        # if os.name == 'nt':     # windows - replace the backslashes.
        #     chgdenrelpath = chgdenrelpath.replace('/', '\\')
        chemels[el] = chemel
        chemelspaths[el] = chgdenrelpath

    nsl.sortByZ(botToTop=True)
    for at in nsl.atlist:
        # realcartpos = np.dot(nsl.ucell, at.pos)
        # use the "real" cartesian system, with Z going up
        for (site, el) in blocks:
            if at in subatlists[(site, el)]:
                chemel = chemels[el]
                chgdenpath = chemelspaths[el]
        output += ("1 "+str(tl.leedbase.periodic_table.index(chemel)+1)
                   + ".  0.  0.  '"+chgdenpath+"'\n")
        ol = ""
        for j in range(0, 3):
            # ol += str(round(realcartpos[j],4))+" "
            ol += str(round(at.cartpos[j], 4))+" "
        output += ol + "     Coordinates(LCunits)\n"
    output += "SCATTERING: \n"
    output += "'"+excosource+"' |exchange-correlation file\n"
    output += str(lmax)+"  |lmax\n"
    output += "'r' |SelectCalculation: 'relativistic'/'nonrelativistic'\n"
    output += "'p' |SelectOutput: 'phaseshift'/'sigma'/'dataflow'\n"
    output += "'n' |phaseshift: print_SpinPhaseShift? 'yes'/no'\n"
    output += ("'n' |phaseshift,sigma: print_log10(DsigmaDomega)? "
               "'yes'/'no'\n")
    output += "'n' |phaseshift,sigma: print_DsigmaDtheta? 'yes'/'no'\n"
    output += "'n' |phaseshift,sigma: print_Sherman? 'yes'/'no'\n"
    output += ("'n' |phaseshift,sigma: print_TotalCrossSection? "
               "'yes'/'no'\n")
    output += "'n' |dataflow: print_RhoPot? 'yes'/'no'\n"
    output += "'n' |dataflow: print_PSvsE? 'yes'/'no'\n"
    output += "'n' |dataflow: print_WaveFunction? 'yes'/'no'\n"
    output += ("0.0," + str(round(float(rp.THEO_ENERGIES[1]) + 20, 1)) # hardcoded lower boundary of 0
               + ","+str(round(float(rp.THEO_ENERGIES[2]), 1))
               + " |'phaseshift'/'dataflow' run: E1->E2,PS_Estep\n")
    output += "100.,100.,1  |'sigma' run: E1->E2,NumVals\n"
    output += "MT OPTIMIZATION (Nelder-Mead):\n"
    output += "'*'           | NMselect: '*' or '+'\n"
    output += "0.05d0   | NMlambda, simplex point shift.\n"
    output += "3            | NMiter, no. of simplex calls.\n"
    output += "1.d-04    ! NMeps,  simplex epsilon.\n"
    output += "0.125d0 | NMepsit, epsilon iteration, eps=eps*epsit.\n"
#    outfilename = 'eeasisss-input-bulk' if bulk else 'eeasisss-input-slab'
    outfilename = 'eeasisss-input'
    try:
        with open(outfilename, 'w') as wf:
            wf.write(output)
    except Exception:
        logger.error("Phaseshift data generation: Failed to write "
                     + outfilename + ". Proceeding with execution...")
    # if os.name == 'nt':
    #     logger.error("Phaseshift generation is currently not "
    #                  "supported on Windows. Use a linux shell to run "
    #                  "phaseshift generation script.")
    #     raise EnvironmentError("Phaseshift generation is currently not "
    #                            "supported on Windows.")
    # else:

    ##################################
    subprocess.run(psgensource, input=output, encoding='ascii') # RUNS phaseshift programm

    # go through all the files that were generated by EEASiSSS and read
    filelist = [filename for filename in os.listdir('.') if
                filename.startswith('PS.r.')]
    rgx = re.compile(r'PS\.r\.[0-9]+\.[0-9]+')
    remlist = []
    for filename in filelist:
        m = rgx.match(filename)
        if not m:
            remlist.append(filename)
        else:
            if m.group(0) != filename:
                remlist.append(filename)
            else:
                try:
                    int(filename.split('.')[-1])
                except Exception:
                    remlist.append(filename)
    for filename in remlist:
        filelist.remove(filename)
    if not filelist:
        logger.error("Phaseshift generation failed: No output files found.")
        raise RuntimeError("Phaseshift generation failed.")

    # sort by atom number
    filelist.sort(key=lambda filename: int(filename.split('.')[-1]))
    firstline = "" # first line contains parameters for the potential
    # data from all the PS files
    atoms_phaseshifts = [[] for i in range(0, len(filelist))]
    for (i, filename) in enumerate(filelist):
        # if bulk or wsl.atlist[i].site not in bulksites:
        if nsl.atlist[i] not in newbulkats: # only atoms that were not added in the new bulk
            psfile = open(filename, 'r')
            reade = -1.
            pslist = []
            ps_of_e = {}
            # dictionary of pslist for energy, where pslist is a list of
            #  floats as read from the PS-file
            for (j, line) in enumerate(psfile):
                if j == 0:
                    # if not bulk:
                    if firstline == "":         # only saved once, since same in all files
                        firstline = line[2:]  # ignore I2 at beginning
                else:   # line should contain data
                    values = [float(s) for s in line.split()]   # put phaseshift floats into values
                    if reade < 0 or values[0] == reade + rp.THEO_ENERGIES[2]: # this takes care of multi-line pshifts
                        # new energy value, start phaseshift value list
                        if reade >= 0:
                            ps_of_e[reade] = pslist
                        reade = values[0]
                        pslist = values[1:]
                    else:
                        pslist.extend(values)
            ps_of_e[reade] = pslist
            psfile.close()
            atoms_phaseshifts[i] = ps_of_e
        os.remove(os.path.join('.', filename))
    for (site, el) in blocks:   # does the averaging over atoms of same element
        writeblock = False
        pssum = None
        for (i, at) in enumerate(nsl.atlist):
            if at in newbulkats or at not in subatlists[(site, el)]: # ignored because added before
                continue
            # if bulk or at.site not in bulksites:
            if pssum is None:
                writeblock = True
                pel = at.el  # POSCAR element for block
                pssum = atoms_phaseshifts[i]
                n = 1
            else:
                for en in atoms_phaseshifts[i]:
                    for j in range(0, len(atoms_phaseshifts[i][en])):
                        pssum[en][j] += atoms_phaseshifts[i][en][j]
                n += 1
        if writeblock:
            # append the values for the whole block to outvals:
            for en in pssum:
                pssum[en] = [v/n for v in pssum[en]]
                if en not in outvals:
                    outvals[en] = []
                outvals[en].append((pel, el, site, pssum[en]))
    # clean up
    # bss = "-bulk" if bulk else "-slab"
    bss = ""
    try:
        os.rename('logfile', 'eeasisss-logfile'+bss)
    except Exception:
        logger.warning("Failed to rename phaseshift generation file "
                       "'logfile' to 'eeasisss-logfile"+bss+"'")
    try:
        os.rename('QMTvsE', 'eeasisss-QMTvsE'+bss)
    except Exception:
        logger.warning("Failed to rename phaseshift generation file "
                       "'QMTvsE' to 'eeasisss-QMTvsE"+bss+"'")
    try:
        os.rename('RMTvsE', 'eeasisss-RMTvsE'+bss)
    except Exception:
        logger.warning("Failed to rename phaseshift generation file "
                       "'RMTvsE' to 'eeasisss-RMTvsE"+bss+"'")
    try:
        os.rename('V0vsE', 'eeasisss-V0vsE'+bss)
    except Exception:
        logger.warning("Failed to rename phaseshift generation file "
                       "'V0vsE' to 'eeasisss-V0vsE"+bss+"'")
    # sort blocks in outvals:
    outvalsSorted = {}
    outvalLength = 0
    # Due to a bug in eeasisss, it does not generate exactly the energies it
    # is asked to; this can result in different energies being present for
    # the bulk and the slab calculations. Workaround: Discard energies where
    # not all sites are present.
    for en in outvals:
        outvalsSorted[en] = []
        if len(outvals[en]) > outvalLength:
            outvalLength = len(outvals[en])
    # first sort by POSCAR elements, same order as POSCAR:
    for el in sl.elements:
        if el in rp.ELEMENT_MIX:
            chemelList = rp.ELEMENT_MIX[el]
        else:
            chemelList = [el]
        # then by sites:
        siteList = [site for site in sl.sitelist if site.el == el]
        for cel in chemelList:
            for site in siteList:
                for en in outvalsSorted:
                    for (o_el, o_cel, o_site, pslist) in outvals[en]:
                        if (o_el == el and o_cel == cel
                                and o_site.isEquivalent(site)):
                            outvalsSorted[en].append(pslist)
    # return:
    # writePHASESHIFTS wants values as list of tuples and energies in Hartree:
    phaseshifts = []
    for en in outvalsSorted:
        if len(outvalsSorted[en]) == outvalLength:
            # drop energies where phaseshift was not calculated for all sites
            phaseshifts.append((en/27.211396, outvalsSorted[en])) # conversion eV to Hartree
    if firstline == "":
        logger.error("Could not find first line for PHASESHIFTS file "
                     "(should contain MUFTIN parameters).")
        firstline = "ERROR: first line not found in EEASiSSS.x output\n"
        rp.setHaltingLevel(2)
    else:
        # add number of blocks to firstline
        nblocks = len(phaseshifts[0][1])
        firstline = str(nblocks).rjust(3)+"  "+firstline
        # remove the "PS.r.**.**"
        firstline = re.sub(r"PS\.r\.[0-9]+\.[0-9]+", "", firstline)
    return (firstline, phaseshifts)

def runPhaseshiftGen(sl, rp, psgensource=os.path.join('tensorleed', 'eeasisss_new', 'eeas')):
    """
    Creates the input for EEASISSS, runs it and read the output from the generated files.
    Returns phaseshift information (but does not yet write PHASSHIFTS file) in addition with the firstline string which
    contains the muffin tin potential paramenters.
    """
    nsl, newbulkats = sl.addBulkLayers(rp)
    outvals = {}
    # dict containing lists of values to output: outvals[energy][block][L]

    blocks = []  # tuples of sites (with poscar elements) and elements
    #        (same as POSCAR if no ELEMENT_MIX, ELEMENT_MIX elements if not)
    for site in nsl.sitelist:
        if site.el in rp.ELEMENT_MIX:
            for el in rp.ELEMENT_MIX[site.el]:
                blocks.append((site, el))
        else:
            blocks.append((site, site.el))
    scsize = 1 # Super cell size

    if len(rp.ELEMENT_MIX) > 0:
        minnum = -1
        for (site, el) in [(site, el) for (site, el) in blocks if site.el
                                                                  in rp.ELEMENT_MIX and (site.occ[el] > 0. or
                                                                                         el in site.mixedEls)]:
            al = [at for at in nsl.atlist if at.site == site]
            atcount = len(al) * site.occ[el]
            if minnum < 0 or (minnum > atcount >= 0):
                minnum = atcount
        # we want at least 2 atoms of each element in each site type:
        if 0 < minnum < 2.0:
            scsize = int(np.ceil(2.0 / minnum))
        elif minnum == 0:
            scsize = 100  # large number, will be decreased below
    if scsize > 1:  # some checks to make sure it doesn't get too large
        maxcells = 20  # upper limit on supercell size
        maxats = 500  # upper limit on atoms in supercell
        if scsize > maxcells:
            scsize = maxcells
            # don't warn - this is a large unit cell either way.
        if len(nsl.atlist) * scsize > maxats:
            logger.debug(
                "Phaseshift generation: Given element "
                "concentrations would require a very large supercell. "
                "Element concentrations for low-occupancy elements will be "
                "increased to avoid this. This only concerns the phaseshifts "
                "calculation and should not cause problems.")
            # determine minimum size to have 2 of each element
            minsize = 1
            for site in [s for s in nsl.sitelist if s.el in rp.ELEMENT_MIX]:
                ats = len([at for at in nsl.atlist if at.site == site])
                els = len([el for el in rp.ELEMENT_MIX[site.el]
                           if site.occ[el] > 0.])
                minsize = max(minsize, int(np.ceil(2 * els / ats)))
            scsize = max(minsize, int(maxats / len(nsl.atlist)))

    subatlists = {}  # atlist per block tuple
    if scsize > 1:  # construct supercell to get enough atoms
        xsize = int(np.ceil(np.sqrt(scsize)))  # if scsize is not prime, try
        while scsize % xsize != 0:  # making it close to square
            xsize += 1
        ysize = int(scsize / xsize)
        cpatlist = nsl.atlist[:]
        for at in cpatlist:
            for i in range(0, xsize):
                for j in range(0, ysize):
                    if i == j == 0:
                        continue
                    tmpat = at.duplicate()
                    tmpat.pos[0] += i
                    tmpat.pos[1] += j
        nsl.getCartesianCoordinates()
        nsl.ucell = np.dot(np.array([[xsize, 0, 0], [0, ysize, 0], [0, 0, 1]]),
                           nsl.ucell)
        nsl.getFractionalCoordinates()

    for site in nsl.sitelist:
        if site.el in rp.ELEMENT_MIX:
            occdict = {}
            for (k, v) in site.occ.items():
                if v > 0.0 or k in rp.ELEMENT_MIX[site.el]:
                    occdict[k] = v
            # sort by occupancy values
            occdict = dict(sorted(occdict.items(),
                                  key=lambda kv: (kv[1], kv[0])))
            al = [at for at in nsl.atlist if at.site == site]
            totats = len(al)
            for el in occdict:
                subatlists[(site, el)] = []
                reqats = int(np.ceil(totats * site.occ[el]))
                reqats = max(2, reqats)
                while reqats > 0 and len(al) > 0:
                    at = random.choice(al)
                    subatlists[(site, el)].append(at)
                    al.remove(at)
                    reqats -= 1
            if len(al) > 0:  # should never happen
                logger.warning("Error in PHASESHIFTS file "
                               "generation: Not all atoms were distributed!")
        else:
            subatlists[(site, site.el)] = [at for at in nsl.atlist
                                           if at.site == site]
    blocks = [(site, el) for (site, el) in blocks
              if len(subatlists[(site, el)]) > 0]

    ###########
    # subdirecory with phaseshifts
    ps_outdir = 'PS_out'

    # Function below takes care of formatting the input so it can be used with eeasisss
    # calls function to format and organize atoms
    atom_types, input_file_txt, lmax, E2, Estep = format_eeasisss_input(nsl, newbulkats, rp)

    # Write input file
    input_file_name = "EEASISSS-input.txt" # hardcoded, change?
    try:
        with open(input_file_name, 'w') as wf:
            wf.write(input_file_txt)
    except Exception:
        logger.error("Phaseshift data generation: Failed to write "
                     + input_file_name + ". Proceeding with execution...")

    #################
    # Call EEASISSS with input file
    #################
    log_filename = "EEASISSS-log.txt" # log file
    psgensource = os.path.join('tensorleed', 'eeasisss_new', 'eeasisss')
    psgensource = os.path.join(rp.sourcedir, psgensource) # otherwise the location would not be known
    atlib_dir = os.path.join('tensorleed', 'eeasisss_new', 'atlib/') # atom density files, by Sernelius
    atlib_dir = os.path.join(rp.sourcedir, atlib_dir)
    outdir_path = os.path.join(".",ps_outdir+"/")
    # create directory for individual phaseshift files if not yet present
    os.makedirs(outdir_path, exist_ok=True)

    # execution of EEASISSS requires the executable eeas to be in the work directory
    # (a shell script is created and executes eeas from the work directory, which is hardcoded)
    # -> copy eeas to work directory beforehand and remove it afterwards
    eeasisss_exec_path = os.path.join('tensorleed', 'eeasisss_new', 'eeasisss')
    eeasisss_exec_path = os.path.join(rp.sourcedir, eeasisss_exec_path)
    try:
        shutil.copy2(eeasisss_exec_path, ".")
    except Exception:
        logger.error("Could not copy file eeas required by EEASISSS. Phaseshift generation will fail if not present.")
        rp.setHaltinglevel(2)

    # We are now ready to call EEASISS
    psgencommand = [psgensource, '-i', input_file_name, '-l', log_filename, '-a', atlib_dir, '-o', outdir_path]
    logger.debug("Now calling EEASISSS...")
    try:
        complete = subprocess.run(psgencommand, capture_output=True, text=True)  # returns Subprocess.complete instance
        complete.check_returncode()  # If returncode is non-zero, raise a CalledProcessError. -> except & finally
        logger.debug("EEASISSS execution finished without error. See EEASISSS-log.txt") # -> only if returncode == 0
    except Exception:  # can subprocess even fail?
        logger.error("Error during EEASISSS execution.")
        raise RuntimeError("Subprocess EEASISSS failed")
    finally:
        if complete.stdout != "":
            logger.debug("EEASISSS stdout:\n" + str(complete.stdout))
        if complete.stderr != "":
            logger.error("EEASISSS stderr:\n" + str(complete.stderr)) # Put error from Fortran into log!

    #when done, remove eeas executable from work directory again:
    try:
        os.remove('eeasisss')
    except Exception:
        logger.warning("Could not remove eeasisss executable from work directory.") # Not a big deal if this happens.

    # Now the results of the phaseshift calculation are read out and formatted for further processing
    firstline, phaseshifts = convert_eeasisss_output(sl, rp, atom_types, lmax, E2, Estep, ps_outdir)

    # Finish up by moving relevant files to work directory and then removing the out_dir
    move_EEASISSS_files(ps_outdir, log_filename, Vxc0_files = ['Vxc0Einc', 'Vxc0EincAprx', 'Vxc0EincAprx_v0coef'])

    return (firstline, phaseshifts)


def format_eeasisss_input(nsl, newbulkats, rp):
    """
    Produces the input file for EEASISSS (version from 2021).
    """
    ##############################
    # Create EEASiSSS input
    ##############################
    logger.debug("Creating EEASISSS input")
    input_file_header = ""  # Header with info about project
    input_file_header += "STRUCTURE:\n"
    input_file_header += "! " + rp.systemName + " " + rp.timestamp + "\n"
    input_file_header += str(os.getcwd()) + "\n"
    #######
    # Unit cell block
    #######
    input_file_uc_block = ""  # unit cell block to be written into input file
    # Unit of length
    input_file_length_line = '1.889727\t!UnitOfLength conversion to Bohr radii\n'  # EEAS converts to Bohr radii internally
    # and thus needs to be given a conversion factor
    # below adapted from old code
    uct = nsl.ucell.transpose()
    for i in range(0, 3):
        ol = ''
        for j in range(0, 3):
            s = '{: .6f}'.format(uct[i, j]) + '\t'
            ol += s
        ol += '          !CoordinatesOfUnitCell(UOL)\n'
        input_file_uc_block += ol

    ############################
    # Organize atoms and create input block with atom positions
    ############################
    atom_types = organize_atoms_by_types(newbulkats, nsl)
    input_file_atom_pos_block, input_file_mt_params_block = make_atoms_input_blocks(atom_types, uct)

    ##############################
    # Options block
    ##############################

    # l_max always set to 18 (not expensive)
    l_max = 18

    # energy range
    E2 = round(float(rp.THEO_ENERGIES[1]) + 20, 2) # add 20 eV to energy range
    E2_str = str(E2)
    Estep = round(float(rp.THEO_ENERGIES[2]), 2)
    Estep_str = str(Estep)
    n_cores_str = str(rp.N_CORES)  # number of threads

    input_options = ""
    input_options += "OPTIONS:" + '\n'
    input_options += "'s' !crystal: 'bulk'/'slab'." + '\n'
    input_options += "'n' !'yes'/'no': SpinUp&Down PhaseShifts calc?\"" + '\n'
    input_options += "'n' !'yes'/'no': Rho print?\"" + '\n'
    input_options += "'n' !'yes'/'no': Pot print?\"" + '\n'
    input_options += "'n' !'yes'/'no': WaveFunction print?\"" + '\n'
    input_options += " 000.00  " + E2_str + "    " + Estep_str + " !energy interval E1,E2,Estep" + '\n'
    input_options += "   " + n_cores_str + "   " + str(l_max) + "                !nthread,lmax" + '\n'
    input_options += "  1.d-06 1.d-09         !relerr abserr" + '\n'    # I hope these are reasonable values?
    input_options += " DIFFERENTIAL EVOLUTION METHOD" + '\n'
    # TODO check params below
    input_options += "    0.80   0.50         !F_XC,CR_XC" + '\n'       # does this need changing?
    input_options += "    0   1   0           !method" + '\n'           # what does this actually do?
    input_options += "    2   0.800000        !strategy,F_CR" + '\n'
    input_options += "  8000    0             !itermax"                 # fix itermax and other params?

    ##############################
    # Assemble input file
    ##############################
    input_file_txt = ""
    input_file_txt += input_file_header
    input_file_txt += input_file_length_line
    input_file_txt += input_file_uc_block
    input_file_txt += str(len(atom_types)) + '          ! # inequivalent atoms\n'
    input_file_txt += input_file_atom_pos_block
    input_file_txt += input_file_mt_params_block
    input_file_txt += input_options

    return atom_types, input_file_txt, l_max, E2, Estep


def organize_atoms_by_types(newbulkats, nsl):
    """
    Takes atoms from the slab (nsl) and groups them for the EEASISSS input. Atoms need to be organized into groups of
    same element and site. Note that a site could have mixed occupation.
    Returns atom_types, which contains the groups of atom types mapped to an index required for identification. Further,
    the input blocks containing the atom positions are produced here.
    """
    atom_types = {} # dict will contain Atom types
    for atom in nsl.atlist:
        new_bulk = True if atom in newbulkats else False
        if (atom.site, atom.el, new_bulk) not in atom_types.keys():
            atom_types[(atom.site, atom.el, new_bulk)] = tl.classes.sitetype.Atom_type(atom.el, str(atom.site), new_bulk)
            atom_types[(atom.site, atom.el, new_bulk)].add_atom(atom)
        else:
            atom_types[(atom.site, atom.el, new_bulk)].add_atom(atom)

    # Now simplify by remapping atom_types to a new dict with unique integer ID...
    re_map = {}
    for type_id, old_ref in enumerate(atom_types.keys(),
                                      1):  # old_ref is (sublayer_id, atom.site, atom.el, new_bulk)
        re_map[old_ref] = type_id
    atom_types = {re_map[id]: atom_types[id] for id in re_map}  # list comprehension magic

    return atom_types


def organize_atoms_by_sublayers(newbulkats, nsl):
    """
    Unused and outdated.
    Does same thing as organize_atoms_by_type but additionally groups them by sublayers as produced by
    slab.createSublayer(...)
    Produces worse results and EEASISSS takes much longer to execute.
    """

    # iterate over all types of atoms (i.e. all atoms that can be in each site types)
    #create the sublayers
    nsl.createSublayers(eps=0.001)
    atom_types = {}
    for sublayer_id, sublayer in enumerate(nsl.sublayers, 1):  # enumerate(..., 1) makes sublayer_id start at 1
        for atom in sublayer.atlist:
            new_bulk = True if atom in newbulkats else False
            if (sublayer_id, atom.site, atom.el, new_bulk) not in atom_types.keys():
                atom_types[(sublayer_id, atom.site, atom.el, new_bulk)] = \
                    tl.classes.sitetype.Atom_type(atom.el, str(atom.site), new_bulk)
            atom_types[(sublayer_id, atom.site, atom.el, new_bulk)].add_atom(atom)
    # We need to go through all this trouble of making Atom types, since we need to group the atoms into sublayers, but
    # we can't mix elements (which could otherwise happen if we have mixed occupation)
    # Now simplify by remapping atom_types to a new dict with unique integer ID...
    re_map = {}
    for type_id, old_ref in enumerate(atom_types.keys(), 1):  # old_ref is (sublayer_id, atom.site, atom.el, new_bulk)
        re_map[old_ref] = type_id
    atom_types = {re_map[id]: atom_types[id] for id in re_map}  # list comprehension magic

    return atom_types


def make_atoms_input_blocks(atom_types, uct):
    input_file_atom_pos_block = ""
    input_file_mt_params_block = ""
    for type_id in atom_types.keys():  # should be ordered
        atom_type = atom_types[type_id]
        type_header = '\t'.join([str(len(atom_type.atoms)), str(atom_type.get_atomic_number()), str(type_id),
                                 str(type_id), '{: .4f}'.format(atom_type.rmtmin),
                                 '{: .4f}'.format(atom_type.rmtmax), '{: .4f}'.format(atom_type.S), atom_type.el])
        if atom_type.new_bulk == True:  # for easy visibility
            type_header += '\t !atom in new bulk'
        input_file_atom_pos_block += type_header + '\n'
        atom_coords_table = ""
        for atom in atom_type.atoms:
            # rescale atomic coordinates from fractional to absolute
            # for this multiply fractional positon by unit cell matrix
            #fractional_vec = np.array([pos for pos in atom.pos])
            #position_vec = uct.dot(fractional_vec)

            atom_coords_table += ' '.join(['{: .5f}'.format(pos) for pos in atom.cartpos])
            atom_coords_table += '\n'
        input_file_atom_pos_block += atom_coords_table

        input_file_mt_params_block += ' '.join([str(type_id), '{: .4f}'.format(atom_type.rmtmin),
                                                '{: .4f}'.format(atom_type.rmtmax), '{: .4f}'.format(atom_type.S),
                                                '{: .4f}'.format(atom_type.fxc), atom_type.el])
        if type_id == 1:  # only in first line
            input_file_mt_params_block += '\t!iA rmtmin rmtmax rmtS fxc elem'
        elif type_id == 2:  # only in second line
            input_file_mt_params_block += '\t!etc ...'
        input_file_mt_params_block += '\n'
    return input_file_atom_pos_block, input_file_mt_params_block


def convert_eeasisss_output(sl, rp, atom_types, lmax, Emax, Estep, ps_outdir):
    ##############################
    # Read and convert EEASiSSS ouput
    ##############################
    """
    EEASiSSS produces a log and files called PS.xx.n.txt where xx is atomic number and n is an index. These files
    contain phase shifts in radiants for the atom types.
    """

    # go through all the files that were generated by EEASiSSS and read
    filelist = [filename for filename in os.listdir('./' + ps_outdir) if
                filename.startswith('PS.')]

    # Sanity check that output files are as expected:
    expected_files =[]
    for type_id, at_type in atom_types.items():
        expected_filename = 'PS.' + str(at_type.atomic_number) + '.' + str(type_id) + '.txt'
        expected_files.append(expected_filename)
    if set(filelist) == set(expected_files):
        pass  # should be the case...
    else:
        logger.error("Unexpected phaseshift files were generated.") # should never happen if input formatted right
        raise RuntimeError("Phaseshift generation failed.")

    ###################
    # Now we work on reshuffling
    sites_to_average = {} # likely an inelegant solution...
    for type_id, at_type in atom_types.items():
        if at_type.new_bulk:  # kick out all atoms that were added to the new bulk
            continue
        site_el = (at_type.label, at_type.el) # average phaseshifts over all atoms of same element in same site type
        if site_el not in sites_to_average.keys():
            sites_to_average[site_el] = []

        # lmax, E2 and Estep were read and set in format_eeasisss_input(...)
        filename = 'PS.' + str(at_type.atomic_number) + '.' + str(type_id) + '.txt'
        dir_prefix = './'+ ps_outdir +'/'
        # get number of lines: # could also be calculated from energy range
        num_lines = sum(1 for line in open(dir_prefix + filename, 'r', encoding='utf-8'))
        Emin = 0
        ps = np.zeros([num_lines, lmax])  # initialize np array with zeros
        energy = np.zeros([num_lines])
        with open(dir_prefix + filename, 'r', encoding='utf-8') as psfile:
            psfile.readline()  # firstline contains V0 parameters; discard read later
            for j, line in enumerate(psfile):
                values = line.split()
                energy[j] = values[0]
                ps[j, :] = values[1:lmax+1]

            sites_to_average[site_el].append(ps) # move the phaseshifts read from file to dict

    # Now average phaseshifts for each element and site:
    phaseshift_averages = {}
    for (label, el), values in sites_to_average.items():
        stack = np.stack(values, axis=-1)  # stack ps into np array
        ps_avg = np.mean(stack, axis=-1)  # take average; np has got to be the most efficient option.

        phaseshift_averages[(label, el)] = ps_avg

    # read inner potential coefficients from file; see Docstring and Wiki
    coef_filename = 'Vxc0EincAprx_v0coef'
    c0, c1, c2, c3 = read_V0_coefficients(os.path.join(ps_outdir,coef_filename))

    # old output format
    # sort first by element, then by site
    # first sort by POSCAR elements, same order as POSCAR:

    phaseshifts = []
    for j in range(int((Emax - Emin) / Estep) + 1):
        #old output format
        ps = []
        for el in sl.elements:
            if el in rp.ELEMENT_MIX:
                chemelList = rp.ELEMENT_MIX[el]
            else:
                chemelList = [el]
            # then by sites:
            siteList = [site for site in sl.sitelist if site.el == el]
            for cel in chemelList:
                for site in siteList:
                    label = cel + '_in_' + site.label
                    ps.append(phaseshift_averages[(label, cel)][j,:].tolist())

        energy_hartree= energy[j]/27.211396 # conversion from eV to Hartree
        phaseshifts.append([energy_hartree,ps])

    # format into old output format – int at beginning of line is skipped in old version too!
    formater = ff.FortranRecordWriter('(4F8.2)')
    firstline = formater.write([c0,c1,c2,c3])
    firstline += "\t" + "PS.r.00.00\txxxx 210930-000000"  # TODO clean?
    # copied from above, old formatting of first line... TODO
    if firstline == "":
        logger.error("Could not find first line for PHASESHIFTS file "
                     "(should contain MUFTIN parameters).")
        firstline = "ERROR: first line not found in EEASiSSS.x output\n"
        rp.setHaltingLevel(2)
    else:
        # add number of blocks to firstline
        nblocks = len(phaseshifts[0][1])
        firstline = str(nblocks).rjust(3) + " " + firstline # technically fails if more then 999 blocks - non-issue?
        # remove the "PS.r.**.**"
        firstline = re.sub(r"PS\.r\.[0-9]+\.[0-9]+", "", firstline)


    return firstline, phaseshifts


def read_V0_coefficients(coef_file_path='PS_out/Vxc0EincAprx_v0coef'):
    """
    Reads coefficients of the innerpotenial from file generated by EEASISSS. This used to be done by taking the first
    line of one of the PS.xx.xx files. While that is still possible, it is a) more elegant to read from the file
    intended for this purpose and b) the coefficients are written to that file with more significant digits.
    Attention: Order of coefficients has changed in the 2021 EEASISSS version. (0123 <-> 1230)
    """
    try:
        #open file
        with open(coef_file_path) as coef_file:
            # file should just be one line with content:
            # v0coef = a b c d
            # where a,b,c,d, are floats coefficient values
            line = coef_file.readline()
            c1, c2, c3, c0 = line.split()[2:6]
        (c1, c2, c3, c0) = [float(c) for c in (c1, c2, c3, c0)]
    except Exception:
        logger.error("Could not open and read file " + coef_file_path)
        raise RuntimeError("Unable to read inner potential coefficients")
    return c0,c1,c2,c3


def move_EEASISSS_files(ps_outdir, log_filename, Vxc0_files = ['Vxc0Einc', 'Vxc0EincAprx', 'Vxc0EincAprx_v0coef']):
    """
    Moves files from PS_out folder to work folder and then removes PS_out.
    """
    outdir_path = os.path.join(".", ps_outdir + "/")
    # move log file
    shutil.move(os.path.join(outdir_path, log_filename), os.path.join(".", log_filename))
    # move muffin tin potential files
    for file in Vxc0_files:
        shutil.move(os.path.join(outdir_path, file), os.path.join(".", file))
    # only remaining files are individual phaseshift files; can be discarded since they are combines in PHASESHIFTS
    # remove PS_out folder with contents
    shutil.rmtree(outdir_path)

    # remove files ulog*, udat* and uinp* created by eeasisss programme
    files_in_directory = os.listdir(".")
    filtered_files = [file for file in files_in_directory if (file.startswith('ulog')
                                                     or file.startswith('udat') or file.startswith('uinp'))]
    for file in filtered_files:
        os.remove(os.path.join(".", file))
