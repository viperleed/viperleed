# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 15:32:12 2020

@author: Florian Kraushofer
"""

import os
import copy
import numpy as np
import logging
import random
import subprocess
import re

import tleedmlib as tl

logger = logging.getLogger("tleedm.psgen")

def runPhaseshiftGen(sl,rp, psgensource=os.path.join('source','EEASiSSS.x'),
                            excosource=os.path.join('source','seSernelius'),
                            atdenssource=os.path.join('source', 
                                                      'atom_density_files')):
    """Creates required input for EEASiSSS.x, then runs it. Reads the output 
    files and extracts information for _PHASESHIFTS file, then returns that 
    information (without writing _PHASESHIFTS)."""
    psgensource = os.path.join(rp.workdir, psgensource)
    excosource = os.path.join(rp.workdir, excosource)
    atdenssource = os.path.join(rp.workdir, atdenssource)
    
    lmax = 16   # this could be a variable, for now set fixed...

    nsl = copy.deepcopy(sl)
    nsl.getFractionalCoordinates()
    # bulk calculation:
    blayers = [l for l in nsl.layers if l.isBulk]
    # construct bulk slab
    # bsl = copy.deepcopy(nsl)    #this makes the site different objects... 
    # bsl.atlist = [at for at in bsl.atlist if at.layer.isBulk]
    # bsl.layers = [l for l in bsl.layers if l.isBulk]
    
    if type(rp.BULK_REPEAT) == np.ndarray:
        bulkc = rp.BULK_REPEAT
        if bulkc[2] < 0:
            bulkc = -bulkc
    else:
        cvec = nsl.ucell[:,2]
        if rp.BULK_REPEAT is None:
            # assume that interlayer vector from bottom non-bulk to top 
            #  bulk layer is the same as between bulk units
            zdiff = (blayers[-1].cartbotz 
                     - nsl.layers[blayers[0].num-1].cartbotz)
        elif type(rp.BULK_REPEAT) == float:
            zdiff = rp.BULK_REPEAT
        bulkc = cvec * zdiff / cvec[2]

    # bsl.ucell[:,2] = bulkc
    # bsl.collapseCartesianCoordinates()
    
    # add a bulk layer to the slab
    nsl.getCartesianCoordinates()
    cfact = (nsl.ucell[2,2] + abs(bulkc[2])) / nsl.ucell[2,2]
    nsl.ucell[:,2] = nsl.ucell[:,2] * cfact
    nsl.getFractionalCoordinates()
    newbulkats = []
    tmplist = nsl.atlist[:]
    for at in tmplist:
        if at.layer.isBulk: 
            newbulkats.append(at.duplicate())
        at.cartpos = at.cartpos - bulkc
    nsl.collapseCartesianCoordinates(updateOrigin=True)
    nsl.sortByZ()
    
    outvals = {}    # dict containing lists of values to output. 
                    #   format outvals[energy][block][L]
                    
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
    blocks = []  #tuples of sites (with poscar elements) and elements 
        #  (same as POSCAR if no ELEMENT_MIX, ELEMENT_MIX elements if not)
    for site in nsl.sitelist:
        if site.el in rp.ELEMENT_MIX:
            for el in rp.ELEMENT_MIX[site.el]:
                blocks.append((site,el))
        else:
            blocks.append((site,site.el))
    scsize = 1 
    if len(rp.ELEMENT_MIX) > 0:
        minnum = -1
        for (site, el) in [(site, el) for (site, el) in blocks if site.el 
                          in rp.ELEMENT_MIX and site.occ[el] > 0.]:
            al = [at for at in nsl.atlist if at.site == site]
            atcount = len(al)*site.occ[el]
            if minnum < 0 or (minnum > atcount > 0):
                minnum = atcount
        # we want at least 2 atoms of each element in each site type:
        if 0 < minnum < 2.0:
            scsize = int(np.ceil(2.0/minnum))
    if scsize > 1:  # some checks to make sure it doesn't get too large
        maxcells = 20  # upper limit on supercell size
        maxats = 500   # upper limit on atoms in supercell
        if scsize > maxcells:
            scsize = maxcells
            # don't warn - this is a large unit cell either way.
        if len(nsl.atlist) * scsize > maxats:
            logger.debug("Phaseshift generation: Given element "
                "concentrations would require a very large supercell. "
                "Element concentrations for low-occupancy elements will be "
                "increased to avoid this. This only concerns the phaseshifts "
                "calculation and should not cause problems.")
            # determine minimum size to have 2 of each element
            minsize = 1
            for site in [s for s in nsl.sitelist if s.el in rp.ELEMENT_MIX]:
                ats = len([at for at in nsl.atlist if at.site == site])
                els = len([el for el in rp.ELEMENT_MIX[site.el] if 
                                   site.occ[el] > 0.])
                minsize = max(minsize, int(np.ceil(2*els / ats)))
            scsize = max(minsize, int(maxats / len(nsl.atlist)))
                
    subatlists = {}     # atlist per block tuple
    if scsize > 1:  #construct supercell to get enough atoms
        xsize = int(np.ceil(np.sqrt(scsize)))  # if scsize is not prime, try 
        while scsize % xsize != 0:             #   making it close to square 
            xsize += 1
        ysize = int(scsize / xsize)
        cpatlist = nsl.atlist[:]
        for at in cpatlist:
            for i in range(0,xsize):
                for j in range(0, ysize):
                    if i == j == 0:
                        continue
                    tmpat = at.duplicate()
                    tmpat.pos[0] += i
                    tmpat.pos[1] += j
        nsl.getCartesianCoordinates()
        nsl.ucell = np.dot(np.array([[xsize,0,0],[0,ysize,0],[0,0,1]]),
                           nsl.ucell)
        nsl.getFractionalCoordinates()
      
    for site in nsl.sitelist:
        if site.el in rp.ELEMENT_MIX:
            occdict = {}
            for (k, v) in site.occ.items(): 
                if v > 0.0: occdict[k] = v
            occdict = dict(sorted(occdict.items(), 
                                  key = lambda kv:(kv[1],kv[0])))  
                                    #sort by occupancy values
            al = [at for at in nsl.atlist if at.site == site]
            totats = len(al)
            for el in occdict:
                subatlists[(site,el)] = []
                reqats = int(np.ceil(totats * site.occ[el]))
                reqats = max(2, reqats)
                while reqats > 0 and len(al) > 0:
                    at = random.choice(al)
                    subatlists[(site,el)].append(at)
                    al.remove(at)
                    reqats -= 1
            if len(al) > 0: #should never happen
                logger.warning("Error in _PHASESHIFTS file "
                        "generation: Not all atoms were distributed!")
        else:
            subatlists[(site,site.el)] = [at for at in nsl.atlist 
                                          if at.site == site]
    blocks = [(site,el) for (site,el) in blocks 
              if len(subatlists[(site,el)]) > 0]

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
    for i in range(0,3):
        ol = ''
        for j in range(0,3):
            s = str(round(uct[i,j],4))+' '
            ol += s.ljust(8)
        if i == 0:
            ol += '      CoordinatesOfUnitCell(LCunits)\n'
        else:
            ol += '\n'
        output += ol

    output += (str(len(nsl.atlist))+"  "+str(len(nsl.atlist))
               +"                  #AtomTypes,#OccupiedAtomTypes\n")
    ptl = [el.lower() for el in tl.leedbase.periodic_table]

    chemels = {}
    chemelspaths = {}
    for (site, el) in blocks:
        if el in rp.ELEMENT_RENAME:
            chemel = rp.ELEMENT_RENAME[el]
        elif el.lower() in ptl:
            chemel = el.capitalize()
        else:
            logger.error("Error generating _PHASESHIFTS file: Could not "
                          "identify "+el+" as a chemical element. Define "
                          "ELEMENT_RENAME or ELEMENT_MIX parameter.")
            raise
        chgdenrelpath = os.path.join(atdenssource,chemel,"chgden"+chemel)
        if os.name == 'nt':     #windows - replace the backslashes.
            chgdenrelpath = chgdenrelpath.replace('/','\\')
        chemels[el] = chemel
        chemelspaths[el] = chgdenrelpath
    
    nsl.sortByZ(botToTop=True)
    for at in nsl.atlist:
        # realcartpos = np.dot(nsl.ucell, at.pos)   
              # use the "real" cartesian system, with Z going up
        for (site,el) in blocks:
            if at in subatlists[(site,el)]:
                chemel = chemels[el]
                chgdenpath = chemelspaths[el]
        output += ("1 "+str(tl.leedbase.periodic_table.index(chemel)+1)
                   + ".  0.  0.  '"+chgdenpath+"'\n")
        ol = ""
        for j in range(0,3):
            # ol += str(round(realcartpos[j],4))+" "
            ol += str(round(at.cartpos[j],4))+" "
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
    output += ("0.0," + str(round(float(rp.THEO_ENERGIES[1]) + 20, 1))
               +","+str(round(float(rp.THEO_ENERGIES[2]),1))
               +" |'phaseshift'/'dataflow' run: E1->E2,PS_Estep\n")
    output += "100.,100.,1  |'sigma' run: E1->E2,NumVals\n"
    # TODO - check the following, do they need to be changed?
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
    except:
        logger.error("Phaseshift data generation: Failed to write "
                      +outfilename+". Proceeding with execution...")
    if os.name == 'nt':
        logger.error("Phaseshift generation is currently not "
                         "supported on Windows. Use a linux shell to run "
                         "phaseshift generation script.")
        raise EnvironmentError("Phaseshift generation is currently not "
                               "supported on Windows.")
    else:
        subprocess.run(psgensource,input=output,encoding='ascii')
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
                except:
                    remlist.append(filename)
    for filename in remlist:
        filelist.remove(filename)
    filelist.sort(key=lambda filename: int(filename.split('.')[-1]))
                                                #now sorted by atom number
    firstline = ""
    atoms_phaseshifts = [[] for i in range(0,len(filelist))]    
                                                #data from all the PS files                                               
    for (i,filename) in enumerate(filelist):
#        if bulk or wsl.atlist[i].site not in bulksites: 
        if nsl.atlist[i] not in newbulkats:
            psfile = open(filename, 'r')
            reade = -1.
            pslist = []
            ps_of_e = {}    #dictionary of pslist for energy, where pslist 
                            # is a list of floats as read from the PS-file
            for (j,line) in enumerate(psfile):
                if j == 0:
#                    if not bulk: 
                    if firstline == "":
                        firstline = line[2:] #ignore I2 at beginning
                else:   #line should contain data
                    values = [float(s) for s in line.split()]
                    if reade < 0 or values[0] == reade + rp.THEO_ENERGIES[2]:     
                        #new energy value, start phaseshift value list
                        if reade >= 0:
                            ps_of_e[reade] = pslist
                        reade = values[0]
                        pslist = values[1:]
                    else:
                        pslist.extend(values)
            ps_of_e[reade] = pslist
            psfile.close()
            atoms_phaseshifts[i] = ps_of_e
        os.remove(os.path.join('.',filename))   #delete files
    for (site,el) in blocks:
        writeblock = False
        pssum = None
        for (i,at) in enumerate(nsl.atlist):
            if at not in newbulkats:
#            if bulk or at.site not in bulksites: 
                if at in subatlists[(site,el)]:
                    if pssum is None:
                        writeblock = True
                        pel = at.el #POSCAR element for block
                        pssum = atoms_phaseshifts[i]
                        n = 1
                    else:
                        for en in atoms_phaseshifts[i]:
                            for j in range(0,
                                           len(atoms_phaseshifts[i][en])):
                                pssum[en][j] += atoms_phaseshifts[i][en][j]
                        n += 1
        if writeblock:
            # append the values for the whole block to outvals:
            for en in pssum:
                pssum[en] = [v/n for v in pssum[en]]
                if not en in outvals: 
                    outvals[en] = []
                outvals[en].append((pel,el,site,pssum[en]))
    # clean up
#    bss = "-bulk" if bulk else "-slab"
    bss = ""
    try:
        os.rename('logfile','eeasisss-logfile'+bss)
    except:
        logger.warning("Failed to rename phaseshift generation file "
                        "'logfile' to 'eeasisss-logfile"+bss+"'")
    try:
        os.rename('QMTvsE','eeasisss-QMTvsE'+bss)
    except:
        logger.warning("Failed to rename phaseshift generation file "
                        "'QMTvsE' to 'eeasisss-QMTvsE"+bss+"'")
    try:
        os.rename('RMTvsE','eeasisss-RMTvsE'+bss)
    except:
        logger.warning("Failed to rename phaseshift generation file "
                        "'RMTvsE' to 'eeasisss-RMTvsE"+bss+"'")
    try:
        os.rename('V0vsE','eeasisss-V0vsE'+bss)
    except:
        logger.warning("Failed to rename phaseshift generation file "
                        "'V0vsE' to 'eeasisss-V0vsE"+bss+"'")
    # sort blocks in outvals:
    outvalsSorted = {}
    outvalLength = 0    # Due to a bug in eeasisss, it does not generate 
            # exactly the energies it is asked to; this can result in 
            # different energies being present for the bulk and the slab
            # calculations. Workaround: Discard energies where not all sites 
            # are present.
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
        siteList = [site for site in sl.sitelist if site.el == el]
        for cel in chemelList:
            for site in siteList:
                for en in outvalsSorted:
                    for (o_el,o_cel,o_site,pslist) in outvals[en]:
                        if (o_el == el and o_cel == cel 
                              and o_site.isEquivalent(site)):
                            outvalsSorted[en].append(pslist)
    # return:
    # writePHASESHIFTS wants values as list of tuples and energies in Hartree:
    phaseshifts = []
    for en in outvalsSorted:
        if len(outvalsSorted[en]) == outvalLength:    # drop energies where
                        # phaseshift was not calculated for all sites
            phaseshifts.append((en/27.2116,outvalsSorted[en]))
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
        firstline = re.sub("PS\.r\.[0-9]+\.[0-9]+", "", firstline)
    
    return (firstline, phaseshifts)