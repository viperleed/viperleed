# -*- coding: utf-8 -*-
"""
Discarded code

"""

def runPhaseshiftGen(sl,rp,psgensource=os.path.join('.','source','EEASiSSS.x'),
                     excosource=os.path.join('.','source','seSernelius'),
                     atdenssource=os.path.join('.','source',
                                               'atom_density_files')):
    lmax = 18   # this could be a variable, for now set fixed...
#    dirpath = os.getcwd()
#    if os.name == 'nt':     #windows - replace the stupid backslashes...
#        dirpath = dirpath.replace('\\','/')
#    folders = dirpath.split('/')
#    if len(folders) > 0:
#        foldername = folders[-1]
#    else:
#        foldername = 'unknown'
    foldername = os.path.basename(os.getcwd())

    nsl = copy.deepcopy(sl)
    nsl.getFractionalCoordinates()
    # bulk calculation:
    blayers = [l for l in nsl.layers if l.isBulk]
    # construct bulk slab
    bsl = copy.deepcopy(nsl)    #this makes the site different objects... 
    bsl.atlist = [at for at in bsl.atlist if at.layer.isBulk]
    newsitelist = [] # ...make site objects the same again
    for site in bsl.sitelist:
        for site2 in nsl.sitelist:
            if site.isEquivalent(site2): newsitelist.append(site2)
    bsl.sitelist = newsitelist
    for at in bsl.atlist:
        for site in bsl.sitelist:
            if at.site.isEquivalent(site): at.site = site
    bsl.layers = [l for l in bsl.layers if l.isBulk]
    
    cvec = bsl.ucell[:,2]
    if rp.ASAZ == tl.DEFAULT:
        # assume that interlayer vector from bottom non-bulk to top bulk layer 
        #  is the same as between bulk units
        zdiff = blayers[-1].cartbotz - nsl.layers[blayers[0].num-1].cartbotz
    else:
        zdiff = blayers[-1].cartbotz - blayers[0].carttopz + rp.ASAZ
    bulkc = cvec * zdiff / cvec[2]
    bsl.ucell[:,2] = bulkc
    bsl.collapseCartesianCoordinates()
    
    # add a bulk layer to the slab
    nsl.getCartesianCoordinates()
    nsl.ucell[:,2] = nsl.ucell[:,2] + bulkc
    nsl.getFractionalCoordinates()
    coff = zdiff / cvec[2]
    newbulkats = []
    tmplist = nsl.atlist[:]
    for at in tmplist:
        if at.layer.isBulk: 
            newbulkats.append(at.duplicate())
        at.pos[2] += coff
    nsl.sortByZ()
    nsl.getCartesianCoordinates()
    
    outvals = {}    # dict containing lists of values to output. 
                    #   format outvals[energy][block][L]
                    
    # The following originally calculated phaseshifts first for only the bulk,
    #   then for the slab, taking phaseshifts from the bulk calculations if the
    #   site is present at all in the bulk. However, this means that different
    #   muffin tin parameters are applied for the different phaseshifts. Now,
    #   take phaseshifts from the slab calculation only. Average only over 
    #   atoms NOT added as "new bulk".
    for bulk in [True,False]:
        wsl = bsl if bulk else nsl  # working slab for this iteration
        # before starting on unit cell, determine whether supercell is needed:
        blocks = []  #tuples of sites (with poscar elements) and elements 
                     #  (same as POSCAR if no ELSPLIT, ELSPLIT elements if not)
        for site in wsl.sitelist:
            if site.el in rp.ELSPLIT:
                for el in rp.ELSPLIT[site.el]:
                    blocks.append((site,el))
            else:
                blocks.append((site,site.el))
        scsize = 1 
        if len(rp.ELSPLIT) > 0:
            minnum = -1
            for (site, el) in [(site, el) for (site, el) in blocks if site.el 
                              in rp.ELSPLIT and site.occ[el] > 0.0]:
                al = [at for at in wsl.atlist if at.site == site]
                atcount = len(al)*site.occ[el]
                if minnum < 0 or (minnum > atcount > 0):
                    minnum = atcount
            # we want at least 2 atoms of each element in each site type:
            if 0 < minnum < 2.0:
                scsize = int(np.ceil(2.0/minnum))
        subatlists = {}     # atlist per block tuple
        if scsize > 1:  #construct N x 1 supercell to get enough atoms
            cpatlist = wsl.atlist[:]
            for at in cpatlist:
                for i in range(1,scsize):
                    tmpat = at.duplicate()
                    tmpat.pos[0] += i
            wsl.getCartesianCoordinates()
            wsl.ucell = np.dot(np.array([[scsize,0,0],[0,1,0],[0,0,1]]),
                               wsl.ucell)
            wsl.getFractionalCoordinates()
          
        for site in wsl.sitelist:
            if site.el in rp.ELSPLIT:
                occdict = {}
                for (k, v) in site.occ.items(): 
                    if v > 0.0: occdict[k] = v
                occdict = dict(sorted(occdict.items(), 
                                      key = lambda kv:(kv[1],kv[0])))  
                                        #sort by occupancy values
                al = [at for at in wsl.atlist if at.site == site]
                totats = len(al)
                for el in occdict:
                    subatlists[(site,el)] = []
                    reqats = int(np.ceil(totats * site.occ[el]))
                    if reqats < 2: reqats = 2
                    while reqats > 0 and len(al) > 0:
                        at = random.choice(al)
                        subatlists[(site,el)].append(at)
                        al.remove(at)
                        reqats -= 1
                if len(al) > 0: #should never happen
                    logging.warning("Error in _PHASESHIFTS file "
                            "generation: Not all atoms were distributed!")
            else:
                subatlists[(site,site.el)] = [at for at in wsl.atlist 
                                              if at.site == site]
        blocks = [(site,el) for (site,el) in blocks 
                  if len(subatlists[(site,el)]) > 0]

        if bulk: 
            bulksites = [site for (site,el) in blocks]
        else:
            # site objects are different in the new slab; copy to equivalent
            oldbulksites = bulksites[:]
            bulksites = []
            for oldsite in oldbulksites:
                for newsite in [site for (site,el) in blocks]:
                    if newsite.isEquivalent(oldsite):
                        bulksites.append(newsite)
        # start writing output, will be input for EEASiSSS code:
        output = ''
        output += "STRUCTURE:\n"
        output += foldername+" "+time.strftime("%Y-%m-%d %H:%M:%S", 
                                               time.gmtime())+"\n"
        if bulk:
            output += ("'b'  1.00 16      BulkOrSlab('b' or 's'), "
                      "LatticeConstant(Angstroms), nshell\n")
        else:
            output += ("'s'  1.00 16      BulkOrSlab('b' or 's'), "
                      "LatticeConstant(Angstroms), nshell\n")
        uct = wsl.ucell.transpose()
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

        output += (str(len(wsl.atlist))+"  "+str(len(wsl.atlist))
                   +"                  #AtomTypes,#OccupiedAtomTypes\n")
        ptl = [el.lower() for el in tl.periodic_table]
    
        chemels = {}
        chemelspaths = {}
        for (site, el) in blocks:
            if el in rp.ELDEF:
                chemel = rp.ELDEF[el]
            elif el.lower() in ptl:
                chemel = el.lower().capitalize()
            else:
                logging.error("Error generating _PHASESHIFTS file: Could not "
                              "identify "+el+" as a chemical element. Define "
                              "ELDEF or ELSPLIT parameter.")
                raise
            chgdenrelpath = os.path.join(atdenssource,chemel,"chgden"+chemel)
            if os.name == 'nt':     #windows - replace the backslashes.
                chgdenrelpath = chgdenrelpath.replace('/','\\')
            chemels[el] = chemel
            chemelspaths[el] = chgdenrelpath
        
        wsl.sortByZ()
        for at in wsl.atlist:
            realcartpos = np.dot(wsl.ucell, at.pos)   
              # use the "real" cartesian system, with Z going up
            for (site,el) in blocks:
                if at in subatlists[(site,el)]:
                    chemel = chemels[el]
                    chgdenpath = chemelspaths[el]
            output += ("1 "+str(tl.periodic_table.index(chemel)+1)
                       + ".  0.  0.  '"+chgdenpath+"'\n")
            ol = ""
            for j in range(0,3):
                ol += str(round(realcartpos[j],4))+" "
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
        output += (str(round(float(rp.TEOENRANGE[0])-4*rp.TEOENRANGE[2],1))
                  +","+str(round(float(rp.TEOENRANGE[1])+2*rp.TEOENRANGE[2],1))
                  +","+str(round(float(rp.TEOENRANGE[2]),1))
                  +" |'phaseshift'/'dataflow' run: E1->E2,PS_Estep\n")
        output += "100.,100.,1  |'sigma' run: E1->E2,NumVals\n"
        # TODO - check the following, do they need to be changed?
        output += "MT OPTIMIZATION (Nelder-Mead):\n"
        output += "'*'           | NMselect: '*' or '+'\n"
        output += "0.05d0   | NMlambda, simplex point shift.\n"
        output += "3            | NMiter, no. of simplex calls.\n"
        output += "1.d-04    ! NMeps,  simplex epsilon.\n"
        output += "0.125d0 | NMepsit, epsilon iteration, eps=eps*epsit.\n"
        outfilename = 'eeasisss-input-bulk' if bulk else 'eeasisss-input-slab'
        try:
            with open(outfilename, 'w') as wf:
                wf.write(output)
        except:
            logging.error("Phaseshift data generation: Failed to write "
                          +outfilename+". Proceeding with execution...")
        if os.name == 'nt':
            logging.critical("Phaseshift generation is currently not "
                             "supported on Windows. Use a linux shell to run "
                             "phaseshift generation script.")
            raise EnvironmentError("Phaseshift generation is currently not "
                                   "supported on Windows.")
        else:
            subprocess.run(psgensource,input=output,encoding='ascii')
        # go through all the files that were generated by EEASiSSS and read
        filelist = [filename for filename in os.listdir('.') if 
                    filename.startswith('PS.r.')]
        for filename in filelist:
            try: 
                int(filename.split('.')[-1])
            except:
                filelist.remove(filename)
        filelist.sort(key=lambda filename: int(filename.split('.')[-1]))    
                                                    #now sorted by atom number
        firstline = ""
        atoms_phaseshifts = [[] for i in range(0,len(filelist))]    
                                                    #data from all the PS files                                               
        for (i,filename) in enumerate(filelist):
            if bulk or wsl.atlist[i].site not in bulksites: 
                psfile = open(filename, 'r')
                reade = -1.
                pslist = []
                ps_of_e = {}    #dictionary of pslist for energy, where pslist 
                                # is a list of floats as read from the PS-file
                for (j,line) in enumerate(psfile):
                    if j == 0:
                        if not bulk: 
                            if firstline == "":
                                firstline = line[2:] #ignore I2 at beginning
                    else:   #line should contain data
                        values = [float(s) for s in tl.linelist(line)]
                        if reade < 0 or values[0] == reade + rp.TEOENRANGE[2]:     
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
            pssum = tl.DEFAULT
            for (i,at) in enumerate(wsl.atlist):
                if bulk or at.site not in bulksites: 
                    if at in subatlists[(site,el)]:
                        if pssum == tl.DEFAULT:
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
        bss = "bulk" if bulk else "slab"
        try:
            os.rename('logfile','eeasisss-logfile-'+bss)
        except:
            logging.warning("Failed to rename phaseshift generation file "
                            "'logfile' to 'eeasisss-logfile-"+bss+"'")
        try:
            os.rename('QMTvsE','eeasisss-QMTvsE-'+bss)
        except:
            logging.warning("Failed to rename phaseshift generation file "
                            "'QMTvsE' to 'eeasisss-QMTvsE-"+bss+"'")
        try:
            os.rename('RMTvsE','eeasisss-RMTvsE-'+bss)
        except:
            logging.warning("Failed to rename phaseshift generation file "
                            "'RMTvsE' to 'eeasisss-RMTvsE-"+bss+"'")
        try:
            os.rename('V0vsE','eeasisss-V0vsE-'+bss)
        except:
            logging.warning("Failed to rename phaseshift generation file "
                            "'V0vsE' to 'eeasisss-V0vsE-"+bss+"'")
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
        if el in rp.ELSPLIT:
            chemelList = rp.ELSPLIT[el]
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
    # append zeroes for vacancy site:
    for en in outvals:
        outvalsSorted[en].append([0.0]*len(outvalsSorted[en][0]))
    # return:
    # writePHASESHIFTS wants values as list of tuples and energies in Hartree:
    phaseshifts = []
    for en in outvalsSorted:
        if len(outvalsSorted[en]) == outvalLength+1:    # drop energies where
                        # phaseshift was not calculated for all sites
            phaseshifts.append((en/27.2116,outvalsSorted[en]))
    if firstline == "":
        logging.critical("Could not find first line for PHASESHIFTS file "
                         "(should contain MUFTIN parameters).")
        firstline = "ERROR: first line not found in EEASiSSS.x output\n"
        rp.halt = 2
    else:
        # add number of blocks to firstline
        nblocks = len(phaseshifts[0][1])
        firstline = str(nblocks).rjust(3)+"  "+firstline
    
    return (firstline, phaseshifts)


def reduceUnitCell_OLD(ab, eps = 0.0001):
    """Takes an obtuse unit cell as a (2x2) matrix and reduces it to minimum 
    circumference, keeping the area constant. This might reduce oblique unit 
    cells to rectangular or hexagonal ones. Returns the modified unit cell, 
    as well as the transformation matrix."""
    u = np.array([[1,0],[0,1]])
    while abs(np.vdot(ab[0],ab[1])) > eps and np.vdot(ab[0],ab[1]) < 0:
        if np.linalg.norm(ab[0]) > np.linalg.norm(ab[1]):
            m = np.array([[1,1],[0,1]])
        else:
            m = np.array([[1,0],[1,1]])
        tmp = np.dot(m,ab)
        if np.vdot(tmp[0],tmp[1]) <= 0 or abs(np.vdot(tmp[0],tmp[1])) < eps:
            u = np.dot(m, u)
            ab = np.dot(m, ab)
        else:
            break
    return ab, u

def collapseCartesianCoordinatesOld(self):
    """Finds atoms outside the parallelogram spanned by the unit vectors a and b and moves them inside."""
    ab = self.ucell[0:2,0:2]
    angleAB = tl.angle(ab[:,0],ab[:,1])
    for at in self.atlist:
        xyCart = at.cartpos[0:2]
        for i in range(0,2):
            proj = (np.dot(xyCart,ab[:,i]) / np.linalg.norm(ab[:,i])**2) * ab[:,i]   #projection to unit vector parallel to the other unit vector
            if angleAB < np.pi/2:    #acute
                projlen = ((np.dot(proj,ab[:,i])/(np.linalg.norm(ab[:,i])))-(np.linalg.norm(xyCart-proj)/np.tan(angleAB))) / np.linalg.norm(ab[:,i])
            else:   #obtuse         there might be a better solution for this than if/else, but this is what I came up with...
                projlen = ((np.dot(proj,ab[:,i])/(np.linalg.norm(ab[:,i])))+(np.linalg.norm(xyCart-proj)/np.tan(np.pi-angleAB))) / np.linalg.norm(ab[:,i])
            while projlen >= 1 or projlen < 0:
                if projlen >= 1:
                    at.cartpos[0:2] -= ab[:,i]
                else:
                    at.cartpos[0:2] += ab[:,i]
                proj = (np.dot(xyCart,ab[:,i]) / np.linalg.norm(ab[:,i])**2) * ab[:,i]
                if angleAB < np.pi/2:    #acute
                    projlen = ((np.dot(proj,ab[:,i])/(np.linalg.norm(ab[:,i])))-(np.linalg.norm(xyCart-proj)/np.tan(angleAB))) / np.linalg.norm(ab[:,i])
                else:   #obtuse 
                    projlen = ((np.dot(proj,ab[:,i])/(np.linalg.norm(ab[:,i])))+(np.linalg.norm(xyCart-proj)/np.tan(np.pi-angleAB))) / np.linalg.norm(ab[:,i])
 

def getFractionalCoordinatesOld(self):
    """Calculates fractional coordinates for all atoms from their cartesian coordinates, using the slab unit cell."""
    self.collapseCartesianCoordinates()
    angleC = tl.angle(self.ucell[:,2],[0,0,1])
    ab = self.ucell[0:2,0:2]
    angleAB = tl.angle(ab[:,0],ab[:,1])
    for at in self.atlist:
        tp = np.copy(at.cartpos)    #temp position, will be simplified on the way
        tp[2] = self.topatOriZ-tp[2] 
        at.pos[2] = (tp[2]/np.cos(angleC))/np.linalg.norm(self.ucell[:,2])
        tp -= at.pos[2]*self.ucell[:,2]     # c contribution is done now
        xyCart = tp[0:2]
        for i in range(0,2):
            proj = (np.dot(xyCart,ab[:,i]) / (np.linalg.norm(ab[:,i])**2)) * ab[:,i]
            at.pos[i] = ((np.dot(proj,ab[:,i])/(np.linalg.norm(ab[:,i])))-(np.linalg.norm(xyCart-proj)/np.tan(angleAB))) / np.linalg.norm(ab[:,i]) % 1.0



pool = mp.Pool(mp.cpu_count())  #finding the candidate positions can be time expensive, but generating sublists per atom can be parallelized
result = pool.starmap_async(getAtomSymposlist, [(lowocclayer.atlist,i,eps,celltype) for i in range(0,len(lowocclayer.atlist))])
pool.close()
pool.join()
resultlist = result.get()
symposlist = [np.array([0.0,0.0])]    #always check the origin, even if it was not found.
for sublist in resultlist:
    symposlist = tl.addUnequalPoints(sublist, symposlist, eps, uniqueLists=True)
    
    
def getAtomSymposlist(atlist,i,eps,celltype):
    # TODO: For slabs with many atoms per layer, this is a bottleneck for speed in symmetry detection. Could be implemented as a compiled C routine
    result = []
    at = atlist[i]
    tmplist = []
    threshold = 100
    tmplist.append(at.cartpos[0:2])
    for (j, atj) in [(j, atj) for (j, atj) in enumerate(atlist) if j>i]:
        tmplist.append((at.cartpos[0:2]+atj.cartpos[0:2])/2)
        if len(tmplist) >= threshold:
            result = tl.addUnequalPoints(tmplist, result, eps)
            threshold = len(result)
            tmplist = []
        if celltype == "hexagonal":
            for (k, atk) in [(k, atk) for (k, atk) in enumerate(atlist) if k>j]:
                tmplist.append((at.cartpos[0:2]+atj.cartpos[0:2]+atk.cartpos[0:2])/3)
                if len(tmplist) >= threshold:
                    result = tl.addUnequalPoints(tmplist, result, eps)
                    threshold = len(result)
                    tmplist = []
    if len(tmplist) > 0: result = tl.addUnequalPoints(tmplist,result,eps)
    return result


############################################
## fom generating candidate positions:    ##
############################################


#            for (i, ati) in enumerate(sl.atlist):
##                p = ati.cartpos[0:2]
##                if not tl.pointIsInList(p,sl.symposlist,eps): sl.symposlist.append(p)
#                endtime = timer()
#                logging.debug("Time elapsed: "+str(endtime-starttime))
#                starttime=endtime
#                logging.debug("Atom "+str(i))
#                tmplist.append(ati.cartpos[0:2])
#                if len(tmplist) >= threshold:
#                    sl.symposlist = tl.addUnequalPoints(tmplist, sl.symposlist, eps)
#                    threshold = len(sl.symposlist)
#                    tmplist = []
#                for (j, atj) in [(j, atj) for (j, atj) in enumerate(sl.atlist) if j>i]:
##                    p = (ati.cartpos[0:2]+atj.cartpos[0:2])/2
##                    if not tl.pointIsInList(p,sl.symposlist,eps): sl.symposlist.append(p)
#                    tmplist.append((ati.cartpos[0:2]+atj.cartpos[0:2])/2)
#                    if len(tmplist) >= threshold:
#                        sl.symposlist = tl.addUnequalPoints(tmplist, sl.symposlist, eps)
#                        threshold = len(sl.symposlist)
#                        tmplist = []
#                    if celltype == "hexagonal":
#                        for (k, atk) in [(k, atk) for (k, atk) in enumerate(sl.atlist) if k>j]:
##                            p = (ati.cartpos[0:2]+atj.cartpos[0:2]+atk.cartpos[0:2])/3
##                            if not tl.pointIsInList(p,sl.symposlist,eps): sl.symposlist.append((ati.cartpos[0:2]+atj.cartpos[0:2]+atk.cartpos[0:2])/3)
#                            tmplist.append((ati.cartpos[0:2]+atj.cartpos[0:2]+atk.cartpos[0:2])/3)
#                            if len(tmplist) >= threshold:
#                                sl.symposlist = tl.addUnequalPoints(tmplist, sl.symposlist, eps)
#                                threshold = len(sl.symposlist)
#                                tmplist = []
#            if len(tmplist) > 0: sl.symposlist = tl.addUnequalPoints(tmplist,sl.symposlist,eps)
#            # weed out duplicates
#            logging.debug("Looking for duplicates...")
#            logging.debug("Positions to check: "+str(len(sl.symposlist)))
            
#            newlist = []
#            for (i,posi) in enumerate(sl.symposlist):
#                found = False
#                if i % 100 == 0: logging.debug(i)
#                for (j,pos) in [(j, posj) for (j, posj) in enumerate(sl.symposlist) if j>i]:
#                    if np.linalg.norm(sl.symposlist[i]-sl.symposlist[j]) < eps: 
#                        found = True
#                        break
#                if not found: newlist.append(sl.symposlist[i])  #cheaper than pop from old list
#            sl.symposlist = newlist
            
#            i = 0
#            while i < len(sl.symposlist):
#                if i % 100 == 0: logging.debug(i)
#                j = i+1
#                while j < len(sl.symposlist):
#                    if np.linalg.norm(sl.symposlist[i]-sl.symposlist[j]) < eps:
#                        sl.symposlist.pop(j)
#                    else:
#                        j += 1
#                i += 1
            
#            i = 0
#            while i < len(sl.symposlist)-1:
#                logging.debug(i)
#                dl = cdist(np.array([list(sl.symposlist[i])]), np.vstack(tuple(sl.symposlist[i+1:])), 'euclidean')
#                j = i+1
#                c = 0
#                while j < len(sl.symposlist):
#                    if dl[0][c] < eps:
#                        del sl.symposlist[j]
#                    else:
#                        j += 1
#                    c += 1
#                i += 1
            
#            i = 0
#            while i < len(sl.symposlist)-1:
#                logging.debug(i)
#                dl = cdist(np.array([list(sl.symposlist[i])]), np.vstack(tuple(sl.symposlist[i+1:])), 'euclidean')
#                zl = list(zip(sl.symposlist[i+1:],dl[0],range(0,len(dl[0]))))
#                zl.sort(key=lambda tu: tu[1])   #sort by distance from point
#                for (j,tu) in enumerate(zl):
#                    if tu[1] > eps:
#                        sl.symposlist = sl.symposlist[:i+1]+list(list(zip(*sorted(zl[j:], key=lambda tu: tu[2])))[0])
#                        break
#                i += 1
#                for d in dl[0]:
#                    if d < eps: sl.symposlist.pop()
#                if min(dl[0]) < eps: newlist.append(posi)