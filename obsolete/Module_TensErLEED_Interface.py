# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

Contains classes and functions for a user interface around the TensorLEED scripts.

##########################################################################################
########   12.09.2019  FILE IS OBSOLETE, REPLACED BY TensErLeedModules PACKAGE    ########
##########################################################################################
"""

#import time
#import sys
import logging
import numpy as np
import re
import fortranformat as ff
import copy

DEFAULT = object()

###############################################
#                CLASSES                      #
###############################################

class Rparams:
    """Stores the parameters found in a PARAMETERS file (default values in __init__)"""
    def __init__(self):
        self.ASAZ = DEFAULT
        self.BLAY = 1           #number of bulk layers
        self.BULKDOUBLING_EPS = 0.001
        self.BULKDOUBLING_ITER = 16
        self.CTRUNC = []        #list of c coordinates separating the layers
        self.ELSPLIT = {}       #if any ELSPLIT is defined, it will be added to the dictionary with the element name as the label and the splitlist as the value
        self.IDEG = 3           # TODO - find out how this is defined / what the default should be! (MIDEG in PARAM)
        self.INPOT_Z = 1.0
        self.INTERLAYER_MAXATTENUATION = 0.0001
        self.PHASESHIFTMIN = 0.05
        self.PHI = 0.0          #from BEAMINCIDENCE
        self.SITEDEF = {}       #labels are the element names, content is dictionaries of format {sitename, list of atom numbers in POSCAR}
        self.SUPERLATTICE = np.array([[1.0,0.0],[0.0,1.0]])
        self.SYMMETRY_EPS = 0.0001
        self.TEOENRANGE = [29, 800, 3]
        self.THETA = 0.0        #from BEAMINCIDENCE
        self.TOUTPUT = []       #defines whether Tensor output is required for a layer. Default = 1.
        
        # script-defined values
        self.ElementAliasDict = {}      #if element names are redefined due to a naming conflict, store the information here
        self.LMAX = 0                   #will be calculated based on PHASESHIFTMIN parameter
        self.phaseshifts = DEFAULT
    
class Sitetype:
    """Site types are identified by (main) element and name, and store vibrational amplitude and occupation"""
    def __init__(self,el,name):
        self.el = el
        self.name = name
        self.label = self.el + '_' + self.name
        self.vibamp = {}    #vibrational amplitude per element
        self.occ = {}       #occupation per element
    
class Atom:
    """To be used for Slab; each atom has an element and a position"""
    def __init__(self,el,pos,oriN,slab):
        self.el = el        #element type
        self.pos = pos      #atom position
        self.oriN = oriN    #original atom number in the POSCAR
        self.slab = slab     #which slab the atom belongs to
        self.layer = DEFAULT     #which layer the atom belongs to. DEFAULT if no layer was assigned yet.
        self.site = DEFAULT   #the site type that the atom is in, supplied by the SITEDEF parameter, assigned by the updateSites routine of Slab
        self.cartpos = DEFAULT    #position in cartesian coordinates, with the highest atom as z = 0, z going down
        self.posInLayer = DEFAULT   #same as cartpos, but from the layer origin

        
class Layer:
    """To be used with Slab; has origin, atoms (a subset of the ones in slab), and a number."""
    def __init__(self,slab,num,isBulk=False):
        self.slab = slab
        self.num = num      # consecutive layer numbering, 0 being highest
        self.isBulk = isBulk    # defined by BLAY in PARAMETERS file
        self.atlist = []    # atoms in this layer
        self.cartori = DEFAULT  # origin: xy from POSCAR with possible displacements, z is highest atom
        self.cartbotz = DEFAULT  # z position of lowest atom
    
    def getLayerPos(self):
        """Gets a cartesian origin coordinate for the layer, using z of the highest atom. x,y are calculated
        from the origin of an a,b unit cell at that height. Also assigns posInLayer for all atoms in this layer."""
        al = self.atlist[:]     #temporary copy
        al.sort(key=lambda atom: atom.pos[2])
        topat = al[-1]
        botat = al[0]
        oripos = np.array([0,0,topat.pos[2]])
        self.cartori = np.dot(self.slab.ucell, oripos)  #this gets x and y correct, but z still in the wrong direction and with origin as POSCAR
        self.cartori[2] = topat.cartpos[2]  #just take the z from the highest atom.
        self.cartbotz = botat.cartpos[2]
#        logging.debug('Layer '+str(self.num)+' origin: '+str(self.cartori))
        for atom in al:
            atom.posInLayer = atom.cartpos - self.cartori
#            logging.debug(atom.el+' '+str(atom.oriN)+' '+str(atom.cartpos)+' '+str(atom.posInLayer))
        
class Slab:
    """Contains unit cell, element information and atom coordinates. After feeding raw data into atpos, call 
    initAtomList to create a list of Atom objects."""
    def __init__(self):
        self.ucell = np.array([])         # unit cell
        self.atpos = []	        # list of atom positions as read from POSCAR (list of arrays)
        self.elements = []      # element labels as read from POSCAR
        self.chemelem = []      # actual chemical elements - including from ELSPLIT!
        self.nelem = 0          # number of different elements - including from ELSPLIT!
        self.nperelem = []      # number of atoms per element
        self.atlist = []        # list of Atom objects; will be sorted by coordinate or by element type
        self.layers = []        # list of Layer objects
        self.sitelist = []      # list of distinct sites as Sitetype elements
        self.uCellTransform = np.array([[1,0],[0,1]])   # stores transformations of the unit cell by multiplication
        self.uCellShift = np.array([0.0,0.0,0.0])       # stores translations of the unit cell
    
    def fullUpdate(self,rparams):
        """readPOSCAR initializes the slab with information from POSCAR; fullUpdate re-initializes the atom list,
        then uses the information from the parameters file to create layers, calculate cartesian coordinates 
        (absolute and per layer), and to update elements and sites."""
        self.initAtomList()
        self.createLayers(rparams)
        self.updateElements(rparams)
        self.updateSites(rparams)
        
    def initAtomList(self):
        """Creates a list of Atom objects based on the data read previously"""
        n = 0
        for nat in self.nperelem:
            n += nat
        self.atlist = []
        elnum = 0
        atcount = 0
        for i in range(0,n):
            self.atlist.append(Atom(self.elements[elnum],self.atpos[i],i+1,self))
            atcount += 1
            if atcount == self.nperelem[elnum]:
                elnum += 1
                atcount = 0
        self.getCartesianCoordinates()
        
    def getCartesianCoordinates(self):
        """Assigns absolute cartesian coordinates to all atoms, with x,y using the unit cell (top plane),
        while z = 0 for the topmost atom and positive going down through the slab."""
        al = self.atlist[:]     #temporary copy
        al.sort(key=lambda atom: atom.pos[2])
        topat = al[-1]
        topcart = np.dot(self.ucell, topat.pos)
        for atom in al:
            atom.cartpos = np.dot(self.ucell, atom.pos)
            atom.cartpos[2] = topcart[2]-atom.cartpos[2]
#            logging.debug(atom.el+' '+str(atom.oriN)+' '+str(atom.pos)+' '+str(atom.cartpos))      
                
    def createLayers(self,rparams):
        """Creates a list of Layer objects based on the BLAY and CTRUNC parameters in rparams.
        If layers were already defined, overwrite."""
        self.layers = []
        tmplist = self.atlist[:]
        self.sortByZ()
        laynum = 0
        b = True if rparams.BLAY > 0 else False
        newlayer = Layer(self,0,b)
        self.layers.append(newlayer)
        for atom in self.atlist:
            if laynum < len(rparams.CTRUNC):        #only check for new layer if we're not in the top layer already
                if atom.pos[2] > rparams.CTRUNC[laynum]:   #if the atom is higher than the next c cutoff, make a new layer
                    laynum += 1
                    b = True if rparams.BLAY > laynum else False
                    newlayer = Layer(self,laynum,b)
                    self.layers.append(newlayer)
                    check = True    #check for empty layer
                    while check:
                        if laynum >= len(rparams.CTRUNC):
                            check = False
                        elif atom.pos[2] <= rparams.CTRUNC[laynum]:
                            check = False
                        else:
                            laynum += 1
                            b = True if rparams.BLAY > laynum else False
                            newlayer = Layer(self,laynum,b)
                            self.layers.append(newlayer)
            atom.layer = newlayer
            newlayer.atlist.append(atom)
        dl = []
        for layer in self.layers:
            if not layer.atlist:
                logging.warning('A layer containing no atoms was found. Layer will be deleted. Check CTRUNC parameter.')
                dl.append(layer)
        for layer in dl:
            if layer.isBulk: self.layers[layer.num+1].isBulk = True
            self.layers.remove(layer)
            del layer
        self.layers.reverse()
        for i, layer in enumerate(self.layers):
            layer.getLayerPos()
            layer.num = i
        self.atlist = tmplist
        
    def updateElements(self,rp):
        """Updates nelem based on the ELSPLIT parameter. Also checks whether any element in ELSPLIT has the same
        name as an element in the POSCAR, and if so, renames the one in the POSCAR."""
        # update nelem
        c = 0
        oldels = self.elements[:]  #copy elements so we don't need to overwrite the list during iteration
        for i, pel in enumerate(oldels):
            if not pel in rp.ELSPLIT:
                c += 1
            else:
                c += len(rp.ELSPLIT[pel])
                # check for overlapping names:
                delkeys = []
                for el in rp.ELSPLIT[pel]:
                    if el in oldels:
                        nn='P'+el
                        logging.warning('element name '+el+' given in ELSPLIT is also an element name in POSCAR. The POSCAR element will be renamed to '+nn+'.')
                        for atom in self.atlist:
                            if atom.el == el: atom.el = nn
                        self.elements[i] = nn
                        rp.ElementAliasDict[el]=nn
                        rp.ELSPLIT[nn] = rp.ELSPLIT[pel]    #copy the ELSPLIT entry to the new element name
                        delkeys.append(pel)         #can't delete old key from ELSPLIT right away since we're iterating over it right now
                for k in delkeys:   #remove old keys
                    del rp.ELSPLIT[k]
        self.nelem = c
        self.chemelem = []
        for el in self.elements:
            if el in rp.ELSPLIT:
                self.chemelem.extend(rp.ELSPLIT[el])
            else:
                self.chemelem.append(el)
        # check if names overlap:
    
    def updateSites(self,rparams):
        """Goes through the atom list and supplies them with appropriate Sitedef objects, based on the SITEDEF parameters 
        from the supplied Rparams."""
        atlist = self.atlist[:]     #copy to not have any permanent changes
        atlist.sort(key=lambda atom: atom.oriN)
        sl = []
        for el, sitedict in rparams.SITEDEF.items():
            aliasfound = False
            for sitename, sitelist in sitedict.items():
                newsite = Sitetype(el, sitename)
                sl.append(newsite)
                for i in sitelist:
                    try:
                        if atlist[i-1].el != el:
                            if atlist[i-1].el == rparams.ElementAliasDict[el]:
                                aliasfound = True
                                atlist[i-1].site = newsite
                            else:
                                logging.warning('SITEDEF tries to assign atom number '+str(i)+' as '+el+', but POSCAR has it as '+atlist[i-1].el+'. Atom will be skipped and left as default site type!')
                        else:
                            atlist[i-1].site = newsite
                    except IndexError:
                        logging.error('SITEDEF: atom number out of bounds.')
            if aliasfound: 
                logging.debug('SITEDEF assignments to element '+el+' were redirected to element '+rparams.ElementAliasDict[el])
                newsite.el=rparams.ElementAliasDict[el]
                newsite.label = newsite.el + '_' + newsite.name
        for el in self.elements:
            newsite = Sitetype(el, 'def')
            found = False
            for at in atlist:
                if at.el == el and at.site == DEFAULT:
                    at.site = newsite
                    found = True
            if found: sl.append(newsite)
        self.sitelist = sl
      
    def sortByZ(self, topToBot=False):
        """Sorts atlist by z coordinate"""
        self.atlist.sort(key=lambda atom: atom.pos[2])
        if topToBot:
            self.atlist.reverse()
        
    def sortByEl(self):
        """Sorts atlist by elements, preserving the element order from the original POSCAR"""
        # unfortunately, simply calling the sort function by element does not preserve the element order from the POSCAR
        esortlist = sorted(self.atlist, key=lambda atom: atom.el)   #create a list sorted by element
        lastel = ''
        tmpElList = []
        isoLists = []
        # generate sub-lists isolated by elements; probably could be done cleaner...
        for at in esortlist:
            if at.el != lastel:
                tmpElList.append(at.el)
                isoLists.append([])
                lastel = at.el
            isoLists[-1].append(at)
        sortedlist = []
        # going through the elements in the order they appear in POSCAR, find the corresponding index in tmpElList and append the atoms of that type to sorted list
        for el in self.elements:
            i = tmpElList.index(el)
            for at in isoLists[i]:
                sortedlist.append(at)
        self.atlist = sortedlist
    
    def sortOriginal(self):
        """Sorts atlist by original atom order from POSCAR"""
        self.atlist.sort(key=lambda atom: atom.oriN)
        
    def projectCToZ(self):
        """makes the c vector of the unit cell perpendicular to the surface, changing 
        all atom coordinates to fit the new base"""
        if self.ucell[0,2] != 0.0 or self.ucell[1,2] != 0.0:
            self.getCartesianCoordinates()
            ab = self.ucell[0:2,0:2]
            alpha = angle(ab[:,0],ab[:,1])
            for at in self.atlist:
                xyCart = at.cartpos[0:2]
                xyCart[0] = xyCart[0] % np.linalg.norm(ab[:,0])
                xyCart[1] = xyCart[1] % np.linalg.norm(ab[:,1])
                for i in range(0,2):
                    proj = (np.dot(xyCart,ab[:,i]) / (np.linalg.norm(ab[:,i])**2)) * ab[:,i]
                    at.pos[i] = ((np.linalg.norm(proj)-(np.linalg.norm(xyCart-proj)/np.tan(alpha))) % np.linalg.norm(ab[:,i])) / np.linalg.norm(ab[:,i])
            self.ucell[0,2] = 0
            self.ucell[1,2] = 0
            self.getCartesianCoordinates()

    def findSymmetry(self, rp):
        """Reduces the unit cell if necessary and finds the point group of the slab. Returns...""" #TODO: complete later
        celltype = "ERROR - not recognized"     #storing this as string is computationally inefficient, but readability seems more important than speed here
        eps = rp.SYMMETRY_EPS
        # reduce surface unit cell
        abst = np.transpose(self.ucell[0:2,0:2])     #surface unit cell, transposed
        usurf = np.array([[1,0],[0,1]])             #track unit cell changes
        dp = np.dot(abst[0],abst[1])
        if abs(dp) < eps:
            if abs(np.linalg.norm(abst[0])-np.linalg.norm(abst[1])) < eps:
                celltype = "square"
            else:
                celltype = "rectangular"
        else:
            #rhombic, hexagonal or oblique
            if dp > 0:      #transform acute to obtuse
                abst = np.dot(np.array([[-1,0],[0,1]]), abst)
                usurf = np.dot(np.array([[-1,0],[0,1]]), usurf)
            (abst, u) = reduceUnitCell(abst, eps)
            usurf = np.dot(u,usurf)
        # reduce bulk unit cell
        abbt = np.dot(np.linalg.inv(rp.SUPERLATTICE),abst)      #bulk ab unit cell, transposed
        ubulk = np.array([[1,0],[0,1]])                       #track unit cell changes
        dp = np.dot(abbt[0],abbt[1])
        if abs(dp) > eps:    #rhombic, hexagonal or oblique
            if dp > 0:      #transform acute to obtuse
                abbt = np.dot(np.array([[-1,0],[0,1]]), abbt)
                ubulk = np.dot(np.array([[-1,0],[0,1]]), ubulk)
            (abbt, u) = reduceUnitCell(abbt, eps)
            ubulk = np.dot(u,ubulk)
        # MODIFY SUPERLATTICE PARAMETER         #TODO: This should probably also be written to the parameters file?
        rp.SUPERLATTICE = np.dot(usurf, np.dot(rp.SUPERLATTICE, np.linalg.inv(ubulk)))
        # MODIFY UNIT CELL
        utr = np.array([[0,0,0],[0,0,0],[0,0,1]])
        utr[0:2,0:2] = usurf
        self.uCellTransform = utr
        if not np.array_equal(utr, np.array([[1,0,0],[0,1,0],[0,0,1]])):
            if np.array_equal(utr, np.array([[-1,0,0],[0,1,0],[0,0,1]])):
                logging.info("The POSCAR unit cell was changed from an acute to an obtuse form.")
            else:
                logging.warning("The POSCAR unit cell was not in its highest symmetry form. The unit cell will be modified with the transformation matrix: \n"+str(utr))
            self.ucell = np.dot(self.ucell, np.transpose(utr))  #same as np.transpose(np.dot(utr,np.transpose(self.ucell)))
            for at in self.atlist:
                at.pos = np.dot(np.linalg.inv(np.transpose(utr)),at.pos)    #inversion and transposition commute
            self.getCartesianCoordinates()
        # check cell type again
        abst = np.transpose(self.ucell[0:2,0:2])
        dp = np.dot(abst[0],abst[1])
        if abs(dp) < eps:  #square or rectangular
            if abs(np.linalg.norm(abst[0])-np.linalg.norm(abst[1])) < eps:
                celltype = "square"
            else:
                celltype = "rectangular"
        elif np.linalg.norm(abst[0]) - np.linalg.norm(abst[1]) >= eps:
            celltype = "oblique"
        elif angle(abst[0],abst[1]) - (2*np.pi/3) < eps:
            celltype = "hexagonal"
        else:
            celltype = "rhombic"
        logging.debug("Unit cell type: "+celltype)
        # FIND HIGHEST SYMMETRY ORIGIN
        ts = copy.deepcopy(self)
        ts.projectCToZ()
        ts.sortByZ(topToBot=True)
        
        return
        
    def revertUnitCell(self):
        """If the unit cell in a and b was transformed earlier, restore the original form and coordinates."""
        self.ucell = np.dot(self.ucell, np.linalg.inv(np.transpose(self.uCellTransform)))     #same as self.ucell = np.transpose(np.dot(np.linalg.inv(self.uCellTransform),np.transpose(self.ucell)))
        for at in self.atlist:
            at.pos = np.dot(np.transpose(self.uCellTransform),at.pos)
        self.getCartesianCoordinates()
        
###############################################
#                FUNCTIONS                    #
###############################################

def parseMathSqrt(s):
    try:
        f = float(s)
    except ValueError:
        f = 1.0
        if '*' in s:
            sl = s.split('*')
        else:
            sl = [s]
        for el in sl:
            try:
                f *= float(el)
            except ValueError:
                if not 'sqrt' in el:
                    logging.error('Could not interpret '+el+' as float value')
                else:
                    p = re.compile(r'\s*sqrt\s*\(\s*(?P<value>[\d.]*)\s*\)')
                    m = p.match(el)
                    if m:
                        try:
                            f *= np.sqrt(float(m.group('value')))
                        except ValueError:
                            logging.error('Could not interpret '+m.group('value')+' as float value')
                    else:
                        logging.error('Could not interpret '+el+' as float value')
    return f

def angle(v1, v2, acute=True):
    angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    if acute == True:
        return angle
    else:
        return 2 * np.pi - angle
    
def readWoodsNotation(s, ucell):
    """Takes a string that should contain the transformation from the bulk to the surface
    unit cell in Wood notation, as well as a bulk unit cell (from which only the surface vectors 
    are read). Returns a 2x2 transformation matrix."""
    p = re.compile(r'\s*(?P<type>[PCpc]*)\s*\(\s*(?P<g1>.+)\s*[xX]\s*(?P<g2>.+)\s*\)\s*[rR]*\s*(?P<alpha>[\d.]*)')
    # this regular expression matches if (any amount of whitespace at any point is ignored):
    #   - optional: first character is p or c (or P or C) (-> type)
    #   - then there is a '('
    #   - then whatever (at least one character), interpret later (-> g1)
    #   - then 'x' or 'X'
    #   - then whatever (at least one character), interpret later (-> g2)
    #   - then ')'
    #   - then (optional) 'r' or 'R'
    #   - then (optional) an integer or float number (-> alpha)
    m = p.match(s)
    if not m:
        logging.critical('Could not read woods notation input '+s)
        return
    if not m.group('type'):
        t = 'p'
    else:
        t = m.group('type').lower()
    if not m.group('alpha'):
        alpha = 0.0
    else:
        try:
            alpha = float(m.group('alpha'))
        except:
            logging.error('Could not read Woods notation angle: '+m.group('alpha')+', setting angle to zero')
            alpha = 0.0
    alpha *= np.pi/180
    g1 = parseMathSqrt(m.group('g1'))
    g2 = parseMathSqrt(m.group('g2'))
    #logging.debug(t+', '+str(g1)+', '+str(g2)+', '+str(alpha*180/np.pi))
    # get surface unit cell vectors from bulk unit cell (has surface periodicity!!):
    r = [ucell[:2,0],ucell[:2,1]]
    #q = np.linalg.norm(r[1])/np.linalg.norm(r[0])  #this would be to get from bulk vectors to surface, we have to reverse
    q = 1/(np.linalg.norm(r[1])/np.linalg.norm(r[0]))
    omega = angle(r[0],r[1])    #this is always constant in Wood notation, no need to reverse.
    #logging.debug(str(alpha)+' '+str(omega)+' '+str(g1)+' '+str(g2)+' '+str(q))
    if t == 'p':    #matrices from: Klaus Hermann; Crystallography and Surface Structure (Second Edition, Wiley)
        mat = (1/np.sin(omega))*np.array([[g1*np.sin(omega-alpha), g1*(1/q)*np.sin(alpha)], [-g2*q*np.sin(alpha), g2*np.sin(omega+alpha)]])
    else:
        mat = (1/(2*np.sin(omega)))*np.array([[g1*np.sin(omega-alpha)-g2*q*np.sin(alpha), g1*(1/q)*np.sin(alpha)+g2*np.sin(omega+alpha)], [-g1*np.sin(omega-alpha)-g2*q*np.sin(alpha),-g1*(1/q)*np.sin(alpha)+g2*np.sin(omega+alpha)]])
    return mat
     
def linelist(line):
    """splits a line at whitespace, deletes empty elements and line breaks, then returns elements as a list"""
    llist1 = line.split()
    llist = []
    for part in llist1:
        if part != "":	# get rid of empty elements
            if part[-1][-1:] != '\n':   # get rid of line breaks
                llist.append(part)
            else:
                if part != '\n':	# don't append if the only element was the line break
                    llist.append(part[:-1])
    return llist

def readPOSCAR(filename='POSCAR'):
    """Reads a POSCAR and returns a Slab object with the information"""

    # open input file
    try:
        rf = open(filename, 'r')
    except FileNotFoundError:
        logging.error("POSCAR not found.")
        raise

    #read POSCAR:
    linenum = 1		# iterates the current line being read
    read = True		# set false to stop reading
    sl = Slab()
    for line in rf:
        if linenum == 1:
            pass
        elif linenum == 2:
            scaling = float(linelist(line)[0])
        elif linenum <= 5:
            if linenum == 3: ucellList = []
            llist = [float(i) for i in linelist(line)]
            ucellList.append(llist)
            if linenum == 5: 
                sl.ucell = scaling * np.transpose(np.array(ucellList))    #I'm not 100% sure if this needs transpose or not. check later!
                if sl.ucell[2,0] != 0.0 or sl.ucell[2,1] != 0.0:
                    logging.error('ERROR: Unit cell a and b vectors must not have an out-of-surface (Z) component!')
                    raise
        elif linenum == 6:
            sl.elements = linelist(line)			# element labels
            sl.nelem = len(sl.elements)					# number of different elements
        elif linenum == 7:
            sl.nperelem = [int(i) for i in linelist(line)]			# number of atoms per element
            if len(sl.nperelem) != sl.nelem: logging.warning('\nPOSSIBLE PROBLEM: lenght of element list does not match length of atoms-per-element list\n')
        elif linenum == 8:
            pass		# skip
        elif linenum == 9:
            # this line might already contain coordinates, or not, depending on whether the "Selective dynamics" line is there
            llist = linelist(line)
            try:
                pos = np.array([float(llist[0]), float(llist[1]), float(llist[2])])
                for i in range(0,2):    # in a and b, make sure the values are between 0 and 1
                    pos[i] = pos[i] % 1.0
                sl.atpos.append(pos)  # if the line does NOT start with a list of 3 floats, this will raise an exception
            except:
                logging.debug("POSCAR contains 'Selective dynamics' line, skipping line 9")
                pass		# exception was raised because of the 'Selective dynamics' line; this is fine, move on.
        elif read:
            llist = linelist(line)
            if len(llist) == 0:
                read = False
                logging.debug("POSCAR: Empty line found; stopping position readout")
            else:
                pos = np.array([float(llist[0]), float(llist[1]), float(llist[2])])
                for i in range(0,2):    # in a and b, make sure the values are between 0 and 1
                    pos[i] = pos[i] % 1.0
                sl.atpos.append(pos)
        linenum += 1
    rf.close()
    sl.initAtomList()
    return(sl)

def writeCONTCAR(sl, filename='CONTCAR', reorder=False, comments=False):
    """Takes a Slab object and writes a CONTCAR based on the atlist. If a POSCAR is in the folder, 
    it will take header data from it, otherwise is will fill in a generic header. Defining a filename 
    other than CONTCAR is optional. If the reorder tag is set as True, atoms will be sorted by element 
    and by z coordinate; otherwise, the original order will be used."""
    if reorder:
        sl.sortByZ()
        sl.sortByEl()
#    else:
#        sl.sortOriginal()
    output = '' 	# store output here temporarily before writing it
    # see if a POSCAR is there
    epos = True
    try:
        rf = open('POSCAR', 'r')
    except:
        epos = False
    #header
    if epos:
        #copy the first two lines:
        output += rf.readline()
        output += rf.readline()
        rf.close()
    else:
        #fill dummy header:
        output += 'unknown\n'
        output += '   1.0\n'
    #unit cell
    m = np.transpose(sl.ucell)
    for vec in m:
        for val in vec:
            s = '%.16f'%val
            output += '{:>22}'.format(s)
        output += '\n'
    #atom types and numbers
    for el in sl.elements:
        output += '{:>5}'.format(el)
    output += '\n'
    for n in sl.nperelem:
        output += '{:>5}'.format(n)
    output += '\n'
    #header line 'Direct'
    ol = 'Direct'
    if comments:
        ol = ol.ljust(60)
        ol += 'N'.rjust(5)+'NperEl'.rjust(9)+'SiteLabel'.rjust(12)+'Layer'.rjust(7)+'Linking'.rjust(9)
    output += ol+'\n'
    NperElCount = 0
    lastEl = sl.atlist[0].el
    for i, at in enumerate(sl.atlist):
        if lastEl != at.el:
            lastEl = at.el
            NperElCount = 1
        else:
            NperElCount += 1
        ol = ''
        for coord in at.pos:
            s = '%.16f'%coord
            ol += '{:>20}'.format(s)
        if comments:
            ol += str(i+1).rjust(5)                         #N
            ol += (at.el+str(NperElCount)).rjust(9)         #NperEl
            ol += at.site.label.rjust(12)                   #SiteLabel
            ol += str(at.layer.num+1).rjust(7)
        output += ol+'\n'
    #write output
    try:
        wf = open(filename, 'w')
        wf.write(output)
        wf.close()
    except:
        logging.error("Failed to write "+filename)
        raise
    logging.info("Wrote to "+filename+" successfully")
    return

def readToExc(llist):
    """For reading PARAMETERS files; takes a list, returns elements until the first one that starts with an 
    exclamation mark."""
    read = True
    newlist = []
    for s in llist:
        if read:
            if s[0] == '!':
                read = False
            else:
                newlist.append(s)
    return newlist        

def splitSublists(llist, sep):
    """Takes a list and a separator, splits strings in the list by the separator, returns results as list of lists"""
    newlist = []
    sublist = []
    for el in llist:
        if not sep in el:
            sublist.append(el)
        else:
            tl = el.split(sep)
            if tl[0]: sublist.append(tl[0])
            newlist.append(sublist)
            sublist = []
            if len(tl) > 1:
                for i in range(1, len(tl)):
                    s = tl[i]
                    if s: sublist.append(s)
                    if not i==len(tl)-1:
                        newlist.append(sublist)
                        sublist = []
    newlist.append(sublist)
    return(newlist)
    
def splitMaxRight(s, sep):
    """Same as s.split(sep, maxsplit=1), but splitting at the first instance from the right."""
    sr = s[::-1]
    l = sr.split(sep, maxsplit=1)
    l.reverse()
    nl = []
    for ns in l: 
        nl.append(ns[::-1])
    return nl

def recombineListElements(llist, com):
    """Takes a list, checks in each element whether the first/last characters
    are the given combination character, and if so, combines list elements with
    the list element before/after."""
    i = 0
    newlist = llist[:]
    while i < len(newlist)-1:
        if newlist[i][-1] == com or newlist[i+1][0] == com:
            newlist[i] += newlist.pop(i+1)
        else:
            i += 1
    return newlist

def readPARAMETERS(filename='PARAMETERS', slab=DEFAULT):
    """Reads a PARAMETERS and returns an Rparams object with the information"""
    # open input file
    try:
        rf = open(filename, 'r')
    except FileNotFoundError:
        logging.error("PARAMETERS not found.")
        raise
    #read PARAMETERS:
    rpars = Rparams()
    for line in rf:
        readvalue = False
        if "=" in line:                         #ignore all lines that don't have an "=" sign at all
            param = line.split('=')[0]          #parameter is defined left of "="
            if param: 
                readvalue = True
                plist = param.split()  #if param is not empty, get rid of spaces and check the leftmost entry.
                if plist: param = plist[0]
                if param[0] == '!': readvalue = False
            if readvalue:
                try:
                    value = line.split('=')[1]
                    llist = linelist(value)            #read the stuff to the right of "="
                except IndexError:
                    logging.warning('PARAMETERS file: ' + param + ' appears to have no value')
                    readvalue = False
                if not llist:
                    logging.warning('PARAMETERS file: ' + param + ' appears to have no value')
                    readvalue = False
            if readvalue:
                if param == 'ASAZ':
                    try:
                        rpars.ASAZ = float(llist[0])
                    except ValueError:
                        logging.warning('PARAMETERS file: ASAZ: Could not convert value to float. Input will be ignored.')
                elif param == 'BEAMINCIDENCE':
                    llist = readToExc(llist)
                    if ',' in value:
                        try:
                            sublists = splitSublists(llist, ',')
                            for sl in sublists:
                                if sl[0].lower() == 'theta':
                                    theta = float(sl[1])
                                    if 0 <= theta and theta <= 90:
                                        rpars.THETA = float(theta)
                                    else:
                                        logging.warning('PARAMETERS file: BEAMINCIDENCE: Unexpected value for theta (should be between 0 and 90). Input will be ignored.')
                                elif sl[0].lower() == 'phi':
                                    rpars.PHI = float(sl[1])%360
                                else:
                                    logging.warning('PARAMETERS file: BEAMINCIDENCE: Unknown flag found. Input will be ignored.')
                        except ValueError:
                            logging.warning('PARAMETERS file: BEAMINCIDENCE: Could not convert value to float. Input will be ignored.')
                    else:
                        try:
                            theta = float(llist[0])
                            if 0 <= theta and theta <= 90:
                                rpars.THETA = float(theta)
                            else:
                                logging.warning('PARAMETERS file: BEAMINCIDENCE: Unexpected value for theta (should be between 0 and 90). Input will be ignored.')
                            rpars.PHI = float(llist[1])
                        except ValueError:
                            logging.warning('PARAMETERS file: BEAMINCIDENCE: Could not convert value to float. Input will be ignored.')
                elif param == 'BLAY':
                    rpars.BLAY = int(llist[0])
                    if rpars.BLAY != 1 and rpars.BLAY != 2:
                        logging.warning('PARAMETERS file: BLAY was set to '+rpars.BLAY+', 1 or 2 expected. Value was set to 1 (default).')
                        rpars.BLAY = 1
                elif param == 'BULKDOUBLING_EPS':
                    try:
                        rpars.BULKDOUBLING_EPS = float(llist[0])
                    except ValueError:
                        logging.warning('PARAMETERS file: BULKDOUBLING_EPS: Could not convert value to float. Input will be ignored.')
                    if rpars.BULKDOUBLING_EPS < 0.001:
                        rpars.BULKDOUBLING_EPS = 0.001
                        logging.warning('PARAMETERS file: BULKDOUBLING_EPS cannot be smaller than 0.001 due to fortran reading it as an F6.3; value was changed to 0.001')
                elif param == 'BULKDOUBLING_ITER':
                    try:
                        i = int(llist[0])
                    except ValueError:
                        logging.warning('PARAMETERS file: BULKDOUBLING_ITER: Could not convert value to integer. Input will be ignored.')
                    if i > 0:
                        rpars.BULKDOUBLING_ITER = i
                    else:
                        logging.warning('PARAMETERS file: BULKDOUBLING_ITER: Unexpected input (0 or negative). Input will be ignored.')
                elif param == 'CTRUNC':
                    for s in readToExc(llist):
                        rpars.CTRUNC.append(float(s))
                    rpars.CTRUNC.sort()
                elif param == 'ELSPLIT':
                    rpars.ELSPLIT[plist[1]]=readToExc(llist)
                elif param == 'INPOT_Z':
                    try:
                        rpars.INPOT_Z = float(llist[0])
                    except ValueError:
                        logging.warning('PARAMETERS file: INPOT_Z: Could not convert value to float. Input will be ignored.')
                elif param == 'INTERLAYER_MAXATTENUATION':
                    try:
                        f = float(llist[0])
                        if f >= 0.0001 and f < 1:
                            rpars.INTERLAYER_MAXATTENUATION = f
                        else:
                            logging.warning('PARAMETERS file: Unexpected input for INTERLAYER_MAXATTENUATION. Input will be ignored.')
                    except ValueError:
                        logging.warning('PARAMETERS file: INTERLAYER_MAXATTENUATION: Could not convert value to float. Input will be ignored.')
                elif param == 'PHASESHIFTMIN':
                    try:
                        f = float(llist[0])
                    except ValueError:
                        s = llist[0].lower()[0]
                        if s == 'r':    #rough
                            f = 0.1
                        elif s == 'n':  #normal
                            f = 0.05
                        elif s == 'f':  #fine
                            f = 0.01
                        elif s == 'e':  #extrafine
                            f = 0.001
                        else:
                            logging.warning('PARAMETERS file: PHASESHIFTMIN: Could not convert value to float. Input will be ignored.')
                    if f > 0 and f < 1:
                        rpars.PHASESHIFTMIN = f
                    else: 
                        logging.warning('PARAMETERS file: PHASESHIFTMIN: Unexpected value (should be between 0 and 1). Input will be ignored.')
                elif param == 'SITEDEF':
                    newdict = {} 
                    sublists = splitSublists(readToExc(llist), ',')
                    for sl in sublists:
                        atnums = []
                        for i in range(1, len(sl)):
                            try:
                                atnums.append(int(sl[i]))
                            except ValueError:
                                if "-" in sl[i]:
                                    spl = sl[i].split("-")
                                    atnums.extend(list(range(int(spl[0]),int(spl[1])+1)))
                                elif ":" in sl[i]:
                                    spl = sl[i].split(":")
                                    atnums.extend(list(range(int(spl[0]),int(spl[1])+1)))
                                elif "top(" in sl[i]:
                                    if slab == DEFAULT:
                                        logging.warning('PARAMETERS file: SITEDEF parameter contains a top() function, but no slab was passed. The atoms will be assigned the default site type instead.')
                                    else:
                                        n = int(sl[i].split('(')[1].split(')')[0])
                                        csatlist = slab.atlist[:]
                                        csatlist.sort(key=lambda atom: atom.pos[2])
                                        while n > 0:
                                            at = csatlist.pop()
                                            if at.el == plist[1]:
                                                atnums.append(at.oriN)
                                                n -= 1
                                else:
                                    logging.error('PARAMETERS file: Problem with SITEDEF input format')
                                    raise
                        newdict[sl[0]] = atnums
                    rpars.SITEDEF[plist[1]]=newdict
                elif param == 'SUPERLATTICE':
                    if not 'M' in plist:
                        if slab == DEFAULT:
                            logging.critical('PARAMETERS file: SUPERLATTICE parameter appears to be in Wood notation, but no slab was passed; cannot calculate bulk unit cell!')
                        else:
                            rpars.SUPERLATTICE = readWoodsNotation(value, slab.ucell)  #don't need to cut off after exclamation marks since those will anyway not match the regex
                    else:
                        sublists = splitSublists(readToExc(llist), ',')
                        if not len(sublists) == 2:
                            logging.critical('PARAMETERS file: error reading SUPERLATTICE matrix: number of lines is not equal 2.')
                        else:
                            write=True
                            nl = []
                            for sl in sublists:
                                if len(sl) == 2:
                                    try:
                                        nl.append([float(s) for s in sl])
                                    except ValueError:
                                        logging.critical('PARAMETERS file: error reading SUPERLATTICE matrix: could not convert '+str(sl)+' to floats.')
                                        write = False
                                else:
                                    logging.critical('PARAMETERS file: error reading SUPERLATTICE matrix: number of columns is not equal 2.')
                                    write = False
                            if write:
                                rpars.SUPERLATTICE = np.array(nl)
                elif param == 'TEOENRANGE':
                    try:
                        fl = [float(s) for s in readToExc(llist)]
                    except ValueError:
                        logging.warning('PARAMETERS file: Failed to convert TEOENRANGE input to floats. Input will be ignored')
                    if fl:
                        if len(fl) == 3 and fl[0]>0 and fl[1]>fl[0] and fl[2]>0 and fl[1]-fl[0]>fl[2]:
                            if (fl[1]-fl[0])%fl[2] != 0: 
                                fl[0] -= fl[2]-(fl[1]-fl[0])%fl[2]     #if the max is not hit by the steps exactly, correct max up to make it so
                                logging.debug('TEOENRANGE parameter: (Eto-Efrom)%Estep != 0, Efrom was corrected to '+str(fl[0]))
                            rpars.TEOENRANGE = fl
                        else:
                            logging.warning('PARAMETERS file: Unexpected input for TEOENRANGE. Input will be ignored.')
                    else:
                        logging.warning('PARAMETERS file: TEOENRANGE appears to have no value.')
                elif param == 'TOUTPUT':
                    nl = recombineListElements(llist, '*')
                    for s in nl:
                        try:
                            v = int(s)
                            if v != 0 and v != 1:
                                logging.warning('PARAMETERS file: Problem with TOUTPUT input format: Found value '+str(v)+', expected 0 or 1. Value will be ignored.')
                            else:
                                rpars.TOUTPUT.append(v)
                        except ValueError:
                            if '*' in s:
                                sl = s.split('*')
                                try:
                                    r = int(sl[0])
                                    v = int(sl[1])
                                except ValueError:
                                    logging.warning('PARAMETERS file: Problem with TOUTPUT input format: could not read value '+s+', value will be ignored.')
                                if v != 0 and v != 1:
                                    logging.warning('PARAMETERS file: Problem with TOUTPUT input format: Found value '+str(v)+', expected 0 or 1. Value will be ignored.')
                                else:
                                    for i in range(0,r):
                                        rpars.TOUTPUT.append(v)
                            else:
                                logging.warning('PARAMETERS file: Problem with TOUTPUT input format: could not read value '+s+', value will be ignored.')
                #elif some other parameter... would go here
    rf.close()
    logging.info("PARAMETERS file was read successfully")
    # INITIALIZE DERIVATIVE PARAMETERS:
    # LMAX
    phaseshifts = []
    try:
        phaseshifts = readPHASESHIFTS()
    except:
        logging.error('Error while reading _PHASESHIFTS file: ', exc_info=True)
    if phaseshifts:
        lmax = 1
        for el in phaseshifts[-1][1]:   #highest energy - go through elements
            for i, val in enumerate(el):
                if val > rpars.PHASESHIFTMIN and (i+1) > lmax:
                    lmax = i+1
        if lmax < 10: logging.warning('The LMAX found based on the PHASESHIFTMIN parameter is very small (LMAX='+str(lmax)+'). Calculation will proceed, but this might cause problems.')
        if lmax > 15:
            lmax = 15
            logging.info('The LMAX found based on the PHASESHIFTMIN parameter is greater than 15, which is currently not supported. LMAX was set to 15.')
        rpars.LMAX = lmax
        rpars.phaseshifts = phaseshifts
    #finished - return
    return(rpars)
    
def readVIBOCCIN(rp, slab, filename='VIBOCCIN'):
    """Reads VIBOCCIN and adds the information to all sites in the slab."""
    # open input file
    try:
        rf = open(filename, 'r')
    except FileNotFoundError:
        logging.error("VIBOCCIN not found.")
        raise
    #read VIBOCCIN:
    mode = 0        #0: not reading; 1: reading vibrational amplitudes, 2: reading occupations
    regex = False   #read regular expressions as-is or not
    for line in rf:
        if '=' in line:
            if line[0] == '=':
                llist = line[1:].split()
                if llist[0][0].lower() == 'v':
                    mode = 1
                elif llist[0][0].lower() == 'o':
                    mode = 2
                elif llist[0][0].lower() == 'r':
                    regex = True
                    if len(llist) >= 2:
                        if llist[1].lower() == 'off':
                            regex = False
                else:
                    logging.warning("VIBOCCIN: Found line starting with '=', but didn't recognize vibrational amplitude or occupation")
            elif mode == 0:
                pass
            else:
                param = line.split('=')[0]
                if param:
                    readvalue = True
                    plist = param.split()  #if param is not empty, get rid of spaces and check the leftmost entry.
                    if plist: param = plist[0]
                    if param[0] == '!': readvalue = False
                if readvalue:
                    try:
                        llist = linelist(line.split('=')[1])            #read the stuff to the right of "="
                    except IndexError:
                        logging.warning('VIBOCCIN file: ' + param + ' appears to have no value')
                        readvalue = False
                    if not llist:
                        logging.warning('VIBOCCIN file: ' + param + ' appears to have no value')
                        readvalue = False
                if readvalue:
                    # first get values on the right
                    sublists = splitSublists(readToExc(llist), ',')
                    # read parameter on the left
                    invDict = {v: k for k, v in rp.ElementAliasDict.items()}    #invert mapping of ElementAliasDict
                    targetsites = []
                    for site in slab.sitelist:
                        if site.el in invDict:     #we have element name substitution, replace to new name
                            p2 = param.replace(invDict[site.el], site.el, 1)
                        else:
                            p2 = param
                        if regex:
                            prep = p2    #if regular expressions are enabled, take the parameter input at face value
                        else:
                            prep = re.escape(p2)    #this double-slashes every non-literal character to escape the regular expressions
                            prep = prep.replace('\\*','.*')  #if regular expressions are not enabled, we want to still interpret * as "any number of any characters"
                        m = re.match(prep, site.label)
                        if m:
                            if m.end(0) == len(site.label):     #if the length of the matched text == the site label, it's a real match
                                targetsites.append(site)
                    if len(targetsites) == 0:
                        logging.warning('VIBOCCIN file: No sites matching '+param+' found, line will be skipped.')
                    else:
                        for site in targetsites:
                            if mode == 1:   #decide which dictionary to write to
                                td = site.vibamp
                            else:
                                td = site.occ
                            for sl in sublists:
                                if len(sl) == 1:    #if there's only one value, it should be a float for the main site element
                                    try:
                                        if site.el in rp.ELSPLIT:
                                            for subel in rp.ELSPLIT[site.el]:
                                                td[subel] = float(sl[0])
                                        else:
                                            td[site.el] = float(sl[0])
                                    except:
                                        logging.error('VIBOCCIN file: Error reading value '+sl[0]+' at parameter '+param)
                                        raise
                                else:
                                    try:
                                        td[sl[0]] = float(sl[1])
                                    except:
                                        logging.error('VIBOCCIN file: Error reading value '+sl[1]+' at parameter '+param)
                                        raise
    logging.info("VIBOCCIN file was read successfully")
    # now fill up default values & do consistency checks:
    for site in slab.sitelist:
        if not site.el in site.vibamp:  #check if the vibrational amplitude is defined for the main element(s)
            if not site.el in rp.ELSPLIT:
                logging.critical('VIBOCCIN file: No '+site.el+' vibrational amplitude defined for site '+site.el+'_'+site.name)
            else:
                v = -1
                for subel in rp.ELSPLIT[site.el]:
                    if subel in site.vibamp:
                        v = site.vibamp[subel]
                if v == -1:
                    logging.critical('VIBOCCIN file: No vibrational amplitude defined for any of the main elements in site '+site.el+'_'+site.name)
                else:
                    # vibrational amplitudes were defined for some of the main elements but not all. Fill the other values
                    for subel in rp.ELSPLIT[site.el]:
                        if not subel in site.vibamp:
                            logging.warning('VIBOCCIN file: No '+subel+' vibrational amplitude defined for site '+site.el+'_'+site.name+'. Using vibrational amplitude from other element in site.')
                            site.vibamp[subel] = v
        if site.el in rp.ELSPLIT:
            mainva = site.vibamp[rp.ELSPLIT[site.el][0]]    #just use first element in ELSPLIT list
        else:
            mainva = site.vibamp[site.el]
        for el in slab.chemelem:
            if not el in site.vibamp:
                site.vibamp[el] = mainva    #for other elements, fill vibrational elements with that of the main element (probably not present)
        # check and fill occupations:
        o = 0.0
        for el in slab.chemelem:
            if el in site.occ:
                o += site.occ[el]
        if not site.el in site.occ:
            if site.el in rp.ELSPLIT:
                found = False
                for el in rp.ELSPLIT[site.el]:
                    if el in site.occ: found = True
                if not found:
                    for el in rp.ELSPLIT[site.el]:
                        site.occ[el] = (1-o)/len(rp.ELSPLIT[site.el])
            else:
                site.occ[site.el] = 1.0
        for el in slab.chemelem:
            if not el in site.occ:
                site.occ[el] = 0.0
        if o > 1.0:
            logging.warning('VIBOCCIN file: Site '+site.el+'_'+site.name+' has a total occupation greater than one ('+str(o)+').')
        # finally, check whether we have any vibrational amplitudes or occupations assigned to non-existant elements, and if so, drop them and warn the user about it:
        dl = []
        for el in site.vibamp:
            if not el in slab.chemelem:
                logging.warning('VIBOCCIN file: Site '+site.el+'_'+site.name+' has a vibrational amplitude defined for an unknown element, which will be dropped ('+el+').')
                dl.append[el]
        for el in dl:
            site.vibamp.pop(el, None)
        dl = []
        for el in site.occ:
            if not el in slab.chemelem:
                logging.warning('VIBOCCIN file: Site '+site.el+'_'+site.name+' has an occupation defined for an unknown element, which will be dropped ('+el+').')
                logging.warning('VIBOCCIN file: Site '+site.el+'_'+site.name+' has an occupation defined for an unknown element, which will be dropped ('+el+').')
                dl.append[el]
        for el in dl:
            site.occ.pop(el, None)
    logging.info("VIBOCCIN value consistency check finished.")                     
            
    
def readIVBEAMS(filename='IVBEAMS'):
    """Reads an IVBEAMS file and returns a list of beams"""
    # open input file
    try:
        rf = open(filename, 'r')
    except FileNotFoundError:
        logging.error("IVBEAMS not found.")
        raise
    #read IVBEAMS:
    linenum = 1		# iterates the current line being read
    beams = []
    for line in rf:
        if linenum == 1:
            pass
        else:
            llist = linelist(line)
            if len(llist) == 1:
                logging.warning('An line with only one element was found in IVBEAMS and will be skipped: '+line)
            elif len(llist) >= 2:
                f = [0,0]
                for i in range(0,2):
                    try:
                        f[i] = float(llist[i])
                    except ValueError:
                        if '/' in llist[i]:
                            f[i] = float(llist[i].split('/')[0])/float(llist[i].split('/')[1])
                        else:
                            logging.error('Error reading IVBEAMS line: '+line)
                            raise
                if not [f[0],f[1]] in beams:
                    beams.append([f[0],f[1]])
        linenum += 1
    rf.close()
    logging.info("IVBEAMS file was read successfully")
    return(beams)
    
def writeAUXBEAMS(beams=[], beamsfile='IVBEAMS', readfile='_BEAMLIST', writefile='AUXBEAMS', write=True):
    """"Reads from a _BEAMLIST file (full list of beams for calculation), finds the beams listed in 
    IVBEAMS (if not passed as a list 'beams', will attempt to call readIVBEAMS directly) and copies 
    the corresponding lines to AUXBEAMS. Returns a list of the corresponding beam numbers in _BEAMLIST."""
    if not readfile == "_BEAMLIST":
        blstr = "Beam list (filename "+readfile+")"
    else:
        blstr = "_BEAMLIST"
        
    if beams == []:         #if 'beams' is empty, try to fill it
        beams = readIVBEAMS(beamsfile)
        
    output = '   1               IFORM\n'      # !!! WHAT IS THIS VALUE? WHERE TO GET IT FROM?
    
    # read BEAMLIST
    try:
        rf = open(readfile, 'r')
    except FileNotFoundError:
        logging.error(blstr+" not found.")
        raise
    
    err = 1e-4           #since beams are saved as floats, give some error tolerance when comparing e.g. 1/3
    #read BEAMLIST - very little error handling here, I'm assuming the beamlists are safe; could be added later
    bc = beams[:]
    numlist = []
    blocks = 0
    totalbeams = 0
    for line in rf:
        llist = linelist(line)
        if len(llist) > 1:
            totalbeams += 1
            for b in bc:
                if abs(b[0]-float(llist[0]))<err and abs(b[1]-float(llist[1]))<err:
                    output += line
                    bc.remove(b)
                    numlist.append(int(line.split('.')[-1]))
        else:
            blocks += 1
    rf.close()
    logging.info(blstr+" was read successfully")
    for beam in bc:
        logging.warning('IVBEAMS contains beam '+str(beam)+', which was not found in '+blstr)
    #write output
    if write:
        if not writefile == 'AUXBEAMS':
            wfstr = 'AUXBEAMS file (filename '+writefile+')'
        else:
            wfstr = 'AUXBEAMS'
        try:
            wf = open(writefile, 'w')
            wf.write(output)
            wf.close()
        except:
            logging.error("Failed to write "+wfstr)
            raise
        logging.info("Wrote to "+wfstr+" successfully.")
    return numlist, blocks, totalbeams
    
def readPHASESHIFTS(readfile='_PHASESHIFTS'):
    """"Reads from a _PHASESHIFTS file, returns the data as a list of tuples (E, enps),
    where enps is a list of lists, containing one list of values (for different L) for each
    element. Therefore, len(phaseshifts) is the number of energies found, len(phaseshifts[0][1]) 
    should match the number of elements, and len(phaseshifts[0][1][0]) is the number of different
    values of L in the phaseshift file."""
    rf74x10 = ff.FortranRecordReader('10F7.4')
    ri3 = ff.FortranRecordReader('I3')
    
    try:
        rf = open(readfile, 'r')
    except FileNotFoundError:
        logging.error("_PHASESHIFTS file not found.")
        raise
        
    filelines = []
    for line in rf:
        filelines.append(line[:-1])
    rf.close()
    
    NEL = ri3.read(filelines[0])[0]
    phaseshifts = []
    
    readline = 1
    linesperblock = 0
    while readline < len(filelines):
        if linesperblock:
            En = rf74x10.read(filelines[readline])[0]
            enps = []
            for i in range(0,NEL):
                elps = []
                for j in range(0,linesperblock):
                    llist = rf74x10.read(filelines[readline+(i*linesperblock)+j+1])
                    llist = [f for f in llist if f is not None]
                    elps.extend(llist)
                enps.append(elps)
            phaseshifts.append((En,enps))
            readline += linesperblock*NEL+1
        else:
            #first check how many lines until the next energy:
            lineit = 1
            llist = rf74x10.read(filelines[readline+lineit])
            llist = [f for f in llist if f is not None]
            longestline = len(llist)
            shortestline = longestline
            lastlen = longestline
            cont = True
            while cont:
                lineit += 1
                llist = rf74x10.read(filelines[readline+lineit])
                llist = [f for f in llist if f is not None]
                if len(llist) == 1:
                    if lastlen == 1 or (shortestline > 1 and shortestline < longestline):   #found next energy
                        cont = False
                    else:
                        shortestline = 1
                        #blocklines += 1
                elif len(llist) != longestline:
                    shortestline = len(llist)
                lastlen = len(llist)
            linesperblock = int((lineit-1)/NEL)
            # don't increase readline -> read the same block again afterwards
    return phaseshifts

def reduceUnitCell(ab, eps = 0.0001):
    """Takes an obtuse unit cell as a (2x2) matrix and reduces it to minimum circumference, keeping the area constant.
    This might reduce oblique unit cells to rectangular or hexagonal ones. Returns the modified unit cell, as well as the
    transformation matrix."""
    u = np.array([[1,0],[0,1]])
    brk = False
    while abs(np.vdot(ab[0],ab[1])) > eps and np.vdot(ab[0],ab[1]) < 0 and not brk:
        if np.linalg.norm(ab[0]) > np.linalg.norm(ab[1]):
            m = np.array([[1,1],[0,1]])
        else:
            m = np.array([[1,0],[1,1]])
        tmp = np.dot(m,ab)
        if np.vdot(tmp[0],tmp[1]) <= 0 or abs(np.vdot(tmp[0],tmp[1])) < eps:
            u = np.dot(m, u)
            ab = np.dot(m, ab)
        else:
            brk = True
    return ab, u

        