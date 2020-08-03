# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

Class accumulating atoms and layers, listing their elements and various other properties. Includes functions for
manipulation of those properties.
"""

from timeit import default_timer as timer

import logging
import numpy as np
import copy
import re
import multiprocessing as mp
from scipy.spatial.distance import cdist

import TensErLeedModules as tl
from TensErLeedModules import DEFAULT


class SymPlane:
    """Candidate plane for a symmetry operation. 'ty' pre-defines a type (mirror or glide), 'index2' allows the (1,2) and (2,1)
    directions if True, and collapse moves pos into the (0,0) unit cell if True."""
    def __init__(self, pos, dr, abt, ty="none", index2=False, collapse=True):
        if collapse:
            self.pos = np.dot(np.transpose(abt), (np.dot(np.linalg.inv(np.transpose(abt)), pos) % 1.0)) #collapsed to (0,0) cell
        else:
            self.pos = pos
        self.dir = dr/np.linalg.norm(dr)   #normalized vector perpendicular to pos = in-plane
        self.type = ty
        self.par = []
        optionlist = [(1,0),(0,1),(1,1),(1,-1)]
        if index2: optionlist.extend([(2,1),(1,2)])
        for (i,j) in optionlist:
            if abs((abs(np.dot(self.dir,(i*abt[0]+j*abt[1])))/(np.linalg.norm(self.dir)*np.linalg.norm(i*abt[0]+j*abt[1])))-1.0) < 0.001:
                self.par = np.array([i,j])

    def isEquivalent(self,pl2,abt,eps=0.001):
        """Checks whether two symmetry planes have the same position and direction (including duplicates in next unit cell)"""
        if not np.array_equal(self.par,pl2.par): return False
        complist = [self.pos]
        fpos = np.dot(np.linalg.inv(np.transpose(abt)), self.pos) % 1.0
        for i in range(0,2):    # if we're close to an edge or corner, also check translations to the other edges
            releps = eps / np.linalg.norm(abt[i])
            if abs(fpos[i]) < releps:
                complist.append(self.pos+abt[i])
            if abs(fpos[i]-1) < releps:
                complist.append(self.pos-abt[i])
        if len(complist) == 3:  # coner - add the diagonally opposed one
            complist.append(complist[1]+complist[2]-complist[0])

        for p in complist:
            if tl.distanceLineThroughPointsFromPoint(pl2.pos,pl2.pos+pl2.dir,p) < eps:
                return True
        return False
     
        
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
        self.sublayers = []     # list of Layer objects, each containing atoms of equal element and Z coordinate
        self.sitelist = []      # list of distinct sites as Sitetype elements
        self.uCellTransform = np.array([[1,0],[0,1]])   # stores a,b transformations of the unit cell by multiplication
        self.uCellShift = np.array([0.0,0.0])       # stores translations of the unit cell
        self.topatOriZ = 0.0      # stores the original position of the topmost atom in cartesian coordinates
        self.planegroup = "unknown"     # symmetry group of the slab, as string
        self.orisymplane = DEFAULT      # only stored if the planegroup is ambigious as to which unit vector the symmetry plane at the origin is parallel to
        self.linklists = []             # list of lists of atoms which are linked by a symmetry operation
    
    def fullUpdate(self,rparams):
        """readPOSCAR initializes the slab with information from POSCAR; fullUpdate re-initializes the atom list,
        then uses the information from the parameters file to create layers, calculate cartesian coordinates 
        (absolute and per layer), and to update elements and sites."""
        #self.initAtomList()
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
            self.atlist.append(tl.Atom(self.elements[elnum],self.atpos[i],i+1,self))
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
        self.topatOriZ = topcart[2]
        for atom in al:
            atom.cartpos = np.dot(self.ucell, atom.pos)
            atom.cartpos[2] = topcart[2]-atom.cartpos[2]
#            logging.debug(atom.el+' '+str(atom.oriN)+' '+str(atom.pos)+' '+str(atom.cartpos))             
            
    def getFractionalCoordinates(self):
        """Calculates fractional coordinates for all atoms from their cartesian coordinates, using the slab unit cell."""
        uci = np.linalg.inv(self.ucell)
        for at in self.atlist:
            tp = np.copy(at.cartpos)
            tp[2] = self.topatOriZ-tp[2]
            at.pos = np.dot(uci,tp)

    def collapseFractionalCoordinates(self):
        """Finds atoms outside the parallelogram spanned by the unit vectors a and b and moves them inside."""
        for at in self.atlist:
            at.pos[0] = at.pos[0] % 1.0
            at.pos[1] = at.pos[1] % 1.0
        
    def collapseCartesianCoordinates(self):
        """Finds atoms outside the parallelogram spanned by the unit vectors a and b and moves them inside."""
        self.getFractionalCoordinates()
        self.collapseFractionalCoordinates()
        self.getCartesianCoordinates()
                
    def createLayers(self,rparams):
        """Creates a list of Layer objects based on the BLAY and CTRUNC parameters in rparams.
        If layers were already defined, overwrite."""
        self.layers = []
        tmplist = self.atlist[:]
        self.sortByZ()
        laynum = 0
        b = True if rparams.BLAY > 0 else False
        newlayer = tl.Layer(self,0,b)
        self.layers.append(newlayer)
        for atom in self.atlist:
            if laynum < len(rparams.CTRUNC):        #only check for new layer if we're not in the top layer already
                if atom.pos[2] > rparams.CTRUNC[laynum]:   #if the atom is higher than the next c cutoff, make a new layer
                    laynum += 1
                    b = True if rparams.BLAY > laynum else False
                    newlayer = tl.Layer(self,laynum,b)
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
                            newlayer = tl.Layer(self,laynum,b)
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
        
    def createSublayers(self,eps=0.001):
        """Sorts the atoms in the slab into sublayers, sorted by element and Z coordinate."""
        self.sortByZ()
        self.sublayers = []  #will be a list of sublayers, using the Layer class, where a sublayer is all atoms of the same element at the same z position (+- eps)
        for el in self.elements:
            currentz = -1
            for at in [a for a in self.atlist if a.el == el]:
                if currentz < 0 or abs(at.cartpos[2]-currentz) > eps:
                    newsl = tl.Layer(self, len(self.sublayers),sublayer=True)
                    self.sublayers.append(newsl)
                    newsl.cartbotz = at.cartpos[2]
                    currentz = at.cartpos[2]
                    newsl.atlist.append(at)
                else:
                    newsl.cartbotz = at.cartpos[2]
                    newsl.atlist.append(at)
        self.sublayers.sort(key=lambda sl: sl.cartbotz)     # reorder sublayers by z
        
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
                newsite = tl.Sitetype(el, sitename)
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
            newsite = tl.Sitetype(el, 'def')
            found = False
            for at in atlist:
                if at.el == el and at.site == DEFAULT:
                    at.site = newsite
                    found = True
            if found: sl.append(newsite)
        self.sitelist = sl
      
    def sortByZ(self, botToTop=False):
        """Sorts atlist by z coordinate"""
        self.atlist.sort(key=lambda atom: atom.pos[2])
        if botToTop:
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
            alpha = tl.angle(ab[:,0],ab[:,1])
            for at in self.atlist:
                xyCart = at.cartpos[0:2]
                for i in range(0,2):
                    proj = (np.dot(xyCart,ab[:,i]) / (np.linalg.norm(ab[:,i])**2)) * ab[:,i]
                    at.pos[i] = ((np.dot(proj,ab[:,i])/(np.linalg.norm(ab[:,i])))-(np.linalg.norm(xyCart-proj)/np.tan(alpha))) / np.linalg.norm(ab[:,i]) % 1.0
            self.ucell[0,2] = 0
            self.ucell[1,2] = 0
            self.getCartesianCoordinates()

    def rotate(self, axis, order):
        """Translates the atoms in the slab to have the axis in the origin, applies an order-fold rotation matrix, then translates back"""
        if order == 2:              #these explicit definitions are likely useless, but sqrts might be marginally more accurate than sin/cos
            m = np.array([[-1,0],[0,-1]])
        elif order == 3:
            m = np.array([[-0.5,-np.sqrt(3)/2],[np.sqrt(3)/2,-0.5]])
        elif order == -3:
            m = np.array([[-0.5,np.sqrt(3)/2],[-np.sqrt(3)/2,-0.5]])
        elif order == 4:
            m = np.array([[0,1],[-1,0]])
        elif order == -4:
            m = np.array([[0,-1],[1,0]])
        elif order == 6:
            m = np.array([[0.5,np.sqrt(3)/2],[-np.sqrt(3)/2,0.5]])
        elif order == -6:
            m = np.array([[0.5,-np.sqrt(3)/2],[np.sqrt(3)/2,0.5]])
        else:
            angle = 2*np.pi/order
            m = np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
        for at in self.atlist:
            at.cartpos[0:2] -= axis    # translate origin to candidate point
            at.cartpos[0:2] = np.dot(m, at.cartpos[0:2])    # rotation
            at.cartpos[0:2] += axis    # undo translation
            
    def mirror(self, symplane, glide=False):
        """Translates the atoms in the slab to have the symplane in the origin, applies a mirror or glide matrix, then translates back"""
        angle = tl.angle(np.array([1,0]),symplane.dir)
        if symplane.dir[1] > 0: angle *= -1
        rotm = np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
        rotmirm = np.dot(np.linalg.inv(rotm),np.dot(np.array([[1,0],[0,-1]]),rotm))     #rotate to have plane in x direction, mirror on x
        if glide:
            abt = np.transpose(self.ucell[0:2,0:2])
            glidevec = (symplane.par[0]*abt[0]+symplane.par[1]*abt[1])/2
        for at in self.atlist:
            at.cartpos[0:2] -= symplane.pos     #translate to plane
            at.cartpos[0:2] = np.dot(rotmirm, at.cartpos[0:2])    #apply mirror
            at.cartpos[0:2] += symplane.pos     #translate back
            if glide: at.cartpos[0:2] += glidevec        # if we're testing for glide plane, add the appropriate vector
            
    def isRotationSymmetric(self,axis,order,eps):
        """Evaluates whether the slab is equivalent to itself when rotated around the axis with the given rotational order"""
#        logging.debug("##### starting rotation symmetry check #####\nOrder: "+str(order))
        if order == 2:              #these explicit definitions are likely useless, but sqrts might be marginally more accurate than sin/cos
            m = np.array([[-1,0],[0,-1]])
        elif order == 3:
            m = np.array([[-0.5,-np.sqrt(3)/2],[np.sqrt(3)/2,-0.5]])
        elif order == -3:
            m = np.array([[-0.5,np.sqrt(3)/2],[-np.sqrt(3)/2,-0.5]])
        elif order == 4:
            m = np.array([[0,1],[-1,0]])
        elif order == -4:
            m = np.array([[0,-1],[1,0]])
        elif order == 6:
            m = np.array([[0.5,np.sqrt(3)/2],[-np.sqrt(3)/2,0.5]])
        elif order == -6:
            m = np.array([[0.5,-np.sqrt(3)/2],[np.sqrt(3)/2,0.5]])
        else:
            angle = 2*np.pi/order
            m = np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
        ab = self.ucell[0:2,0:2]
        abt = np.transpose(ab)
        releps = np.array([0.0,0.0])
        for j in range(0,2):
            releps[j] = eps / np.linalg.norm(abt[j])
        shiftv = axis.reshape(2,1)
        for sl in self.sublayers:
            coordlist = []
            for at in sl.atlist:
                coordlist.append(at.cartpos[0:2])
            shiftm = np.tile(shiftv,len(coordlist))             # matrix to shift all coordinates by axis
            oricm = np.array(coordlist)             # original cartesian coordinate matrix
            oripm = np.dot(np.linalg.inv(ab), oricm.transpose()) % 1.0  # collapse (relative) coordinates to base unit cell
            oricm = np.dot(ab, oripm).transpose()       # original cartesian coordinates collapsed to base unit cell
            tmpcoords = np.copy(oricm).transpose()                          # copy of coordinate matrix to be rotated
            tmpcoords -= shiftm
            tmpcoords = np.dot(m,tmpcoords)
            tmpcoords += shiftm
            tmpcoords = np.dot(ab, (np.dot(np.linalg.inv(ab),tmpcoords) % 1.0))  # collapse coordinates to base unit cell
            # for every point in matrix, check whether is equal
            for (i,p) in enumerate(oripm.transpose()):     # get extended comparison list for edges/corners
                addlist = []
                for j in range(0,2):
                    if abs(p[j]) < releps[j]:
                        addlist.append(oricm[i]+abt[j])
                    if abs(p[j]-1) < releps[j]:
                        addlist.append(oricm[i]-abt[j])
                if len(addlist) == 2:
                    addlist.append(addlist[0]+addlist[1]-oricm[i])    # coner - add the diagonally opposed one
                for v in addlist:
                    oricm = np.concatenate((oricm,v.reshape(1,2)))
            distances = cdist(tmpcoords.transpose(), oricm, 'euclidean')
            for sublist in distances:
                if min(sublist) > eps:
                    return False
        return True
                
    def isMirrorSymmetric(self, symplane, eps, glide=False):
        """Evaluates whether the slab is equivalent to itself when applying a mirror or glide operation at a given plane"""
        angle = tl.angle(np.array([1,0]),symplane.dir)
        if symplane.dir[1] > 0: angle *= -1
        rotm = np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
        rotmirm = np.dot(np.linalg.inv(rotm),np.dot(np.array([[1,0],[0,-1]]),rotm))     #rotate to have plane in x direction, mirror on x
        ab = self.ucell[0:2,0:2]
        abt = np.transpose(ab)
        releps = np.array([0.0,0.0])
        for j in range(0,2):
            releps[j] = eps / np.linalg.norm(abt[j])
        shiftv = symplane.pos.reshape(2,1)
        if glide:
            glidev = ((symplane.par[0]*abt[0]+symplane.par[1]*abt[1])/2).reshape(2,1)
        for sl in self.sublayers:
            coordlist = []
            for at in sl.atlist:
                coordlist.append(at.cartpos[0:2])
            shiftm = np.tile(shiftv,len(coordlist))             # matrix to shift all coordinates by 
            if glide: glidem = np.tile(glidev,len(coordlist))   # matrix to shift all coordinates by glidev
            oricm = np.array(coordlist)             # original cartesian coordinate matrix
            oripm = np.dot(np.linalg.inv(ab), oricm.transpose()) % 1.0  # collapse (relative) coordinates to base unit cell
            oricm = np.dot(ab, oripm).transpose()       # original cartesian coordinates collapsed to base unit cell
            tmpcoords = np.copy(oricm).transpose()              # copy of coordinate matrix to be rotated
            tmpcoords -= shiftm
            tmpcoords = np.dot(rotmirm,tmpcoords)
            tmpcoords += shiftm
            if glide: tmpcoords += glidem
            tmpcoords = np.dot(ab, (np.dot(np.linalg.inv(ab),tmpcoords) % 1.0))  # collapse coordinates to base unit cell
            # for every point in matrix, check whether is equal
            for (i,p) in enumerate(oripm.transpose()):     # get extended comparison list for edges/corners
                addlist = []
                for j in range(0,2):
                    if abs(p[j]) < releps[j]:
                        addlist.append(oricm[i]+abt[j])
                    if abs(p[j]-1) < releps[j]:
                        addlist.append(oricm[i]-abt[j])
                if len(addlist) == 2:
                    addlist.append(addlist[0]+addlist[1]-oricm[i])    # coner - add the diagonally opposed one
                for v in addlist:
                    oricm = np.concatenate((oricm,v.reshape(1,2)))
            distances = cdist(tmpcoords.transpose(), oricm, 'euclidean')
            for (i,sublist) in enumerate(distances):    # TODO: remove i once debugging is done
                if min(sublist) > eps:
                    return False
        return True
            
            
    def isEquivalent(self,slab,eps=0.001):
        """Compares the slab to another slab, returns True if all atom cartpos match (with at least one 
        other atom, if there are duplicates), False if not. Both slabs are copied and collapsed to the (0,0) cell before."""
        slab1 = copy.deepcopy(self)
        slab2 = copy.deepcopy(slab)
        slab1.collapseCartesianCoordinates()
        slab2.collapseCartesianCoordinates()
        # reorder sublayers by Z to then compare by index
        slab1.sublayers.sort(key=lambda sl: sl.cartbotz)
        slab2.sublayers.sort(key=lambda sl: sl.cartbotz)
        ab = self.ucell[0:2,0:2]
        for (i,sl) in enumerate(slab1.sublayers):
            if len(sl.atlist) != len(slab2.sublayers[i].atlist) or abs(sl.cartbotz-slab2.sublayers[i].cartbotz) > eps or sl.atlist[0].el != slab2.sublayers[i].atlist[0].el: return False
            for at1 in sl.atlist:
                complist = [at1.cartpos[0:2]]    # if we're close to an edge or corner, also check translations to the other edges
                for j in range(0,2):
                    releps = eps / np.linalg.norm(ab[:,j])
                    if abs(at1.pos[j]) < releps:
                        complist.append(at1.cartpos[0:2]+ab[:,j])
                    if abs(at1.pos[j]-1) < releps:
                        complist.append(at1.cartpos[0:2]-ab[:,j])
                if len(complist) == 3:
                    complist.append(complist[1]+complist[2]-complist[0])    # coner - add the diagonally opposed one
                found = False
                for at2 in slab2.sublayers[i].atlist:
                    for p in complist:
                        if np.linalg.norm(p-at2.cartpos[0:2]) < eps: 
                            found = True
                            break
                    if found: break
                if not found:
                    return False
        return True
    
    def revertUnitCell(self):
        """If the unit cell in a and b was transformed earlier, restore the original form and coordinates."""
        # TODO - implement the option to do this before writing as a parameter
        # TODO - test
        for at in self.atlist:
            at.cartpos += self.uCellShift
        self.getFractionalCoordinates()
        self.ucell = np.dot(self.ucell, np.linalg.inv(np.transpose(self.uCellTransform)))     #same as self.ucell = np.transpose(np.dot(np.linalg.inv(self.uCellTransform),np.transpose(self.ucell)))
        for at in self.atlist:
            at.pos = np.dot(np.transpose(self.uCellTransform),at.pos)
        self.getCartesianCoordinates()

    def findSymmetry(self, rp):
        """Reduces the unit cell if necessary and finds the plane group of the slab. Stores the plane group and
        the higher-symmetry direction of the unit cell, if there is one."""
        celltype = "ERROR - not recognized"     #storing this as string is computationally inefficient, but readability seems more important than speed here
        planegroup = ""     #plane group will be stored in Hermann-Mauguin notation
        eps = rp.SYMMETRY_EPS
        epsz = rp.SYMMETRY_EPS_Z
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
            (abst, u) = tl.reduceUnitCell(abst, eps)
            usurf = np.dot(u,usurf)
        # reduce bulk unit cell
        abbt = np.dot(np.linalg.inv(rp.SUPERLATTICE),abst)      #bulk ab unit cell, transposed
        ubulk = np.array([[1,0],[0,1]])                       #track unit cell changes
        dp = np.dot(abbt[0],abbt[1])
        if abs(dp) > eps:    #rhombic, hexagonal or oblique
            if dp > 0:      #transform acute to obtuse
                abbt = np.dot(np.array([[-1,0],[0,1]]), abbt)
                ubulk = np.dot(np.array([[-1,0],[0,1]]), ubulk)
            (abbt, u) = tl.reduceUnitCell(abbt, eps)
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
        elif tl.angle(abst[0],abst[1]) - (2*np.pi/3) < eps:
            celltype = "hexagonal"
        else:
            celltype = "rhombic"
        logging.info("Found unit cell type: "+celltype)
        # FIND HIGHEST SYMMETRY ORIGIN
        logging.debug("Initializing symmetry search...")
        self.collapseCartesianCoordinates()
        # create a testslab: C projected to Z
        ts = copy.deepcopy(self)
        ts.projectCToZ()
        ts.sortByZ()
        ts.collapseCartesianCoordinates()
        bigslab = copy.deepcopy(ts)     # will have atoms duplicated and shifted to 4 unit cells ([0,0], [0,1], [1,0], [1,1])
        tmplist = bigslab.atlist[:]
        for at in tmplist:
            for i in range(0,2):
                for j in range(0,2):
                    if not (i==0 and j==0):
                        tmpat = at.duplicate()
                        tmpat.pos[0] += i
                        tmpat.pos[1] += j
        bigslab.getCartesianCoordinates()
        bigslab.fullUpdate(rp)
        bigslab.createSublayers(epsz)
        
        # find the lowest occupancy sublayer; comparing candidate axes / planes to this one will be fastest
        minlen = len(bigslab.sublayers[0].atlist)
        lowocclayer = bigslab.sublayers[0]
        for sl in bigslab.sublayers:
            if len(sl.atlist) < minlen:
                lowocclayer = sl
                minlen = len(sl.atlist)
        
        # find candidate positions for symmetry points / planes:
        
#        c = 0
#        starttime = timer()
#        for sl in bigslab.sublayers:
#            endtime = timer()
#            logging.debug("Time elapsed: "+str(endtime-starttime))
#            starttime=endtime
#            logging.debug("Sublayer "+str(c)+" ("+str(len(sl.atlist))+" atoms)")
#            c += 1
        if not rp.SYMMETRY_USEORI:
            logging.debug("Generating candidate high-symmetry positions from layer with "+str(minlen)+" atoms...")
            pool = mp.Pool(mp.cpu_count())  #finding the candidate positions can be time expensive, but generating sublists per atom can be parallelized
            result = pool.starmap_async(getAtomSymposlist, [(lowocclayer.atlist,i,eps,celltype) for i in range(0,len(lowocclayer.atlist))])
            pool.close()
            pool.join()
            resultlist = result.get()
            symposlist = [np.array([0.0,0.0])]    #always check the origin, even if it was not found.
            for sublist in resultlist:
                symposlist = tl.addUnequalPoints(sublist, symposlist, eps, uniqueLists=True)
        else:
            symposlist = [np.array([0.0,0.0])]  #only check origin, planes are defined explicitly for this case
   
        # we're done with the bigger slab, actually testing symmetry operations can be done just on the basic one.
        ts.createSublayers(epsz)
        lowocclayer = ts.sublayers[bigslab.sublayers.index(lowocclayer)]
        del bigslab
        
#        logging.debug("Sublayers found: "+str(len(ts.sublayers)))
#        logging.debug(symposlist)
        
        #find potential rotation axes:
        if not rp.SYMMETRY_USEORI:
            logging.debug("Checking for rotation axes: "+str(len(symposlist))+" candidates...")
#        spshared = lowocclayer.symposlist[:]
#        for sl in ts.sublayers:
#            if not sl == lowocclayer:
#                for p1 in spshared:
#                    found = False
#                    complist = [p1]
#                    for i in range(0,2):    # if we're close to an edge or corner, also check translations to the other edges
#                        releps = eps / np.linalg.norm(abst[i])
#                        if abs(p1[i]) < releps:
#                            complist.append(p1+abst[i])
#                        if abs(p1[i]-1) < releps:
#                            complist.append(p1-abst[i])
#                    if len(complist) == 3:   # coner - add the diagonally opposed one
#                        complist.append(complist[1]+complist[2]-complist[0])
#                    for p2 in sl.symposlist:
#                        for p3 in complist:
#                            if np.linalg.norm(p2-p3) < eps:     # counts as same position
#                                found = True
#                                break
#                        if found: break
#                    if not found: 
#                        i = 0
#                        while(i<len(spshared)):     # there might be multiple occurances of the same point in one sublayer
#                            if np.array_equal(p1,spshared[i]):
#                                spshared.pop(i)
#                            else:
#                                i += 1
    
        #test potential rotation axes:
        toprotsym = 0   # keep track of the highest rotational symmetry so far
        topsympoint = symposlist[0]

        for p in symposlist:
            rotsymorder = 0
            if ts.isRotationSymmetric(p,2,eps):
                rotsymorder = 2
                if celltype == "square":
                    if ts.isRotationSymmetric(p,4,eps):
                        rotsymorder = 4
            if celltype == "hexagonal":
                if ts.isRotationSymmetric(p,3,eps):
                    rotsymorder = 3
                    if ts.isRotationSymmetric(p,6,eps):
                        rotsymorder = 6
            if rotsymorder > toprotsym:     #new best point found
                topsympoint = p
                toprotsym = rotsymorder
            elif rotsymorder != 0 and rotsymorder == toprotsym:
                if np.linalg.norm(p-abst[0]-abst[1]) < np.linalg.norm(topsympoint-abst[0]-abst[1]):     # for points with equal rotational symmetry, prioritize the one closer to the origin (of cell 1,1) to avoid shifting the unit cell randomly every time
                    topsympoint = p
        if toprotsym > 0:
            #shift origin to highest rotation axis
            for at in self.atlist:
                at.cartpos[0:2] -= topsympoint
            for at in ts.atlist:
                at.cartpos[0:2] -= topsympoint
            self.uCellShift -= topsympoint
            self.getFractionalCoordinates()
            ts.getFractionalCoordinates()
            logging.debug('Highest rotation axis has order '+str(toprotsym))
            
        if toprotsym == 0:
            logging.debug("Checking for mirror/glide planes...")
            #check for mirror/glides
            mirror = False
            glide = False
            symplanelist = []
#            counter = 0
#            counter2 = 0
            
            if rp.SYMMETRY_USEORI:
                for (i,j) in [(1,0),(0,1),(1,1),(1,-1)]:
                    symplanelist.append(SymPlane(np.array([0.0,0.0]),i*abst[0]+j*abst[1],abst))
            else:
                for (i,sposi) in enumerate(symposlist):     
                    #checking only lowest occupation layer - this means more symmetry operations have to be tested later, but we don't 
                    #have to check candidates in every layer. Should scale better, maybe worth testing the other option.
                    posplanelist = []   #planes that the first point is on; others don't need checking
                    pos = np.dot(self.ucell[0:2,0:2], (np.dot(np.linalg.inv(self.ucell[0:2,0:2]), sposi) % 1.0))    #collapse
                    for osp in symplanelist:
                        if tl.distanceLineThroughPointsFromPoint(osp.pos,osp.pos+osp.dir,pos) < eps:
                            posplanelist.append(osp)
                    for (j, sposj) in [(j, sposj) for (j, sposj) in enumerate(symposlist) if j>i]:
                        if np.linalg.norm(sposj-sposi) > eps:
                            dr = sposj-sposi
                            found = False
                            for osp in posplanelist:
                                if abs((abs(np.dot(dr,osp.dir))/np.linalg.norm(dr))-1) < 0.001:     #same direction; eps is too large here!
                                    found = True    #same direction, crossing the same point -> same plane
                                    break
                            if not found:
    #                            counter += 1
                                newsymplane = SymPlane(pos,dr,abst)     #now do the more expensive check of the actual plane...
                                if len(newsymplane.par) == 2:   #checks whether the plane is parallel to any of the relevant axes
    #                                counter2 += 1
                                    for osp in symplanelist:
                                        if osp.isEquivalent(newsymplane,abst,eps):
                                            found = True
                                            break
                                    if not found: symplanelist.append(newsymplane)
#            logging.debug(str(counter)+" candidates tested more closely")
#            logging.debug(str(counter2)+" candidates parallel")
            logging.debug(str(len(symplanelist))+" candidates for mirror/glide planes found...")
            for spl in symplanelist:    #test the candidates
#                tmpslab = copy.deepcopy(ts)
#                tmpslab.mirror(spl)
#                if tmpslab.isEquivalent(ts,eps):
                if ts.isMirrorSymmetric(spl,eps):
                    spl.type = "mirror"
                    mirror = True
                elif ts.isMirrorSymmetric(spl,eps,glide=True):
                    spl.type = "glide"
                    glide = True
#                else:   #test for glide
#                    glidevec = (spl.par[0]*abst[0]+spl.par[1]*abst[1])/2
#                    for at in tmpslab.atlist:
#                        at.cartpos[0:2] += glidevec
#                    if tmpslab.isEquivalent(ts,eps):
#                        spl.type = "glide"
#                        glide = True
            i = 0
            while i < len(symplanelist):
                if symplanelist[i].type == "none":
                    symplanelist.pop(i)
                else:
                    i += 1
            #logging.debug([(sp.type, sp.pos) for sp in symplanelist])
            # identify group:
            oriplane = DEFAULT
            if not mirror:
                if not glide:
                    planegroup = "p1"   #no origin shift required.
                else:
                    planegroup = "pg"
                    for spl in symplanelist:
                        if oriplane == DEFAULT:
                            oriplane = spl
                        elif np.linalg.norm(spl.pos) < np.linalg.norm(oriplane.pos) or abs(np.linalg.norm(spl.pos)-np.linalg.norm(spl.par[0]*abst[0]+spl.par[1]*abst[1])) < np.linalg.norm(oriplane.pos):
                            oriplane = spl  #prioritize planes close to the origin
            else:
                if not glide:
                    planegroup = "pm"
                else:
                    planegroup = "cm"
                for spl in symplanelist:
                    if spl.type == "mirror":
                        if oriplane == DEFAULT:
                            oriplane = spl
                        elif np.linalg.norm(spl.pos) < np.linalg.norm(oriplane.pos) or abs(np.linalg.norm(spl.pos)-np.linalg.norm(spl.par[0]*abst[0]+spl.par[1]*abst[1])) < np.linalg.norm(oriplane.pos):
                            oriplane = spl  #prioritize planes close to the origin
                if planegroup == "cm":  #both mirrors and glides. if parallel to unit vectors -> rcm 
                    if tuple(oriplane.par) in [(1,0),(0,1)]:
                        planegroup = "rcm"
            if oriplane != DEFAULT:
                # shift to closest point on oriplane
                shiftv = np.array([oriplane.dir[1], -oriplane.dir[0]])*tl.distanceLineThroughPointsFromPoint(oriplane.pos,oriplane.pos+oriplane.dir,np.array([0,0]))
                if tl.distanceLineThroughPointsFromPoint(oriplane.pos,oriplane.pos+oriplane.dir,shiftv) > eps: shiftv = -1*shiftv
                for at in self.atlist:
                    at.cartpos[0:2] -= shiftv
#                for at in ts.atlist:
#                    at.cartpos[0:2] -= oriplane.pos
       #ts is not used any more in this case, otherwise those atoms would have to be shifted as well.
                self.uCellShift -= shiftv
                self.getFractionalCoordinates()
#                ts.getFractionalCoordinates()
                oriplane.pos = np.array([0,0])
                self.orisymplane = oriplane

        if not planegroup:
            if toprotsym == 2:  #start by checking special case: in cmm, there are two inequivalent 2fold axes, one of which would not have been found yet -> shift there (potentially), test
                shiftslab = copy.deepcopy(ts)
                for at in shiftslab.atlist:
                    at.cartpos[0:2] -= abst[0]/2
                shiftslab.getFractionalCoordinates()
                #test diagonal mirror at shifted origin
                spl = SymPlane(np.array([0,0]), (abst[0]+abst[1]), abst)
#                tmpslab = copy.deepcopy(shiftslab)
#                tmpslab.mirror(spl)
#                if tmpslab.isEquivalent(shiftslab,eps):
                if shiftslab.isMirrorSymmetric(spl,eps):
                    planegroup = "cmm"
                    ts = shiftslab
                    #correct origin
                    for at in self.atlist:
                        at.cartpos[0:2] -= abst[0]/2
                    self.uCellShift -= abst[0]/2
                    self.getFractionalCoordinates()

        if not planegroup:
            efftype = ""    #effective cell type
            if celltype == "hexagonal":
                if not (toprotsym == 3 or toprotsym == 6):
                    efftype = "rhombic"
                else:
                    #test mirror plane along unit vector at origin
                    spl = SymPlane(np.array([0,0]), abst[0], abst)
#                    tmpslab = copy.deepcopy(ts)
#                    tmpslab.mirror(spl)
#                    if tmpslab.isEquivalent(ts,eps):
                    if ts.isMirrorSymmetric(spl,eps):
                        if toprotsym == 6:
                            planegroup = "p6m"
                        else:
                            planegroup = "p31m"
                    else:
                        if toprotsym == 6:
                            planegroup = "p6"
                        else:
                            #test mirror plane along both diagonals at origin
                            found = False
                            for i in [+1,-1]:
                                spl = SymPlane(np.array([0,0]), (abst[0]+(i*abst[1])), abst)
#                                tmpslab = copy.deepcopy(ts)
#                                tmpslab.mirror(spl)
#                                if tmpslab.isEquivalent(ts,eps):
                                if ts.isMirrorSymmetric(spl,eps):
                                    found = True
                                    break
                            if found:
                                planegroup = "p3m1"
                            else:
                                planegroup = "p3"
            elif celltype == "square":
                if toprotsym == 4:
                    #test mirror plane along unit vector at origin
                    spl = SymPlane(np.array([0,0]), abst[0], abst)
#                    tmpslab = copy.deepcopy(ts)
#                    tmpslab.mirror(spl)
#                    if tmpslab.isEquivalent(ts,eps):
                    if ts.isMirrorSymmetric(spl,eps):
                        planegroup = "p4m"
                    else:
                        #test glide plane along diagonal at origin
                        spl = SymPlane(np.array([0,0]), (abst[0]+abst[1]), abst)
#                        tmpslab = copy.deepcopy(ts)
#                        tmpslab.mirror(spl, glide=True)
#                        if tmpslab.isEquivalent(ts,eps):
                        if ts.isMirrorSymmetric(spl,eps, glide=True):
                            planegroup = "p4g"
                        else:
                            planegroup = "p4"
                else:   # p1 p2 pm pg pmm pmg pgg
                    # test mirror plane along both diagonals at origin
                    found = False
                    for i in [+1,-1]:
                        spl = SymPlane(np.array([0,0]), (abst[0]+i*abst[1]), abst)
#                        tmpslab = copy.deepcopy(ts)
#                        tmpslab.mirror(spl)
#                        if tmpslab.isEquivalent(ts,eps):
                        if ts.isMirrorSymmetric(spl,eps):
                            found = True
                            break
                    if found:
                        if toprotsym == 2:
                            planegroup = "cmm"
                        else:
                            planegroup = "cm"
                    else:
                        efftype = "rectangular"
            if celltype == "rhombic" or efftype == "rhombic":
                # test mirror plane along both diagonals at origin
                found = False
                for i in [+1,-1]:
                    spl = SymPlane(np.array([0,0]), (abst[0]+i*abst[1]), abst)
#                    tmpslab = copy.deepcopy(ts)
#                    tmpslab.mirror(spl)
#                    if tmpslab.isEquivalent(ts,eps):
                    if ts.isMirrorSymmetric(spl,eps):
                        found = True
                        break
                if not found:
                    efftype = "oblique"     #since we shouldn't be here if there is no 2fold rotation, this should simply be p2...
                else:
                    if toprotsym == 2:
                        planegroup = "cmm"
                    else:
                        planegroup = "cm"
                        logging.warning("Unexpected point encountered in findSymmetry routine: FS001")
            if celltype == "rectangular" or efftype == "rectangular":
                # test mirror plane along both unit vectors at origin
                mirs = [False,False]
                for i in range(0,2):
                    spl = SymPlane(np.array([0,0]), abst[i], abst)
#                    tmpslab = copy.deepcopy(ts)
#                    tmpslab.mirror(spl)
#                    if tmpslab.isEquivalent(ts,eps):
                    if ts.isMirrorSymmetric(spl,eps):
                        mirs[i] = True
                if toprotsym == 2:
                    if mirs[0] and mirs[1]:
                        planegroup = "pmm"
                    else:
                        # test glide planes at origin
                        gldplane = DEFAULT
                        for i in range(0,2):
                            spl = SymPlane(np.array([0,0]), abst[i], abst)
#                            tmpslab = copy.deepcopy(ts)
#                            tmpslab.mirror(spl, glide=True)
#                            if tmpslab.isEquivalent(ts,eps):
                            if ts.isMirrorSymmetric(spl,eps, glide=True):
                                spl.type = "glide"
                                gldplane = spl
                        if gldplane != DEFAULT:
                            planegroup = "pmg"
                            self.orisymplane = gldplane
                        else:
                            # test glide plane at a/4:
                            spl = SymPlane(abst[0]/4, abst[1], abst)
#                            tmpslab = copy.deepcopy(ts)
#                            tmpslab.mirror(spl, glide=True)
#                            if tmpslab.isEquivalent(ts,eps):
                            if ts.isMirrorSymmetric(spl,eps, glide=True):
                                planegroup = "pgg"
                            else:
                                planegroup = "p2"
                    if planegroup in ["pmm","pmg","pgg"]: # each of these might be a mis-identified rcmm
#                        tmpslab = copy.deepcopy(ts)
#                        tmpslab.rotate((abst[0]+abst[1])/4, 2)  # check for 2fold rotation at diagonal/4
#                        if tmpslab.isEquivalent(ts,eps):
                        if ts.isRotationSymmetric((abst[0]+abst[1])/4, 2, eps):
                            planegroup = "rcmm"
                else:
                    logging.warning("Unexpected point encountered in findSymmetry routine: FS002")
            if celltype == "oblique" or efftype == "oblique":
                if toprotsym == 2:
                    planegroup = "p2"
                else:
                    planegroup = "p1"
                    logging.warning("Unexpected point encountered in findSymmetry routine: FS003")
        
        self.planegroup = planegroup
        if planegroup in ["pm", "pg", "cm", "rcm", "pmg"]:
            logging.info("Found plane group: "+planegroup+str(self.orisymplane.par))
        else:
            logging.info("Found plane group: "+planegroup)
        if planegroup == "rcm" and not rp.SYMMETRY_FIX:
            logging.warning("The given unit cell could be reduced to half the size in a centered representation. Consider reducing the unit cell size, or using SYMMETRY_FIX to set the symmetry to p1, pm or pg.")
        if planegroup == "rcmm" and not rp.SYMMETRY_FIX:
            logging.warning("The given unit cell could be reduced to half the size in a centered representation. Consider reducing the unit cell size, or using SYMMETRY_FIX to set the symmetry to p1, p2, pm, pg, pmm, pmg, or pgg.")

        # CHECK IF USER WANTS TO MANUALLY REDUCE THE SLAB SYMMETRY, AND WHETHER THE GIVEN REDUCTION IS LEGAL
#        if rp.SYMMETRY_FIX:
#            # DICTIONARY FOR ALLOWED SYMMETRY REDUCTIONS:
#            # TODO - if unit cell size should stay constant, kick out the relevant reductions
#            pgr = {'p1':[], 'p2':['p1'], 'pm':['p1'], 'pg':['p1'], 'cm':['p1','pm','pg'], 'pmm':['p1','p2','pm'], 
#                   'pmg':['p1','p2','pm','pg'],'pgg':['p1','p2','pg'],'cmm':['p1','p2','pm','pg','cm','pmm','pmg','pgg'],
#                   'p4':['p1','p2'],'p4m':['p1','p2','pm','cm','pmm','cmm','p4'],'p4g':['p1','p2','pg','cm','pgg','p4'],
#                   'p3':['p1'],'p3m1':['p1','p3','cm'],'p31m':['p1','p3','cm'],'p6':['p1','p2','p3'],
#                   'p6m':['p1','p2','cm','cmm','p3','p3m1','p31m','p6']}
#            if not '[' in rp.SYMMETRY_FIX:
#                targetsym = rp.SYMMETRY_FIX
#                if targetsym in ['pm','pg'] and not planegroup in ['cm','pmg']:
#                    logging.warning("Symmetry reduction from "+planegroup+" to "+targetsym+" requires a direction, which was not given. Input will be ignored, proceeding without symmetry reduction.")
#                    targetsym = planegroup
#            else:
#                rgx = re.compile(r'\s*(?P<group>[pmgcm]{2,3})\s*\[\s*(?P<i1>[-012]+)\s+(?P<i2>[-012]+)\s*\]')
#                m = rgx.match(rp.SYMMETRY_FIX)
#                targetsym = m.group('group')
#                tspar = [int(m.group('i1')), int(m.group('i2'))]
#            if targetsym != planegroup:
#                if not targetsym in pgr[planegroup]:
#                    logging.warning("The symmetry for the slab was defined as "+targetsym+" in the PARAMETERS file. This is not a valid symmetry reduction from the detected symmetry of "+planegroup+". Symmetry "+planegroup+" will be used.")
#                else:
#                    if targetsym not in ['pm','pg','cm','pmg']:
#                        self.planegroup = targetsym
#                    else:
#                        #### NOW NEED TO GO THROUGH ALL THE MORE TRICKY TRANSFORMS ####
#                        if planegroup == 'cm':  #reducing to pg or pm, no direction required    - TODO: this might get kicked out
#                            if targetsym == 'pg':   #shift origin to glide plane
#                                for at in self.atlist:
#                                    at.cartpos[0:2] -= abst[0]/2
#                                self.uCellShift -= abst[0]/2
#                                self.getFractionalCoordinates()
#                                self.orisymplane.type = "glide"     #since the origin shifts and the direction stays the same, nothing else needs to be changed about the symplane
#                            # double unit cell (both for pg and pm)
#                            self.collapseCartesianCoordinates()
#                            tmplist = self.atlist[:]
#                            for at in tmplist:
#                                tmpat = at.duplicate()
#                                for i in range(0,2):
#                                    if at.pos[i] >= 0.5:
#                                        at.cartpos -= abst[i]
#                                    else:
#                                        at.cartpos += abst[i]
#                            self.ucell = np.dot(np.array([1,1,0],[-1,1,0],[0,0,1]), self.ucell)
#                            self.uCellTransform = np.dot(np.array([1,1,0],[-1,1,0],[0,0,1]), self.uCellTransform)
#                            self.getFractionalCoordinates()
#                            self.planegroup = targetsym
#                        elif planegroup == 'pmm':
#                            if targetsym == 'pm':
#                                if (tspar[0],tspar[1]) in [(1,0),(0,1),(-1,0),(0,-1)]:  #no harm in allowing negative directions
#                                    self.planegroup = targetsym
#                                    self.orisymplane = SymPlane(np.array([0,0]),np.dot(tspar,abst),abst)
#                                else:
#                                    logging.warning("Invalid direction given for symmetry reduction from "+planegroup+" to "+targetsym+". Input will be ignored, proceeding without symmetry reduction.")
#                            else:
#                                self.planegroup = targetsym
#                        elif planegroup == 'pmg':
#                            if targetsym == 'pm': #needs origin shift
#                                sdir = np.dot(self.orisymplane.par,abst)
#                                for at in self.atlist:
#                                    at.cartpos[0:2] -= sdir/4
#                                self.uCellShift -= sdir/4
#                                self.getFractionalCoordinates()
#                                self.orisymplane.type = SymPlane(np.array([0,0]),np.array(self.orisymplane.dir[1],-self.orisymplane.dir[0]),abst)
#                                
#                            self.planegroup = targetsym
#                        elif planegroup == 'pgg':
#                            pass # TODO
#                        elif planegroup == 'cmm':
#                            pass #TODO
#                        elif planegroup == 'p4m':
#                            pass #TODO
#                        elif planegroup == 'p4g':
#                            pass #TODO
#                        elif planegroup == 'p3m1':
#                            pass #TODO
#                        elif planegroup == 'p31m':
#                            pass #TODO
#                        elif planegroup == 'p6m':
#                            pass #TODO
#                        else:
#                            logging.warning("Unexpected point encountered in findSymmetry routine: FS004")
#                if targetsym == self.planegroup:
#                    logging.info("The symmetry for the slab was reduced to "+targetsym+", as requested in the PARAMETERS file.")
        return planegroup
        
    def enforceSymmetry(self, rp, planegroup="fromslab"):
        """Finds how atoms are linked to each other based on the planegroup. If the planegroup argument is not given,
        the planegroup assigned to the slab will be used. Otherwise, the given planegroup has to be a subgroup of the 
        highest symmetry planegroup found for the slab."""
        if planegroup == "fromslab":
            if self.planegroup != "unknown":
                planegroup = self.planegroup
            else:
                logging.warning("Call to enforceSymmetry before findSymmetry; running findSymmetry now.")
                planegroup = self.findSymmetry(rp)
        eps = rp.SYMMETRY_EPS
        epsz = rp.SYMMETRY_EPS_Z
        abst = np.transpose(self.ucell[0:2,0:2])     #surface unit cell, transposed
        
        # FIND ATOM LINKING - HERE WORK WITH self INSTEAD OF ts, SINCE WE WANT TO ASSIGN PROPERTIES TO INDIVIDUAL ATOMS
        for at in self.atlist:  #first put all atoms in a list of their own
            at.linklist = [at]
        if not planegroup == "p1":  #p1 has no symmetry to check for
            self.createSublayers(epsz)
            self.sortOriginal()
            self.collapseCartesianCoordinates()
            # TEST ROTATION AT ORIGIN - TESTING ONLY HIGHEST ROTATIONAL ORDER IS ENOUGH
            if not planegroup in ["p1","pm","pg","cm"]:
                if planegroup in ["p2","pmm","pmg","pgg","cmm"]:
                    toprotsym = 2
                elif planegroup in ["p3","p3m1","p31m"]:
                    toprotsym = 3
                elif planegroup in ["p4","p4m","p4g"]:
                    toprotsym = 4
                elif planegroup in ["p6","p6m"]:
                    toprotsym = 6
                else:
                    logging.warning("Unexpected point encountered in enforceSymmetry routine: ES001")
                tmpslab = copy.deepcopy(self)
                tmpslab.rotate(np.array([0,0]),toprotsym)
                tmpslab.collapseCartesianCoordinates()
                angle = 2*np.pi/toprotsym
                m = np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
                for (sli,sl1) in enumerate(self.sublayers):
                    for (ati,at1) in enumerate(sl1.atlist):
                        for (atj,at2) in enumerate(tmpslab.sublayers[sli].atlist):
                            if not sl1.atlist[atj] in at1.linklist:     #don't check atoms that are already linked
                                if at1.isSameXY(at2.cartpos[0:2],eps):   #combine the two linklists
                                    at1.linklist.extend(sl1.atlist[atj].linklist)
                                    for at3 in sl1.atlist[atj].linklist:
                                        at3.symrefm = np.dot(m,np.dot(at1.symrefm,at3.symrefm))
                                        if not at3 == sl1.atlist[atj]: at3.linklist = at1.linklist
                                    sl1.atlist[atj].linklist = at1.linklist
                                    break
            # TEST MIRROR AND GLIDE PLANES
            if not planegroup in ["p2","p3","p4","p6"]:
                if planegroup in ["pm","pg","cm","pmg"]:
                    testplane = self.orisymplane                            # two possibilities, stored earlier
                elif planegroup in ["pmm","p4m","p6m"]:
                    testplane = SymPlane(np.array([0,0]), abst[0], abst)    # mirror at 0
                elif planegroup in ["pgg","p4g"]:
                    testplane = SymPlane(abst[1]/4, abst[0], abst, ty="glide")  # glide plane parallel to a at b/4
                elif planegroup in ["cmm","p31m"]:
                    testplane = SymPlane(np.array([0,0]), abst[0]+abst[1], abst)  # mirror a+b
                elif planegroup == "p3m1":
                    testplane = SymPlane(np.array([0,0]), abst[0]-abst[1], abst)  # mirror a-b
                else:
                    logging.warning("Unexpected point encountered in eindSymmetry routine: ES001")
                g = True if testplane.type == "glide" else False
                tmpslab = copy.deepcopy(self)
                tmpslab.mirror(testplane, glide=g)
                tmpslab.collapseCartesianCoordinates()
                angle = tl.angle(np.array([1,0]),testplane.dir)
                if testplane.dir[1] > 0: angle *= -1
                rotm = np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
                m = np.dot(np.linalg.inv(rotm),np.dot(np.array([[1,0],[0,-1]]),rotm))
                for (sli,sl1) in enumerate(self.sublayers):
                    for (ati,at1) in enumerate(sl1.atlist):
                        for (atj,at2) in enumerate(tmpslab.sublayers[sli].atlist):
                            if not sl1.atlist[atj] in at1.linklist:     #don't check atoms that are already linked
                                if at1.isSameXY(at2.cartpos[0:2],eps):   #combine the two linklists
                                    at1.linklist.extend(sl1.atlist[atj].linklist)
                                    for at3 in sl1.atlist[atj].linklist:
                                        at3.symrefm = np.dot(m,np.dot(at1.symrefm,at3.symrefm))
                                        if not at3 == sl1.atlist[atj]: at3.linklist = at1.linklist
                                    sl1.atlist[atj].linklist = at1.linklist
                                    break
        self.linklists = []     #re-create linklists
        for at in self.atlist:
            if len(at.linklist) > 1 and not at.linklist in self.linklists:  #don't keep the linklists of length 1
                self.linklists.append(at.linklist)
        
        # FIND ALLOWED MOVE DIRECTIONS FOR ATOMS
        lockpoints = []     #full list of rotation points (no duplication)
        ori = np.array([0,0])
        if planegroup in ["p2","pmm","pmg","pgg","cmm","p4","p4m","p4g","p6","p6m"]:
            lockpoints = [ori,np.array([0,0.5]),np.array([0.5,0]), np.array([0.5,0.5])]
            if planegroup in ["p6","p6m"]:
                lockpoints.extend([np.array([2/3,1/3]),np.array([1/3,2/3])])
        elif planegroup in ["p3","p3m1","p31m"]:
            lockpoints = [ori,np.array([2/3,1/3]),np.array([1/3,2/3])]
        lockplanes = []     #full list of mirror planes; no glide planes, since those don't restrict single atom movement
        if planegroup == "pm":   #include duplication of planes at ori+a / ori+b to avoid having to check atom positions +- a/b
            if np.array_equal(self.orisymplane.par, [1,0]):
                lockplanes = [self.orisymplane, SymPlane(abst[1]/2,abst[0],abst), SymPlane(abst[1],abst[0],abst,collapse=False)]  
            else:
                lockplanes = [self.orisymplane, SymPlane(abst[0]/2,abst[1],abst), SymPlane(abst[0],abst[1],abst,collapse=False)]
        if planegroup == "cm":
            if np.array_equal(self.orisymplane.par, [1,1]):
                lockplanes = [self.orisymplane, SymPlane((abst[0]-abst[1])/2,abst[0]+abst[1],abst,collapse=False), SymPlane((abst[1]-abst[0])/2,abst[0]+abst[1],abst,collapse=False)]
            else:
                lockplanes = [self.orisymplane, SymPlane((abst[0]+abst[1])/2,abst[0]-abst[1],abst), SymPlane((abst[0]+abst[1]),abst[0]-abst[1],abst,collapse=False)]
        if planegroup in ["pmm","p4m"]:
            for i in range(0,2):
                lockplanes.append(SymPlane(ori,abst[i],abst))
                lockplanes.append(SymPlane(abst[abs(i-1)]/2,abst[i],abst))
                lockplanes.append(SymPlane(abst[abs(i-1)],abst[i],abst,collapse=False))
        if planegroup in ["cmm","p4m"]:
            lockplanes.append(SymPlane(ori,abst[0]+abst[1],abst))
            lockplanes.append(SymPlane((abst[0]+abst[1])/2,abst[0]-abst[1],abst))
            if planegroup == "cmm":
                lockplanes.extend([SymPlane(ori,abst[0]-abst[1],abst),SymPlane(abst[0]+abst[1],abst[0]-abst[1],abst,collapse=False),SymPlane((abst[0]-abst[1])/2,abst[0]+abst[1],abst,collapse=False),SymPlane((abst[1]-abst[0])/2,abst[0]+abst[1],abst,collapse=False)])
        if planegroup == "pmg":
            if np.array_equal(self.orisymplane.par, [1,0]):
                lockplanes = [SymPlane(abst[0]/4,abst[1],abst),SymPlane(abst[0]*3/4,abst[1],abst)]
            else:
                lockplanes = [SymPlane(abst[1]/4,abst[0],abst),SymPlane(abst[1]*3/4,abst[0],abst)]
        if planegroup == "p4g":
            lockplanes = [SymPlane((abst[0]+abst[1])/4,abst[0]-abst[1],abst),SymPlane((abst[0]+abst[1])*3/4,abst[0]-abst[1],abst),SymPlane((abst[0]-abst[1])/4,abst[0]+abst[1],abst,collapse=False),SymPlane((abst[1]-abst[0])/4,abst[0]+abst[1],abst,collapse=False)]
        if planegroup in ["p3m1","p6m"]:
            lockplanes.extend([SymPlane((abst[0]+abst[1])/2,abst[0]-abst[1],abst),SymPlane(ori,abst[0]+2*abst[1],abst,index2=True),SymPlane(ori,2*abst[0]+abst[1],abst,index2=True),SymPlane(abst[0]+abst[1],abst[0]+2*abst[1],abst,index2=True,collapse=False),SymPlane(abst[0]+abst[1],2*abst[0]+abst[1],abst,index2=True,collapse=False)])
        if planegroup in ["p31m","p6m"]:
            lockplanes.extend([SymPlane(ori,abst[0],abst),SymPlane(ori,abst[1],abst),SymPlane(abst[0],abst[1],abst,collapse=False),SymPlane(abst[1],abst[0],abst,collapse=False),SymPlane(ori,abst[0]+abst[1],abst)])
        ts = copy.deepcopy(self)
        ts.projectCToZ()
        ts.collapseCartesianCoordinates()
        for at in ts.atlist:
            # first check points
            for p in lockpoints:
                if at.isSameXY(p,eps):
                    at.freedir = 0  #lock completely
                    break
            # then if not locked yet, check planes
            if not at.freedir == 0:
                for pl in lockplanes:
                    if tl.distanceLineThroughPointsFromPoint(pl.pos,pl.pos+pl.dir,at.cartpos[0:2]) < eps:
                        at.freedir = pl.par
                        break
        for at in self.atlist:
            for at2 in ts.atlist:
                if at2.oriN == at.oriN:
                    at.freedir = at2.freedir
                    break
        return
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    