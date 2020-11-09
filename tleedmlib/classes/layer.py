# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

Class storing position and atom list of a layer
"""

import numpy as np
# from tleedmlib import DEFAULT
        
class Layer:
    """To be used with Slab; has origin, atoms (a subset of the ones in slab), 
    and a number. Sublayers also use the "Layer" class."""
    def __init__(self,slab,num,isBulk=False,sublayer=False):
        self.slab = slab
        self.num = num      # consecutive layer numbering, 0 being highest
        self.isBulk = isBulk    # defined by BLAY in PARAMETERS file
        self.atlist = []    # atoms in this layer
        self.cartori = None     # origin: xy from POSCAR with possible 
                                #   displacements, z is highest atom
        self.cartbotz = None    # z position of lowest atom
        self.carttopz = None    # z position of highest atom
        if sublayer:  #only used by sublayers in the symmetry detection routine
            self.symposlist = []    #list of candidate positions for 
                                    #  rotation / mirror / glide planes
    
    def getLayerPos(self):
        """Gets a cartesian origin coordinate for the layer, using z of the 
        highest atom. x,y are calculated from the origin of an a,b unit cell at
        that height. Also assigns posInLayer for all atoms in this layer."""
        al = self.atlist[:]     #temporary copy
        al.sort(key=lambda atom: atom.pos[2])
        topat = al[-1]
        botat = al[0]
        oripos = np.array([0,0,topat.pos[2]])
        self.cartori = np.dot(self.slab.ucell, oripos)  
            #this gets x and y correct, but z still in the wrong direction and
            #  with origin as POSCAR
        self.cartori[2] = topat.cartpos[2]  
            #just take the z from the highest atom.
        self.cartbotz = botat.cartpos[2]
        self.carttopz = topat.cartpos[2]
        for atom in al:
            atom.posInLayer = atom.cartpos - self.cartori
    
