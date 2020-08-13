# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

Contains generic functions used in the TensErLEED scripts.
"""

import logging
import numpy as np
import re
# import math
import scipy.spatial as sps
import itertools
import subprocess
import multiprocessing
import os
import shutil

logger = logging.getLogger("tleedm.base")

###############################################
#                 GLOBALS                     #
###############################################
periodic_table = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 
    'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 
    'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 
    'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 
    'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 
    'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 
    'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

elementCovalentRadii = {"H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96,
    "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05,
    "Cl": 1.02, "Ar": 1.06, "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60,
    "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24,
    "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20,
    "Br": 1.20, "Kr": 1.16, "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75,
    "Nb": 1.64, "Mo": 1.54, "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39,
    "Ag": 1.45, "Cd": 1.44, "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38,
    "I": 1.39, "Xe": 1.40, "Cs": 2.44, "Ba": 2.15, "La": 2.07, "Ce": 2.04,
    "Pr": 2.03, "Nd": 2.01, "Pm": 1.99, "Sm": 1.98, "Eu": 1.98, "Gd": 1.96,
    "Tb": 1.94, "Dy": 1.92, "Ho": 1.92, "Er": 1.89, "Tm": 1.90, "Yb": 1.87,
    "Lu": 1.87, "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44,
    "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "Tl": 1.45, "Pb": 1.46,
    "Bi": 1.48, "Po": 1.40, "At": 1.50, "Rn": 1.50, "Fr": 2.60, "Ra": 2.21,
    "Ac": 2.15, "Th": 2.06, "Pa": 2.00, "U": 1.96, "Np": 1.90, "Pu": 1.87,
    "Am": 1.80, "Cm": 1.69}
                        # from Cordero et al., 2008 (DOI: 10.1039/B801115J)

elementAtomicMass = {"H": 1.00797, "He": 4.00260, "Li": 6.941, "Be": 9.01218,
    "B": 10.81, "C": 12.011, "N": 14.0067, "O": 15.9994, "F": 18.998403, 
    "Ne": 20.179, "Na": 22.98977, "Mg": 24.305, "Al": 26.98154, "Si": 28.0855,
    "P": 30.97376, "S": 32.06, "Cl": 35.453, "K": 39.0983, "Ar": 39.948, 
    "Ca": 40.08, "Sc": 44.9559, "Ti": 47.90, "V": 50.9415, "Cr": 51.996, 
    "Mn": 54.9380, "Fe": 55.847, "Ni": 58.70, "Co": 58.9332, "Cu": 63.546, 
    "Zn": 65.38, "Ga": 69.72, "Ge": 72.59, "As": 74.9216, "Se": 78.96, 
    "Br": 79.904, "Kr": 83.80, "Rb": 85.4678, "Sr": 87.62, "Y": 88.9059, 
    "Zr": 91.22, "Nb": 92.9064, "Mo": 95.94, "Tc": 98, "Ru": 101.07, 
    "Rh": 102.9055, "Pd": 106.4, "Ag": 107.868, "Cd": 112.41, "In": 114.82, 
    "Sn": 118.69, "Sb": 121.75, "I": 126.9045, "Te": 127.60, "Xe": 131.30, 
    "Cs": 132.9054, "Ba": 137.33, "La": 138.9055, "Ce": 140.12, 
    "Pr": 140.9077, "Nd": 144.24, "Pm": 145, "Sm": 150.4, "Eu": 151.96, 
    "Gd": 157.25, "Tb": 158.9254, "Dy": 162.50, "Ho": 164.9304, "Er": 167.26, 
    "Tm": 168.9342, "Yb": 173.04, "Lu": 174.967, "Hf": 178.49, "Ta": 180.9479, 
    "W": 183.85, "Re": 186.207, "Os": 190.2, "Ir": 192.22, "Pt": 195.09, 
    "Au": 196.9665, "Hg": 200.59, "Tl": 204.37, "Pb": 207.2, "Bi": 208.9804, 
    "Po": 209, "At": 210, "Rn": 222, "Fr": 223, "Ra": 226.0254, "Ac": 227.0278,
    "Pa": 231.0359, "Th": 232.0381, "Np": 237.0482, "U": 238.029}

###############################################
#                 CLASSES                     #
###############################################

class BackwardsReader:
  """Simple class for reading a large file in reverse without having to read 
  the entire file to memory. 
  From http://code.activestate.com/recipes/120686/, adapted for python 3 and 
  changed somewhat to fit requirements.
  """
  def readline(self):
    while len(self.data) == 1 and ((self.blkcount * self.blksize) < self.size):
      self.blkcount = self.blkcount + 1
      line = self.data[0]
      try:
        self.f.seek(-self.blksize * self.blkcount, 2) # read from end of file
        self.data = ((self.f.read(self.blksize).decode(self.encoding)) 
                     + line).split("\n")
      except IOError:  # can't seek before the beginning of the file
        self.f.seek(0)
        self.data = (((self.f.read(self.size - (self.blksize 
                                                  *(self.blkcount-1))))
                      .decode(self.encoding) + line).split("\n"))

    if len(self.data) == 0:
      return ""

    line = self.data.pop()
    return line + '\n'

  def close(self):
    self.f.close()

  def __init__(self, file, blksize=4096, encoding="utf-8"):
    """initialize the internal structures"""
    # get the file size
    self.size = os.stat(file)[6]
    # how big of a block to read from the file...
    self.blksize = blksize
    # how many blocks we've read
    self.blkcount = 1
    # what encoding to use
    self.encoding = encoding
    self.f = open(file, 'rb')
    # if the file is smaller than the blocksize, read a block,
    # otherwise, read the whole thing...
    if self.size > self.blksize:
      self.f.seek(-self.blksize * self.blkcount, 2) # read from end of file
    self.data = (self.f.read(self.blksize)).decode(self.encoding).split("\n")
    # strip the last item if it's empty... a byproduct of the last line having
    # a newline at the end of it
    if not self.data[-1]:
      self.data.pop()

###############################################
#                FUNCTIONS                    #
###############################################

def getYfunc(ivfunc, v0i):
    """
    Returns the Y function for a given function. The derivative is 
    approximated by the slope between one point to the left and one to the 
    right, so the Y function does not have values for the highest and lowest 
    energies that are passed.
    Y = L/sqrt(1 + v0i^2 * L^2)
    L = ivfunc' / ivfunc

    Parameters
    ----------
    ivfunc : iterable
        Should be a 2D array, list of lists, or list of tuples containing 
        energies in ivfunc[:,0] and values in ivfunc[:,1], both as float.

    v0i : float
        Imaginary part of the inner potential

    Returns
    -------
    yfunc : numpy array
        2D array of the Y-function energies and values

    """
    if len(ivfunc) < 3:
        # cannot calculate derivative, return empty
        return np.array([[]])
    yfunc = None
    for i in range(1, len(ivfunc)-1):
        vals = ivfunc[i-1:i+2, 1]
        l = (vals[2]-vals[0])/(2*vals[1] + 1e-100) # +1e-100 to avoid div by 0
        y = l / (1 + (l*v0i)**2)
        try:
            yfunc = np.append(yfunc, [[ivfunc[i,0], y]], axis=0)
        except ValueError:
            if not yfunc:
                yfunc = np.array([[ivfunc[i,0], y]])
            else:
                raise
    return yfunc
    

def rotMatrix(order):
    """Returns a (2x2) matrix for in-plane rotation of the given rotation 
    order."""
    #these explicit definitions are likely useless, but sqrts might be
    #  marginally more accurate than sin/cos
    if order == 2:
        return np.array([[-1,0],[0,-1]])
    elif order == 3:
        return np.array([[-0.5,-np.sqrt(3)/2],[np.sqrt(3)/2,-0.5]])
    elif order == -3:
        return np.array([[-0.5,np.sqrt(3)/2],[-np.sqrt(3)/2,-0.5]])
    elif order == 4:
       return np.array([[0,1],[-1,0]])
    elif order == -4:
        return np.array([[0,-1],[1,0]])
    elif order == 6:
        return np.array([[0.5,np.sqrt(3)/2],[-np.sqrt(3)/2,0.5]])
    elif order == -6:
        return np.array([[0.5,-np.sqrt(3)/2],[np.sqrt(3)/2,0.5]])
    else:
        angle = 2*np.pi/order
        return np.array([[np.cos(angle),-np.sin(angle)],
                       [np.sin(angle),np.cos(angle)]])

def getTLEEDdir(home=""):
    """Finds directories in the 'source' folder that have names starting with
    'TensErLEED', then picks the one with the highest version number. Returns 
    a relative path to that directory, eg './source/TensErLEED-v1.6'."""
    sd = os.path.join(home,'source')
    l = [dn for dn in os.listdir(sd) if (os.path.isdir(os.path.join(sd,dn)) 
                                         and dn.startswith('TensErLEED'))]
    highest = 0.0
    founddir = ''
    for dn in l:
        try:
            f = float(dn.split('v')[-1])
            if f > highest:
                highest = f
                founddir = dn
        except:
            pass
    if founddir != '':
        return os.path.join(sd,dn)
    else:
        return ''

def getMaxTensorIndex(home=""):
    """Checks the Tensors folder for the highest Tensor index there, 
    returns that value, or zero if there is no Tensors folder or no valid 
    Tensors zip file."""
    if not os.path.isdir(os.path.join(home,".","Tensors")):
        return 0
    else:
        indlist = []
        rgx = re.compile(r'Tensors_[0-9]{3}\.zip')
        for f in [f for f in os.listdir(os.path.join(home,".","Tensors")) 
                  if (os.path.isfile(os.path.join(home,".","Tensors",f))
                      and rgx.match(f))]:
            m = rgx.match(f)
            if m.span()[1] == 15:  # exact match
                indlist.append(int(m.group(0)[-7:-4]))
        rgx = re.compile(r'Tensors_[0-9]{3}')
        for f in [f for f in os.listdir(os.path.join(home,".","Tensors")) 
                  if (os.path.isdir(os.path.join(home,".","Tensors",f))
                      and rgx.match(f))]:
            m = rgx.match(f)
            if m.span()[1] == 11:  # exact match
                indlist.append(int(m.group(0)[-3:]))
        if indlist:
            return max(indlist)
    return 0

def getTensors(index, required=True):
    """Fetches Tensor files from Tensors or archive with specified tensor 
    index. If required is set True, an error will be printed if no Delta files 
    are found."""
    dn = "Tensors_"+str(index).zfill(3)
    if not os.path.isdir(os.path.join(".","Tensors",dn)):
        if os.path.isfile(os.path.join(".","Tensors",dn+".zip")):
            try:
                logger.info("Unpacking {}.zip...".format(dn))
                os.mkdir(os.path.join(".","Tensors",dn))
                shutil.unpack_archive(os.path.join(".","Tensors",
                                                   dn+".zip"),
                                      os.path.join(".","Tensors",dn))
            except:
                logger.error("Failed to unpack {}.zip".format(dn))
                raise
        else:
            logger.error("Tensors not found")
            return ("Tensors not found")
    return 0

def getDeltas(index, required=True):
    """Fetches Delta files from Deltas or archive with specified tensor index. 
    If required is set True, an error will be printed if no Delta files are 
    found."""
    dn = "Deltas_"+str(index).zfill(3)
    if os.path.isdir(os.path.join(".","Deltas",dn)):
        for f in [f for f in os.listdir(os.path.join(".","Deltas",dn))
                  if (os.path.isfile(os.path.join(".","Deltas",dn,f)) 
                      and f.startswith("DEL_"))]:
            try:
                shutil.copy2(os.path.join(".","Deltas",dn,f), ".")
            except:
                logger.error("Could not copy existing delta files to "
                              "work directory")
                raise
    elif os.path.isfile(os.path.join(".","Deltas",dn+".zip")):
        try:
            logger.info("Unpacking {}.zip...".format(dn))
            shutil.unpack_archive(os.path.join(".","Deltas",dn+".zip"),
                                  ".")
        except:
            logger.error("Failed to unpack {}.zip".dn)
            raise
    elif required:
        logger.error("Deltas not found")
        return ("Deltas not found")
    return 0

def fortranContLine(s):
    """Takes a sting that might be too long to fit on a single line of fortran 
    code and splits it into continuation lines as necessary. Returns a string 
    with ampersands and line breaks."""
    limit = 72  #fortran line length limit
    if len(s) <= limit:
        return s
    o = s[:6]
    s = s[6:]
    while len(s) > (limit-6):   # 6 characters for beginning of line
        o += s[:(limit-6)]
        o += "&\n     &"
        s = s[(limit-6):]
    o += s
    return o

def fortranCompile(pre="", filename="", post="", 
                   logname = "fortran-compile.log"):
    """Assembles pre+filename+post to a filename, tries to execute via 
    subprocess.run, raises an exception if it fails."""
    fc = pre+" "+filename+" "+post
    fcl = linelist(fc)
    sep = ""
    if os.path.isfile(logname):
        sep = "\n\n"
    try:
        with open(logname, "a") as log:
            log.write(sep + "############\n# COMPILING: " + fc 
                      + "\n############\n\n")
        with open(logname, "a") as log:
            r = subprocess.run(fcl, stdout=log, stderr=log)
    except:
        logger.error("Error compiling "+filename)
        raise
    return r.returncode

def readIntRange(s):
    """Takes a string, returns a list of integers. If the string is a single 
    int or space-separated ints, the return value is a list containing only 
    those int. If the string contains i1-i2 or i1:i2, the list is 
    range(i1, i2+1), i.e. contains both i1 and i2."""
    l = []
    sl = s.split()
    for ss in sl:
        try:
            l.append(int(ss))
        except ValueError:
            if "-" in ss:
                spl = ss.split("-")
                try:
                    l.extend(list(range(int(spl[0]),int(spl[1])+1)))
                except:
                    return []
            elif ":" in ss:
                spl = ss.split(":")
                try:
                    l.extend(list(range(int(spl[0]),int(spl[1])+1)))
                except:
                    return []
    return list(set(l))

def readVector(s, ucell=None, defRelaltive=False):
    """Takes a string 'xyz[f1 f2 f3]', 'abc[f1 f2 f3]' or just '[f1 f2 f3]' 
    and returns the corresponding vector in cartesian coordinates, or None if 
    the string cannot be parsed."""
    m = re.match(r'\s*(xyz|abc)?\[\s*(?P<v1>[-0-9.]+)\s+'
                r'(?P<v2>[-0-9.]+)\s+(?P<v3>[-0-9.]+)\s*\]', s)
    if not m:
        return None
    try:
        v1 = float(m.group("v1"))
        v2 = float(m.group("v2"))
        v3 = float(m.group("v3"))
    except:
        return None
    if "abc" in s:
        if ucell is None:
            return None
        uct = ucell.transpose()
        vec = (v1*uct[0] + v2*uct[1] + v3*uct[2])
    else: #xyz
        vec = np.array([v1, v2, v3])
    #vec[2] *= -1  # to LEED coordinates
    return vec

def readIntLine(line, width=3):
    """
    Reads an (arbitrary length) line of integers with fixed width. Will try 
    to interpret everything as integers until the line ends.

    Parameters
    ----------
    line : str
        The line to interpret
    width : integer, optional
        The width of each integer. The default is 3.

    Returns
    -------
    List of integers

    """
    line = line.rstrip()
    l = []
    try:
        while len(line) > 0:
            l.append(int(line[:width]))
            line = line[width:]
    except:
        raise
    return l

def cosvec(x,y):
    """Returns the cosine of the angle between two vectors"""
    return np.dot(x,y)/(np.linalg.norm(x)*np.linalg.norm(y))

def dict_equal(d1,d2):
    """Returns true the dictionaries contain the same keys with the same 
    values"""
    if len({k: d1[k] for k in d1 if k in d2 and d1[k] == d2[k]})-len(d1) == 0:
        return True
    return False

def lcm(a,b):
    "Calculate the lowest common multiple of two integers a and b"
    return a*b//np.gcd(a,b)

def parseMathSqrt(s):
    try:
        f = float(s)
    except ValueError:
        f = 1.0
        if '*' in s:
            sl = s.split('*')
        else:
            sl = [s]
        sl2 = []
        for el in sl:
            if '/' in el:
                sl2.append(el.split('/')[0])
                for s in el.split('/')[1:]:
                    try:
                        sl2.append(1/parseMathSqrt(s))
                    except:
                        raise
            else:
                sl2.append(el)
        for el in sl2:
            try:
                f *= float(el)
            except ValueError:
                if not 'sqrt' in el:
                    logger.error('Could not interpret '+el+' as float value')
                    raise
                else:
                    p = re.compile(r'\s*sqrt\s*\(\s*(?P<value>[\d.]*)\s*\)')
                    m = p.match(el)
                    if m:
                        try:
                            f *= np.sqrt(float(m.group('value')))
                        except ValueError:
                            logger.error('Could not interpret '+
                                          m.group('value')+' as float value')
                            raise
                    else:
                        logger.error('Could not interpret '+el
                                      +' as float value')
                        raise
    return f

def angle(v1, v2, acute=True):
    """Returns the angle between two vectors"""
    angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1)*np.linalg.norm(v2)))
    if acute == True:
        return angle
    else:
        return 2 * np.pi - angle
    
def distanceLineThroughPointsFromPoint(p1,p2,r):
    """Gives the distance of point r from the line defined by p1 and p2 
    (in 2 or 3 dimensions)"""
    if len(p1) == 2:
        return (abs((p2[1]-p1[1])*r[0] - (p2[0]-p1[0])*r[1] 
                   + p2[0]*p1[1] - p2[1]*p1[0])
                / np.sqrt((p2[1]-p1[1])**2 + (p2[0]-p1[0])**2))
    elif len(p1) == 3:
        return np.linalg.norm(np.cross((r-p1),(p2-p1)))/np.linalg.norm(p2-p1)
    else:
        return False
        
def writeWoodsNotation(ucell):
    """Takes a unit cell (as a (2x2) matrix) and attempts to write it in Woods 
    Notation. Returns empty string if no Woods notation is found."""
    # !!! VERY INCOMPLETE, should at least detect simple c(a x b) cases
    # !!! Same functionality exists in guilib; replace at some point
    if ucell[1,0] == 0 and ucell[0,1] == 0:
        return("(" + str(int(ucell[0,0])) + "x" + str(int(ucell[1,1])) + ")")
    else:
        return ""
    
def readWoodsNotation(s, ucell):
    """Takes a string that should contain the transformation from the bulk to 
    the surface unit cell in Wood notation, as well as a bulk unit cell (from 
    which only the surface vectors are read). Returns a 2x2 transformation 
    matrix."""
    p = re.compile(r'\s*(?P<type>[PCpc]*)\s*\(\s*(?P<g1>.+)\s*[xX]\s*'
                   +r'(?P<g2>.+)\s*\)\s*[rR]*\s*(?P<alpha>[\d.]*)')
    # this regular expression matches if (any amount of whitespace at any 
    #                                                     point is ignored):
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
        logging.error('Could not read woods notation input '+s)
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
            logger.error('Could not read Woods notation angle: '
                          +m.group('alpha')+', setting angle to zero')
            alpha = 0.0
    alpha *= np.pi/180
    g1 = parseMathSqrt(m.group('g1'))
    g2 = parseMathSqrt(m.group('g2'))
    # get surface unit cell vectors from bulk unit cell (has surface 
    #  periodicity!!):
    if alpha == 0.0 and t == 'p':
        mat = np.array([[g1, 0.], [0., g2]], dtype=float)
    else:
        r = [ucell[:2,0],ucell[:2,1]]
        #q = np.linalg.norm(r[1])/np.linalg.norm(r[0])  
        # this would be to get from bulk vectors to surface, we have to reverse
        q = 1/(np.linalg.norm(r[1])/np.linalg.norm(r[0]))
        omega = angle(r[0],r[1])
                 #this is always constant in Wood notation, no need to reverse.
        if t == 'p':    
            #matrices from: Klaus Hermann; Crystallography and Surface 
            #                                Structure (Second Edition, Wiley)
            mat = ((1/np.sin(omega))
                  * np.array([[g1*np.sin(omega-alpha), g1*(1/q)*np.sin(alpha)], 
                             [-g2*q*np.sin(alpha), g2*np.sin(omega+alpha)]],
                             dtype=float))
        else:
            mat = ((1/(2*np.sin(omega)))
                  *np.array([[g1*np.sin(omega-alpha)-g2*q*np.sin(alpha), 
                              g1*(1/q)*np.sin(alpha)+g2*np.sin(omega+alpha)], 
                             [-g1*np.sin(omega-alpha)-g2*q*np.sin(alpha),
                              -g1*(1/q)*np.sin(alpha)+g2*np.sin(omega+alpha)]],
                            dtype=float))
    warn = False
    for i in range(0,2):
        for j in range(0,2):
            if abs(mat[i,j] - round(mat[i,j])) < 1e-4:
                mat[i,j] = round(mat[i,j])
            else:
                warn = True
    if warn:
        logger.warning("SUPERLATTICE matrix from Woods notation was "
                        "identified as:\n"+str(mat))
        logger.warning("SUPERLATTICE values do not round to "
                        "integer values. Check SUPERLATTICE parameter.")
    return mat
     
def linelist(line):
    """Splits a line at whitespace, deletes empty elements and line breaks, 
    then returns elements as a list"""
    llist1 = line.split()
    llist = []
    for part in llist1:
        if part != "":    # get rid of empty elements
            if part[-1][-1:] != '\n':   # get rid of line breaks
                llist.append(part)
            else:
                if part != '\n':
                    llist.append(part[:-1])
    return llist

def readToExc(llist):
    """For reading PARAMETERS files; takes a list, returns elements until the 
    first one that starts with an exclamation mark."""
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
    """Takes a list and a separator, splits strings in the list by the 
    separator, returns results as list of lists"""
    newlist = []
    sublist = []
    for el in llist:
        if not sep in el:
            sublist.append(el)
        else:
            pl = el.split(sep)
            if pl[0]: sublist.append(pl[0])
            newlist.append(sublist)
            sublist = []
            if len(pl) > 1:
                for i in range(1, len(pl)):
                    s = pl[i]
                    if s: sublist.append(s)
                    if not i==len(pl)-1:
                        newlist.append(sublist)
                        sublist = []
    newlist.append(sublist)
    return(newlist)
    
def splitMaxRight(s, sep):
    """Same as s.split(sep, maxsplit=1), but splitting at the first instance 
    from the right."""
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

def checkLattice(ab, eps = 0.001):
    """Takes unit vectors a,b as a 2x2 matrix, returns (lat,t), where lat is 
    a string "square", "rectangular", "hexagonal", "rhombic" or "oblique", and 
    t is a 2x2 transformation matrix which will transform the cell to obtuse 
    if it is rhombic or hexagonal."""
    # Author: Michele Riva; slightly modified for consistency 
    #   by Florian Kraushofer
    t = np.array([[1,0],[0,1]])
    c = cosvec(ab[0], ab[1])
    d = (np.linalg.norm(ab[0]) / np.linalg.norm(ab[1])) - 1
    if abs(c) < eps: #angle is 90°
        if np.abs(d)<eps:
            lat="square"
        else:
            lat="rectangular"
    elif np.abs(d)<eps: #rhombic or hex
        if c > eps: #angle is acute -> redefine to make it obtuse
            t = np.dot(np.array([[0,-1],[1,0]]), t) #this keeps the handedness
            ab = np.dot(t,ab)
        c = cosvec(ab[0], ab[1])
        if abs(c+1/2)<eps:
            lat = "hexagonal"
        else:
            lat = "rhombic"
    else:
        lat = "oblique"
    return (lat, t)

def reduceUnitCell(ab, eps = 0.001):
    """Takes an obtuse unit cell as a (2x2) matrix and reduces it to minimum 
    circumference, keeping the area constant. This might reduce oblique unit 
    cells to rectangular or hexagonal ones. Returns (ab, t, celltype), where 
    ab is the modified unit cell, t is the transformation matrix, and celltype 
    is a string describing the unit cell ("square", "rectangular", 
    "hexagonal", "rhombic" or "oblique")."""
    # Author: Michele Riva; slightly modified for consistency 
    #   by Florian Kraushofer
    (lat,t) = checkLattice(ab)
    ab = np.dot(t,ab)
    if lat == "oblique":
#     Transform lattice to have the shortest two vectors, with angle closest to
#       90°. This might bring it to rect, hex or rhombic.
#     If neither, will anyway transform to have the closest to rect.
#     ALGORITHM for reduction to closest to rect:
#     This is a discrete version of Gram-Schmidt's algorithm to find orthogonal
#        bases
#     At each iteration:
#       - order vectors by norm, the shortest first
#       - determine the projection of the second on the first, and calculate 
#           the nearest integer kk
#       - subtract from the second the projection calculated above
#       - check whether now the second is the smallest. If yes, repeat, 
#           otherwise finished.
        swap = np.array([[1,0],[0,1]]) 
        # matrix that keeps track of whether a and b are swapped
        while True: # Swap vectors if needed to get the shortest first
            if np.linalg.norm(ab[0])>np.linalg.norm(ab[1]):
                t0 = np.array([[0,1],[1,0]])
            else:
                t0 = np.array([[1,0],[0,1]])
            swap = np.dot(t0,swap)
            t = np.dot(t0,t)
            ab = np.dot(t0,ab)
            kk = int(np.round(np.dot(ab[0],ab[1])/np.dot(ab[0],ab[0])))
            t0 = np.array([[1,0],[-kk,1]])
            t = np.dot(t0,t)
            ab = np.dot(t0,ab)
            if np.linalg.norm(ab[0]) <= np.linalg.norm(ab[1]):
                break
        # Swap vectors back if they were overall swapped
        t = np.dot(swap,t)
        ab = np.dot(swap,ab)
#         END OF ALGORITHM. Now the lattice ab is closest to rectangular. It 
#           might be still any shape (square, rect, hex, rhombic, oblique)
        lat,t0=checkLattice(ab)
        t=np.dot(t0,t)
        ab=np.dot(t0,ab)
#       If ab is still oblique, try to see if it can be transformed to hex or
#         rhombic by choosing "a" not to be the shortest vector of all.
#       If possible, keep the new transformation. Otherwise, stick to the one 
#         that makes it closest to rectangular.
        if lat=="oblique": #lattice is still oblique
            # Re-swapping guarantees that that the matrix has on the first line
            #   the shortest possible vector, and on the second line the second
            #   shortest possible vector.
            # The only possible combinations that can lead to a rhombic/hex are
            #   a'=b+a or a'=b-a, depending on whether the angle is acute or 
            #   obtuse, respectively
            t2 = swap
            ab = np.dot(swap,ab)
            
            c = cosvec(ab[0],ab[1])
            t0 = [[-int(np.sign(c)),1],[0,1]]
            t2 = np.dot(t0,t2)
            ab = np.dot(t0,ab)
            
            lat, t0 = checkLattice(ab)
            t2 = np.dot(t0,t2)
            
            if lat=="oblique": 
                # lattice is still oblique, no transformation is needed (will 
                #   keep the one closest to rect)
                t2 = np.array([[1,0],[0,1]])
        else:
            t2 = np.array([[1,0],[0,1]])
        t = np.dot(t2,t)
    return ab, t, lat

def addUnequalPoints(l1,l2,eps,uniqueLists=False):
    """Adds all points from l1 to l2, if they are not already in l2 
    (+- epsilon)."""
    nl2 = l2[:]
    nl1 = l1[:]
    if len(l2) == 0:
        nl2 = nl1
    else:
        if not uniqueLists:
            #first get rid of duplicates in nl1
            tree = sps.KDTree(nl1)
            usepoint = [True]*len(nl1)
            for (i,p) in enumerate(nl1):
                if usepoint[i]:
                    for j in tree.query_ball_point(p, eps)[1:]: 
                        usepoint[j] = False
            nl1 = list(itertools.compress(nl1,usepoint))
        #then add remaining elements to l2:
        dl = sps.distance.cdist(np.vstack(tuple(l1)), np.vstack(tuple(l2)), 
                                'euclidean')
        for (i,sublist) in enumerate(dl):
            if min(sublist) >= eps:
                nl2.append(l1[i])
    return nl2

def pointIsInList(p,l,eps):
    """Checks whether a point is contained in a list of points, given an 
    epsilon."""
    if len(l) == 0: 
        return False
    dl = sps.distance.cdist(np.array([list(p)]), np.vstack(tuple(l)), 
                            'euclidean')
    if min(dl[0]) < eps:
        return True
    else:
        return False
    
def available_cpu_count():
    """ Number of available virtual or physical CPUs on this system, i.e.
    user/real as output by time(1) when called with an optimally scaling
    userspace-only program"""

    # cpuset
    # cpuset may restrict the number of *available* processors
    try:
        with open('/proc/self/status') as f:
            m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$', f.read())
        if m:
            res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
            if res > 0:
                return res
    except IOError:
        pass

    # Python 2.6+
    try:
        return multiprocessing.cpu_count()
    except NotImplementedError:
        pass

    # POSIX
    try:
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))

        if res > 0:
            return res
    except (AttributeError, ValueError):
        pass

    # Windows
    try:
        res = int(os.environ['NUMBER_OF_PROCESSORS'])
        if res > 0:
            return res
    except (KeyError, ValueError):
        pass

    # BSD
    try:
        sysctl = subprocess.Popen(['sysctl', '-n', 'hw.ncpu'],
                                  stdout=subprocess.PIPE)
        scStdout = sysctl.communicate()[0]
        res = int(scStdout)
        if res > 0:
            return res
    except (OSError, ValueError):
        pass

    # Linux
    try:
        with open('/proc/cpuinfo') as f:
            res = f.read().count('processor\t:')
        if res > 0:
            return res
    except IOError:
        pass

    # Solaris
    try:
        pseudoDevices = os.listdir('/devices/pseudo/')
        res = 0
        for pd in pseudoDevices:
            if re.match(r'^cpuid@[0-9]+$', pd):
                res += 1
        if res > 0:
            return res
    except OSError:
        pass

    # Other UNIXes (heuristic)
    try:
        try:
            with open('/var/run/dmesg.boot') as f:
                dmesg = f.read()
        except IOError:
            dmesgProcess = subprocess.Popen(['dmesg'], stdout=subprocess.PIPE)
            dmesg = dmesgProcess.communicate()[0]
        res = 0
        while '\ncpu' + str(res) + ':' in dmesg:
            res += 1
        if res > 0:
            return res
    except OSError:
        pass

    return -1