import numpy as np
from fractions import Fraction
from scipy import spatial as sp

#import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import cm #color maps
from matplotlib import colors as mpl_colors
from matplotlib.widgets import Button as mpl_button
from matplotlib.widgets import TextBox as mpl_text

#import timeit
import time

def cosvec(x,y):
    return np.dot(x,y)/(np.linalg.norm(x)*np.linalg.norm(y))

def CheckLattice(L):
    '''
    Input: L -- 2x2 lattice matrix to check

    Returns: (lat,T)
            - lat: string, ="sq","re","he","rh","ob" for square, rectangular, hexagonal, rhombic or oblique
            - T: transformation matrix that brings the input lattice matrix L into conventional representation when multiplied from the left.
                 The convention is: obtuse for rhombic and hexagonal, same as input for all others
    '''
    eps=1e-5 # Tolerance factor within which things are assumed to be equal
    T=E
    c=cosvec(L[0],L[1])
    d=np.linalg.norm(L[0])/np.linalg.norm(L[1])-1
    if abs(c)<eps: #angle is 90°
        if np.abs(d)<eps:
            print("Lattice is square")
            lat="sq"
        else:
            print("Lattice is rectangular")
            lat="re"
    elif np.abs(d)<eps: #rhombic or hex
        print("cosine=",c)
        if c>eps: #angle is acute -> redefine to make it obtuse
            T0=[[0,-1],[1,0]] #this keeps the handedness
        else:
            T0=E
        T=np.dot(T0,T)
        L=np.dot(T,L)
        c=cosvec(L[0],L[1])
        if abs(c+1/2)<eps:
            print("Lattice is hexagonal")
            lat="he"
        else:
            print("Lattice is rhombic")
            lat="rh"
    else:
        lat="ob"
    return (lat,T)
    


####################################################
###   THIS IS THE UNIT CELL REDUCTION FUNCTION   ###
####################################################
    
def GetHighSymLattice(L):
    '''
    Input: L -- 2x2 lattice matrix to transform to high symmetry
    
    Returns: T -- 2x2 transformation matrix that brings the input lattice to the highest possible symmetry representation when multiplied to L from the left
    '''

    lat,T=CheckLattice(L)
    L=np.dot(T,L)
    
    '''
    In what follows, T0 is used to define a specific elementary operation to be performed on the lattice.
    This is left-multiplied to T at each elementary step, so that T contains the overall transformation
    '''

    if lat=="ob": #oblique;
        print("Lattice is oblique, will try to see if it has higher symmetry")
        '''
        # Transform lattice to have the shortest two vectors, with angle closest to 90°.
        # This might bring it to rect, hex or rhombic.
        #
        # If neither, will anyway transform to have the closest to rect.
        '''
        
        '''
        # ALGORITHM for reduction to closest to rect:
        # This is a discrete version of Gram-Schmidt's algorithm to find orthogonal bases
        # At each iteration:
        #   - order vectors by norm, the shortest first
        #   - determine the projection of the second on the first, and calculate the nearest integer kk
        #   - subtract from the second the projection calculated above
        #    - check whether now the second is the smallest. If yes, repeat, otherwise finished.
        '''
        Swap=E #Matrix that keeps track of whether the first and second vectors are swapped at the end of this passage
        while True: # Swap vectors if needed to get the shortest first
            if np.linalg.norm(L[0])>np.linalg.norm(L[1]):
                T0=[[0,1],[1,0]]
            else:
                T0=E
            Swap=np.dot(T0,Swap)
            T=np.dot(T0,T)
            L=np.dot(T0,L)
            kk=int(np.round(np.dot(L[0],L[1])/np.dot(L[0],L[0])))
            T0=[[1,0],[-kk,1]]
            T=np.dot(T0,T)
            L=np.dot(T0,L)
            if np.linalg.norm(L[0])<=np.linalg.norm(L[1]):
                break
        # Swap vectors back if they were overall swapped
        T=np.dot(Swap,T)
        L=np.dot(Swap,L)
        
        # END OF ALGORITHM. Now the lattice L is closest to rectangular. It might be still any shape (square, rect, hex, rhombic, oblique)
        
        # Check which shape it has
        lat,T0=CheckLattice(L)
        T=np.dot(T0,T)
        L=np.dot(T0,L)
        
        # If it's still oblique, try to see if it can be transformed to hex or rhombic by choosing "a" not to be the shortest vector of all.
        #
        # If possible, keep the new transformation. Otherwise, stick to the one that makes it closest to rectangular
        #
        # All the operations that follow are stored in a matrix T2, to be later left-multiplied to T to get the full transformation
        #
        if lat=="ob": #lattice is still oblique, even if closest to rectangular
            print("After transformation the lattice is still oblique. Will try to reduce it to hexagonal/rhombic...")
            # Re-swapping guarantees that that the matrix has on the first line the shortest possible vector,
            # and on the second line the second shortest possible vector.
            #
            # The only possible combinations that can lead to a rhombic/hex are a'=b+a or a'=b-a,
            # depending on whether the angle is acute or obtuse, respectively
            #
            T2=Swap
            L=np.dot(Swap,L)
            
            c=cosvec(L[0],L[1])
            T0=[[-int(np.sign(c)),1],[0,1]]
            T2=np.dot(T0,T2)
            L=np.dot(T0,L)
            
            lat,T0=CheckLattice(L)
            T2=np.dot(T0,T2)
            
            if lat=="ob": #lattice is still oblique, no transformation is needed (will keep the one closest to rect)
                T2=E
                print("Lattice is oblique")
        else:
            T2=E
        T=np.dot(T2,T)
    else:
        T=E
    return T

def GetDomains(_group, M):
    '''
    Input:      _group: dictionary entry of PlaneGroup
             M:    2x2 matrix; matrix of superlattice. It MUST be already the one representing the highest 
                            symmetry lattice (will not be checked)
    
    Returns: _ops: list of the operations in _group that, applied to M, give all the distinct domains
    '''
    
    _ops=[]
    AllOps=[E,C2,C4,Cm4,C3,Cm3,C6,Cm6,Mx,My,M45,Mm45,M21,M12,M01,M10,M11,M1m1]
    
    if _group not in PlaneGroup.keys():
        raise ValueError("Group not identified")

    for op in PlaneGroup[_group]:
        if any(np.array_equal(np.dot(M,op),np.dot(Ei,np.dot(M,x))) for Ei in AllOps for x in _ops):
            pass
        else:
            _ops.append(op)
    
    print("%d domain(s), %d unique sub-pattern(s). Matrix/Matrices of the unique one(s):" %(len(PlaneGroup[_group]),len(_ops)))
    for op in _ops:
        print(np.dot(M,op))
        
    return _ops
    

# TODO - use for rotations
def RotateVtoX(v):
    '''
    Takes a generic 2D vector v as an input and outputs the transformation matrix T that rotates it to the x axis
    
    Input: v -- 1D array-like with 2 components (row vector)
    
    Returns: T -- 2x2 transformation matrix that rotates v to the x axis when multiplied to the RIGHT
    '''
    v=np.array(v)
    vnorm=np.linalg.norm(v)
    if len(v)>2:
        if abs(v[2])/vnorm>1e-5:
            raise ValueError("The vector to be rotated has a nonzero component along z")
        else:
            v=v[0:2]
    
    T=np.array([[v[0],-v[1]],[v[1],v[0]]])/vnorm
    
    return T

def GetMirrorToXY45(_group,L):
    '''
    Find the rotation matrix R that aligns mirror planes of lattice _L to x, y or +-45°
    
    Input:     _group, string; symmetry group of lattice
            _L, 2x2 matrix of lattice (i.e., Cartesian components of the lattice vectors
    
    Returns: R, 2x2 rotation matrix, to be right-multiplied to _L to get correct rotation
    '''
    if _group in ['p1','p2','pg[1 0]','pg[0 1]','pgg','p4','p3','p6']:
        R=E
    elif _group in ['pm[1 0]','pm[0 1]','cm[1 1]','cm[1-1]','rcm[1 0]','rcm[0 1]','pmg[1 0]','pmg[0 1]','pmm','cmm','rcmm']:
        # In this case, rotate by the smallest angle
        if _group in ['cm[1 1]','cm[1-1]','cmm']:
            vec=L[0]+L[1] #the mirror is along one diagonal
        else:
            vec=L[0] #the mirror is along one of the vectors
        R=RotateVtoX(vec) #first find matrix that rotates the first vector to x
        c=cosvec(vec,[1,0]) #find angle between unrotated first vector and x axis
        if abs(abs(c)-1/np.sqrt(2))>1e-5: #vector is closer to x axis than to y -> rotate to x axis
            pass
        else: #rotate vector to y axis instead
            R=np.dot(R,[[-1,0],[0,-1]])
    elif _group in ['p4m','p4g','p3m1','p31m','p6m']:
        #in this case it's sufficient to bring the first vector to the x axis
        R=RotateVtoX(L[0])
    else:
        raise ValueError("Group unrecognized")
    
    return R

def GenerateFullLattice(_Basis, _Type, _Lim):
    '''
    Generates list of lattice points given a _Basis.
    Input:      - _Basis, 2x2 array-like; contains the unit vectors as _Basis[0] and _Basis[1]
             - _Type, string; can be 'real' or 'reciprocal' to build real or reciprocal lattices.
                              This only affects the behavior of _Lim: real lattices are drawn up to a square,
                              Reciprocal ones within a circle
             - _Lim, scalar; determines which portion of the lattice is generated.
                             For _Type='real',  the lattice is generated up to a radius of 1.5*_Lim.
                                                One should then plot from -_Lim to +_Lim. This should 
                                                cover any post-rotation of the lattice that the user 
                                                might later request.
                             for _Type='reciprocal', the lattice is generated up to a radius of _Lim
    
    Returns: (Lat,hk)
             - Lat, 1D np.array of lattice points
             - hk, 1D np.array of indices. The lattice points are Lat=h*_Basis[0]+k*_Basis[1]
    '''
    
    if not np.array_equal(np.shape(_Basis),[2,2]):
#        raise ValueError("The input matrix for generating a lattice should be exactly 2x2. Input has shape ", shape(_Basis))
        raise ValueError()
    if _Type not in ['real','reciprocal']:
        raise ValueError("Unknown lattice type. Input must be 'real' or 'reciprocal'.")
    if hasattr(_Lim, '__len__'):
        raise ValueError("_Lim should be a scalar.")
    
    if _Type=='real':
        _Lim*=1.5
    
    # get limit for the loops that follow
    shortest=min(np.linalg.norm(_Basis[0]),np.linalg.norm(_Basis[1]),np.linalg.norm(_Basis[0]+_Basis[1])/2,np.linalg.norm(_Basis[0]-_Basis[1])/2)
    H=int(np.ceil(_Lim/shortest))
    
    # create grid of indices
    hk=np.array([(i,j) for i in range(-H,H+1) for j in range(-H,H+1)])
    '''
    Create lattice:
    Notice that, given a row vector of indices (h, k) and a basis in matrix form B=(a1,a2;b1,b2), the corresponding lattice
    point can be obtained as the row vector
        L=(L1,L2)=(h,k)B
    '''
    Lat=np.dot(hk,_Basis)
    
    # Now find all those lattice points that lie within _Lim, and use this as a mask for the output
    Mask = np.linalg.norm(Lat,axis=1)<=_Lim

    return (Lat[Mask],hk[Mask])

def GetPatternSymmetry(_group,hk):
    '''
    Starting from a list of integer beam indices and from the symmetry group of the real lattice,
    generates a list of integer grouping labels, where symmetry-equivalent beams have the same label
    Input:     _group, entry of dictionary PlaneGroup; symmetry group of the real lattice
            hk, array; list of integer beam indices
    Returns: Labels, np.array of integers, same length as hk. Equal labels correspond to symmetry equivalent beams
             CrossRefs, list of np.arrays, same length as hk, each np.array contains the positional indices of the beams 
                        equivalent to the current one (including self-reference)
    '''
    if _group not in PlaneGroup.keys():
        raise ValueError("Symmetry group unknown")
    if not all((hk.ravel())%1 == 0):
        raise ValueError("Some of the input indices are not integers")
    
    N=len(hk)
    Labels=np.arange(1,N+1)
    
    CrossRefs=[set([i]) for i in range(N)] # at first, use sets to avoid duplicates

    for op in PlaneGroup[_group]:
        hkT=np.dot(op,hk.transpose()).transpose().tolist() # transform hk according to the operation of the group
        for i in range(N):
            idx=hkT.index(hk[i].tolist()) # find which transformed hk is equal to the current beam
            CrossRefs[i]|=set([idx]) # add this to the set of cross references for the current beam
            Labels[idx]=Labels[i] # and make labels the same
            
    CrossRefs=np.array([np.array(list(x),dtype='int') for x in CrossRefs]) # then convert each set to a np.array of int, so that it's iterable and it's easy to process later
    
    if 'g' in _group: #group has glide
        if '[1 0]' in _group: #pg or pmg
            idx=np.where([x[1]==0 and abs(x[0]%2)==1 for x in hk])
        elif '[0 1]' in _group: #pg or pmg
            idx=np.where([x[0]==0 and abs(x[1]%2)==1 for x in hk])
        else:# pgg or p4g    
            idx=np.where([(x[1]==0 and abs(x[0]%2)==1) or (x[0]==0 and abs(x[1]%2)==1) for x in hk])
        Labels[idx]*=-1
    
    return (Labels,CrossRefs)

def FormatFractionalIndices(hk,M):
    '''
    Function that generates fractional indices of superstructure M whose integer indices are hk.
    Input:    - hk, 2d array of indices
            - M, superlattice matrix
    Returns:- names, 1D array of strings with formatted names for the indices 
    '''

    # The beam indices [hb,kb] for the surface indices [h,k] are
    #    [hb,kb]=np.dot([h,k],G),
    # with G=np.linalg.inv(M).transpose()
    hkbulk=np.dot(np.linalg.inv(M),hk.transpose()).transpose()
    
    #now get the formatting done
    mu=np.linalg.det(M).astype(int) #all fractional indices are of type nn/mu with nn integer
    
    names=[str(tuple([str(Fraction(hh).limit_denominator(abs(mu))) for hh in hkb])).replace('\'','') for hkb in hkbulk]
    
    return names

def mscatter(x,y,ax=None, m=None, **kw):
    '''
    This function provides an extension to matplotlib.pyplot.scatter allowing to define different marker styles for each point
    Input
    -----
     x: iterable, abscissas
     
     y: iterable, ordinates
     
     ax: axes to plot on
     
     m: iterable, markers, len(m)==len(x)
     
     **kw: keyword arguments, will be passed as they are to matplotlib.scatter
    
    Returns
    -------
    reference to the scatter plot, same as matplotlib.scatter
    '''
    import matplotlib.markers as mmarkers
    if not ax: ax=plt.gca()
    sc = ax.scatter(x,y,**kw)
    if (m is not None) and (len(m)==len(x)):
        paths = []
        for marker in m:
            if isinstance(marker, mmarkers.MarkerStyle):
                marker_obj = marker
            else:
                marker_obj = mmarkers.MarkerStyle(marker)
            path = marker_obj.get_path().transformed(marker_obj.get_transform())
            paths.append(path)
        sc.set_paths(paths)
    else:
        raise ValueError("Too many/too few markers")
    return sc

#Definition of matrix operations for groups, in 'fractional' coordinates:
#
#These two are good for all cells
E=np.array([[1,0],[0,1]])
C2=np.array([[-1,0],[0,-1]])
#
#These two are good for rectangular and square cells
Mx=np.array([[1, 0],[0,-1]])
My=np.array([[-1,0],[0, 1]])
#
#These are good for square cells
C4=np.array([[0,-1],[ 1,0]])
Cm4=np.array([[0, 1],[-1,0]])
M45=np.array([[0,-1],[-1,0]])
Mm45=np.array([[0, 1],[ 1,0]])
#
#These are good for rhombic and hex (both obtuse)
M1m1=np.array([[0,-1],[-1,0]])
M11=np.array([[0,1],[1, 0]])
M01=np.array([[-1,-1],[0,1]])
M10=np.array([[1,0],[-1,-1]])
#
#And these are good for hex only (obtuse)
C6=np.array([[1,1],[-1,0]])
Cm6=np.array([[0, -1],[1,1]])
C3=np.array([[0,1],[-1,-1]])
Cm3=np.array([[-1, -1],[1,0]])
M21=np.array([[1,1],[0,-1]])
M12=np.array([[-1,0],[1,1]])

#Definition of groups. Glides are replaced by mirrors.
PlaneGroup={
'p1':[E],
'p2':[E,C2],
'pm[1 0]':[E,Mx],
'pm[0 1]':[E,My],    
'pg[1 0]':[E,Mx],
'pg[0 1]':[E,My],
'cm[1 1]':[E,M11],
'cm[1-1]':[E,M1m1],
'rcm[1 0]':[E,Mx],
'rcm[0 1]':[E,My],
'pmm':[E,Mx,My,C2],
'pmg[1 0]':[E,Mx,My,C2],
'pmg[0 1]':[E,Mx,My,C2],
'pgg':[E,C2,Mx,My],
'cmm':[E,C2,M11,M1m1],
'rcmm':[E,C2,Mx,My],
'p4':[E,C2,C4,Cm4],
'p4m':[E,C2,C4,Cm4,Mx,My,M45,Mm45],
'p4g':[E,C2,C4,Cm4,Mx,My,M45,Mm45],
'p3':[E,C3,Cm3],
'p3m1':[E,C3,Cm3,M12,M21,M1m1],
'p31m':[E,C3,Cm3,M10,M11,M01],
'p6':[E,C6,C3,C2,Cm3,Cm6],
'p6m':[E,C6,C3,C2,Cm3,Cm6,M10,M11,M01,M12,M21,M1m1]
}

#Definition of lattice types.
LatType={
'sq':'Square',
're':'Rectangular',
'he':'Hexagonal',
'rh':'Rhombic',
'ob':'Oblique'
}

######################################################################################################

''' THIS PART IS JUST FOR MY CONVENIENCE DURING MANUAL TESTING - WILL NOT BE IN YOUR CODE '''
#Definitions of lattices
norm_a=3.905 #Angstrom

#oblique
alpha=50
a=norm_a*np.array([np.cos(alpha*np.pi/180),np.sin(alpha*np.pi/180)])
b=norm_a*np.array([0,1])
obl=np.array([a,b])

#square
a=norm_a*np.array([1,0])
b=norm_a*np.array([0,1])
sq=np.array([a,b])

#rectangular
a=2*norm_a*np.array([1,0])
b=np.sqrt(2)*norm_a*np.array([0,1])
rect=np.array([a,b])

#centered
alpha=110
a=norm_a*np.array([np.sin(alpha*np.pi/360),np.cos(alpha*np.pi/360)])
b=norm_a*np.array([-np.sin(alpha*np.pi/360),np.cos(alpha*np.pi/360)])
ctrd=np.array([a,b])

#hex
a=norm_a*np.array([1,0])
b=norm_a*np.array([np.cos(2*np.pi/3),np.sin(2*np.pi/3)])
hex=np.array([a,b])

bulk=hex#rect#
bulkGroup='p6m'#'p6m'#'p4m'#'pmm'#
surfGroup='pgg'#'p2'#'cm[1-1]'#'p6'#'pgg'#'p6'#'p6m'#'pm[1 0]'#

SUPERLATTICE=np.array([[3,2],[-2,0]])
surf=np.dot(SUPERLATTICE,bulk)

########################################################################################################
tstart=time.perf_counter()
''' YOUR NORMAL CODE WOULD START HERE '''

M=SUPERLATTICE            # This is the matrix imported from SUPERLATTICE
surf=np.dot(M,bulk)        # This is the 2x2 surface unit cell, as imported from the POSCAR

print("Modifying the SURFACE lattice for highest symmetry.......",end='')
Qsurf=GetHighSymLattice(surf)   # Qsurf now contains the matrix that brings surf to the highest symmetry
surf=np.dot(Qsurf,surf)
M=np.dot(Qsurf,M)
print("The SUPERLATTICE matrix that realizes the highest symmetry is\n",M)

'''
HERE ONE WOULD LOOK FOR THE SYMMETRY GROUP OF THE SURFACE WITH YOUR TOOL
OR USE THE ONE SUPPLIED BY THE USER (if compatible)
'''
surfGroup=surfGroup

# now get the bulk lattice from the surface one
bulk=np.dot(np.linalg.inv(M),surf)

print("Modifying the BULK lattice for highest symmetry.......",end='')
Qbulk=GetHighSymLattice(bulk) # Get the transformation matrix that brings also the bulk into the highest-symmetry representation (I guess in general the user should have done this already, but better safe than sorry)
bulk=np.dot(Qbulk,bulk) # transform bulk to high symmetry
M=np.dot(M,Qbulk)     # and also transform the matrix M so that in following iterations of the LEED program one gets the correct bulk lattice
                    # This will be saved in the OUTPARAMETERS file as a new SUPERLATTICE matrix
                    # Notice that it's multiplied to the right in this case, so that np.dot(np.linalg.inv(M),surf) gives the high-symmetry bulk now

'''
NOW YOUR CODE FOR DETECTING SYMMETRY RUNS ON THE BULK
It outputs the list of group symmetry operations into bulkGroup.
Notice that the list should be the one defined above, i.e., fractional coordinates, 
with glides replaced by mirrors, and discarding all translations
Probably one needs a function to generate the list of these operations from the group name? Or maybe just a dictionary is enough
TODO: add detection of screw axis perp. to surface for taking into account the presence of terraces rotated with respect to each other.
    It should apply a rotation that depends on the in-plane shape of the bulk unit cell (C2 for rectangular and rhombic; C4, C2 for square;
    C6,C3,C2 for hex; nothing for oblique) + a translation by a integer fraction of the out-of-plane bulk lattice constant (WHICH FRACTION?).
    DOES IT ALSO HAVE TO CHECK MIRRORS (i.e., glide planes)?
    Then fold coordinates back to the unit cell, and check if the untransformed and transformed are the same or not.
    It should also check whether the transformed one can be obtained from the untransformed one by an in-plane translation.
    In this latter case, the screw axis does not give additional domains in the LEED pattern.
'''
bulkGroup=bulkGroup

'''
NOW FIND A ROTATION MATRIX R TO BRING MIRRORS PARALLEL TO x/y/+-45°.
Notice that his applies to the LATTICEs (surface or bulk), while the matrix M stays the same
!!! WE SHOULD ALSO CHANGE THE BEAMINCIDENCE PARAMETER
    I DON'T THINK WE HAVE THIS IMPLEMENTED RIGHT NOW !!!
'''
print("Finding rotation that brings mirror planes to special directions....")
if surfGroup in ['p1','p2','pg[1 0]','pg[0 1]','pgg','p4','p3','p6']:
    '''
    These cases do not have mirror, so no preferred direction exists for 
    the orientation of the surface lattice. I think we can then chose the 
    orientation based on the bulk lattice in this case
    '''
    R=GetMirrorToXY45(bulkGroup,bulk)
else:
    R=GetMirrorToXY45(surfGroup,surf)

bulk=np.dot(bulk,R)
surf=np.dot(surf,R)
print(surf)

''' here the BEAMINCIDENCE rotation
if THETA != 0
    PHI=PHI-GetRotAngle(R)
'''

'''
Now build the unit vectors of the reciprocal lattice of the first surface domain
'''
as3D=np.array([*surf[0],0])    # 1st surface unit cell vector, including rotation. Size=(3x1)
bs3D=np.array([*surf[1],0])    # 2nd surface unit cell vector, including rotation. Size=(3x1)
z=np.array([0,0,1])            # z unit vector, useful for calculating the reciprocal space vectors
Energy=700                    # This is Emax from TEOENRANGE

Ss=np.vdot(as3D,np.cross(bs3D,z)) # Area of surface unit cell

asR3D=(2*np.pi/Ss)*np.cross(bs3D,z)
bsR3D=(2*np.pi/Ss)*np.cross(z,as3D)
surfR=np.array([asR3D[0:2],bsR3D[0:2]]) # This contains the reciprocal space unit vectors for the first domain

'''
And the unit vectors of the reciprocal lattice of the bulk
'''
bulkR=np.dot(M.transpose(),surfR)

'''
Now prepare the full reciprocal lattice of the first domain, and construct the integer beam indices (i.e., not those relative to the bulk)
'''
R=0.41*np.sqrt(Energy)    # Will limit the reciprocal space view to this radius, calculated from the typical geometry of LEED optics
LsurfR,hksurfR=GenerateFullLattice(surfR,'reciprocal',R)
LbulkR,dummy=GenerateFullLattice(bulkR,'reciprocal',R)
Nbeams=len(hksurfR)

'''
For representation, it's maybe good to also show the real space lattice of the bulk and of the first domain

Prepare the lattices here, so one can also plot them.
'''
FOV=4.2*max(surf.ravel()) # Half field of view. Make it such that there are at least 4 unit cells of the surface
Lsurf,dummy=GenerateFullLattice(surf,'real',FOV)
Lbulk,hkbulk=GenerateFullLattice(bulk,'real',FOV)

'''
NOW HANDLE DOMAINS
1)    Find the list DomOps of operations that generate the unique domains compatible with
    the symmetry of the bulk. This function will probably need to be modified a little after we find
    out how to handle screw axes (i.e., terraces related by rotations).
    PROBABLY ALSO NEEDS TO BE EDITED: unrotated square pg on square p4m should give two domains

    Notice that DomOps[0]=E always
'''
DomOps=GetDomains(bulkGroup,M)

Doms=[np.dot(np.dot(hksurfR,np.linalg.inv(np.dot(M,op)).transpose()),bulkR) for op in DomOps]
# Doms contains, for each domain, all the kx ky coordinates of the beams
NDoms=len(Doms)

print('bulkR=\n',bulkR)

'''    
2) Prepare the formatted names of the fractional indices
'''
Names=[]
for i in range(NDoms):
    Names.extend(FormatFractionalIndices(hksurfR,np.dot(M,DomOps[i])))
Names=np.array(Names)

'''
3) Prepare the full pattern of the surface
   LEEDx and LEEDy contain the horizontal and vertical coordinates of all the spots of the surface LEED pattern
'''
LEEDx=np.array([Doms[i][:,0] for i in range(NDoms)]).flatten() #horizontal component
LEEDy=np.array([Doms[i][:,1] for i in range(NDoms)]).flatten() #vertical component
LEED=np.array(list(zip(LEEDx,LEEDy)))

print(LEED)

t1=time.perf_counter()
print("Everything before plot:",t1-tstart,flush=True)

t0=time.perf_counter()
############# Prepare plots
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(13,6))

[s.get_xaxis().set_visible(False) for s in [ax1,ax2]]
[s.get_yaxis().set_visible(False) for s in [ax1,ax2]]

blat,dummy=CheckLattice(bulk)
slat,dummy=CheckLattice(surf)
fig.suptitle('Slab: %s, %s. Bulk: %s, %s.' %(LatType[slat],surfGroup,LatType[blat],bulkGroup), fontsize=12,y=0.07,zorder=-1)
                     # zorder decides which objects are background (smaller values) and which are foreground (larger values)

def EnableButton(butt,enable):
    if enable:
        butt.set_active(True)
        butt.color='0.85'
        butt.hovercolor='0.95'
        butt.label.set_color('0')
        [butt.ax.spines[i].set_color('0') for i in ['top','bottom','left','right']]
    else:
        butt.set_active(False)
        butt.color='0.5'
        butt.hovercolor='0.5'
        butt.label.set_color('0.6')
        [butt.ax.spines[i].set_color('0.6') for i in ['top','bottom','left','right']]
    fig.canvas.draw_idle()
    return

# add button for toggling visibility of domains
show_doms_text=['Show domains','Hide domains']
ax_button=plt.axes([0.863,0.08,0.09,0.045],zorder=-1)
ToggleDoms=mpl_button(ax_button,label=show_doms_text[1])
if NDoms==1: # hide it if there are no domains to toggle
    ax_button.set_visible(False)
Doms_no_text="%d inequivalent domain(s)"%(NDoms)
ax_button.text(x=0.5,y=1.2,s=Doms_no_text,ha='center',va='bottom')

# add textbox so the user can change the energy. Min = 10eV, Max=Eend of TEOENRANGE
ax_en=plt.axes([0.885,0.83,0.055,0.05],zorder=-1)
EnergyBox=mpl_text(ax_en, 'Energy (eV)', initial=str(Energy))
EnergyBox.label.set_ha('left')
EnergyBox.label.set_va('bottom')
EnergyBox.label.set_position((0,1.05))
En_Lims_text="Min: 10\nMax: %d"%(Energy)
ax_en.text(x=0,y=-0.1,s=En_Lims_text,ha='left',va='top',fontsize=8.5)

# also add two buttons next to it, so the user can step the energy by 10eV at a time
ax_enUp=plt.axes([0.942,0.857,0.011,0.023],zorder=-1)
EnUp=mpl_button(ax_enUp,label="▲")
EnUp.label.set_ha('center')
ax_enDown=plt.axes([0.942,0.83,0.011,0.023],zorder=-1)
EnDown=mpl_button(ax_enDown,label="▼")
EnDown.label.set_ha('center')

# add textbox for changing the rotation of the real and reciprocal space plots
# NB: the whole rotation control is on the left side of the figure since, if a button has
# zorder < axes zorder, it will not be accessible in the region where the two
# overlap. When placing the whole "rotation" block in between the real and 
# reciprocal space, the [0 1] "Vert." button is completely below the LEED pattern
# and cannot be operated.
ax_rot=plt.axes([0.03,0.83,0.065,0.05],zorder=-1)
RotBox=mpl_text(ax_rot, 'Rotation (°)', initial=str(0))
RotBox.label.set_ha('left')
RotBox.label.set_va('bottom')
RotBox.label.set_position((0,1.05))

# add two buttons to rotate stepwise, 10° at a time
ax_rotCW=plt.axes([0.097,0.857,0.011,0.023],zorder=-1)
RotCW=mpl_button(ax_rotCW,label="\u21bb")
RotCW.label.set_ha('center')
ax_rotCCW=plt.axes([0.097,0.83,0.011,0.023],zorder=-1)
RotCCW=mpl_button(ax_rotCCW,label="\u21ba")
RotCCW.label.set_ha('center')

# and a few more to place low-index directions horizontal/vertical
#[1 0]
ax_10tox=plt.axes([0.056,0.8,0.025,0.025],zorder=-1)
ax_10tox.text(x=-.05,y=0.5,s="[1 0]",ha='right',va='center')
OneZeroToX=mpl_button(ax_10tox,label="Hor.")
OneZeroToX.label.set_fontsize(8)

ax_10toy=plt.axes([0.083,0.8,0.025,0.025],zorder=-1)
OneZeroToY=mpl_button(ax_10toy,label="Vert.")
OneZeroToY.label.set_fontsize(8)

#[0 1]
ax_01tox=plt.axes([0.056,0.77,0.025,0.025],zorder=-1)
ax_01tox.text(x=-.05,y=0.5,s="[0 1]",ha='right',va='center')
ZeroOneToX=mpl_button(ax_01tox,label="Hor.")
ZeroOneToX.label.set_fontsize(8)

ax_01toy=plt.axes([0.083,0.77,0.025,0.025],zorder=-1)
ZeroOneToY=mpl_button(ax_01toy,label="Vert.")
ZeroOneToY.label.set_fontsize(8)

t1=time.perf_counter()
print("Prep figure:",t1-t0,flush=True)

def PlotRealSpace(Rotation=E):
    t0=time.perf_counter()
    ############# Plot real space lattice
    # Input: 2x2 Rotation matrix to apply to the whole lattice
    
    ax1.set_title("Real Space Lattice", fontsize=15,y=1.03)
    
    #Rotate the lattices
    RotBulk=np.dot(bulk,Rotation)
    RotSurf=np.dot(surf,Rotation) 
    RotSurfL=np.dot(Lsurf,Rotation)
    
    #Plot bulk as lines. A line is p+x*v, where p is any point along the line, and v a vector parallel to the line
    x=np.array([-FOV,FOV])
    K=max(hkbulk.ravel())
    for i in [0,1]:
        slope=[x*v for v in RotBulk[i]]
        for k in range(-K,K):
            ax1.plot(slope[0]+k*RotBulk[(i+1)%2,0],slope[1]+k*RotBulk[(i+1)%2,1],'gray',alpha=0.2)
    ax1.annotate("", xy=(RotBulk[0,0], RotBulk[0,1]), xytext=(0, 0), arrowprops=dict(arrowstyle="->",color='gray'))
    ax1.annotate("", xy=(RotBulk[1,0], RotBulk[1,1]), xytext=(0, 0), arrowprops=dict(arrowstyle="->",color='gray'))
    
    #surface as scatter plot
    ax1.scatter(RotSurfL[:,0],RotSurfL[:,1],s=8,c='k')
    ax1.set_xlim([-FOV,FOV])
    ax1.set_ylim([-FOV,FOV])
    ax1.annotate("", xy=(RotSurf[0,0], RotSurf[0,1]), xytext=(0, 0), arrowprops=dict(arrowstyle="->"))
    ax1.annotate("", xy=(RotSurf[1,0], RotSurf[1,1]), xytext=(0, 0), arrowprops=dict(arrowstyle="->"))
    
    
    t1=time.perf_counter()
    print("Draw real:",t1-t0,flush=True)    
    
    return

def _initLEEDPlot(Limit):
    '''
    Prepare graphical part of LEED pattern
    '''
    ax2.set_title("LEED Pattern", fontsize=15,zorder=-1,y=1.03)
    ax2.patch.set_alpha(0)
    ax2.axis('equal')
    [s.set_visible(False) for s in ax2.spines.values()]
    Screen=plt.Circle((0,0),Limit,color='k',fill=False)
    ax2.add_patch(Screen)
    axis_lim=[-Limit*1.02,Limit*1.02]
    ax2.set_xlim(axis_lim)
    ax2.set_ylim(axis_lim)
    
    return

def PlotLEED(NDomains=NDoms,LEEDEnergy=Energy,Rotation=E):
    '''
    This function processes and plots the LEED pattern of surface and bulk
    Input:      - ax: axes to plot on
             - NDomains: int, number of domains to plot. Can be either 1 or equal to NDoms. Used to toggle visibility of domains.
             - Limit: float, sets which portion of the LEED pattern to plot. All spots closer to (0,0) than Limit will be plotted.
             - Rotation: 2D np.array, Rotation.shape()==(2,2). Used to rotate the pattern.
    Returns: - NDomains: number of domains plotted
             - Limit: Radius of LEED pattern
             - Rotation: current pattern rotation matrix (with respect to standard convention)
    '''
    Limit=0.41*np.sqrt(LEEDEnergy)
    MarkScale=np.sqrt(Energy/LEEDEnergy) # this is the scaling factor for markers
    _initLEEDPlot(Limit)
    
    t0=time.perf_counter()
    '''
    Start from the pattern of the bulk, by selecting all those spots that are within Limit, then rotate the pattern
    '''
    Mask=np.linalg.norm(LbulkR,axis=1)<=Limit
    LEEDBulk=LbulkR[Mask]
    LEEDBulk=np.dot(LEEDBulk,Rotation)
    ax2.scatter(0,0,s=100*MarkScale,c='k', zorder=1e10) # Big dot at 0,0, drawn on top of everything
    ax2.scatter(LEEDBulk[:,0],LEEDBulk[:,1],s=25*MarkScale,facecolors='none', edgecolors='k') # Empty circles at integer-orders
    t1=time.perf_counter()
    print("Draw bulk LEED:",t1-t0,flush=True)
    
    t0=time.perf_counter()
    '''
    Now prepare stuff needed for LEED pattern of the slab
    '''
    Mask = np.full(len(LEED),True)    # Initialize mask to all true values
    
    # First limit the number of beams to the correct number of domains
    if NDomains==1:
        Mask[Nbeams:]=False
    Mask = np.logical_and(Mask,np.linalg.norm(LEED,axis=1)<=Limit) # And then include the radial masking
    
    # All quantities related to the LEED pattern need to be masked so that the indexing of beams is consistent throughout the rest
    LEEDSurf = LEED[Mask]
    LEEDSurf = np.dot(LEEDSurf,Rotation) # Rotate pattern
    LEEDNames=Names[Mask]
    
    Nbeams_tot = len(LEEDSurf) # total number of spots in the pattern
    Nbeams_dom = int(Nbeams_tot/NDomains) # number of beams in each domain
    
    '''
    Here would come the grouping of beams into symmetry equivalent ones, taking into account the 
    symmetry group of the slab. I suppose the output BeamGroupLabels to be an array with the same size as Lsurf,
    and containing 'grouping indices', while the output BeamGroups contains cross-referencing (at each index, 
    a list of all beams equivalent to the one at the specific index). For testing what follows, I implemented my 
    own version, that only looks at the [h,k] indices. It also gives negative indices to the extinct beams.
    '''
    BeamGroupLabels,BeamGroupsFirstDom=GetPatternSymmetry(surfGroup, hksurfR[Mask[0:Nbeams]]) # Do it at first on the first domain 
    
    '''
    And repeat the BeamGroupsFirstDom as many times as there are domains to have all the cross-references
    (beam indices for other domains are offset by Nbeams_dom)
    '''
    LEEDBeamGroups=[]
    for i in range(NDomains):
        LEEDBeamGroups.extend([bg+i*Nbeams_dom for bg in BeamGroupsFirstDom])
    LEEDBeamGroups=np.array(LEEDBeamGroups)
    
    # Group indices of all spots so that each group contains superimposing spots
    Tree=sp.cKDTree(LEEDSurf)
    LEEDSuperposed=np.array(Tree.query_ball_point(x=LEEDSurf,r=1e-8))
    
    
    # Get colors right for the domains
    if NDomains>1:
        colors=cm.gnuplot(np.linspace(0.1, 0.9, NDomains)) # initialize colors, one for each domain
        colorscale=[colors[i//Nbeams_dom] for i in range(Nbeams_tot)] # now initialize beam coloring by repeating 'colors' for each beam of each domain
        idx=np.where(np.greater([len(sp) for sp in LEEDSuperposed],1))[0] # find the indices of spots that are superimposed to others
        for i in idx: # make them gray
            colorscale[i]='gray'
    else:
        colorscale=['k']*Nbeams_dom
        
    # Get markers right: all spots will be 'o', except for those that are extinct due to
    # glide, that will be 'x' (a bit bigger than the others). The glide spots will also be semi-transparent.
    # Markers will also scale with energy with a scaling factor MarkScale
    mark=['o']*Nbeams_dom # create lists for the first domain
    size=[5*MarkScale]*Nbeams_dom
    colorscale=[mpl_colors.to_rgba(x, alpha=1) for x in colorscale]
    for beam in range(Nbeams_dom):
        if BeamGroupLabels[beam]<0:
            mark[beam]='x'
            size[beam]=10*MarkScale
            for k in range(NDomains):
                colorscale[beam+k*Nbeams_dom]=mpl_colors.to_rgba(colorscale[beam+k*Nbeams_dom],alpha=0.2)
    #and create as many copies as there are domains
    mark=np.array([mark]*NDomains).flatten()
    size=np.array([size]*NDomains).flatten()
    
    t2=time.perf_counter()
    
    # And finally plot the pattern
    Scatt=mscatter(LEEDSurf[:,0], LEEDSurf[:,1], ax=ax2, m=mark, s=size, c=colorscale)
    
    t3=time.perf_counter()
    print("Draw surf LEED:",t3-t2,flush=True)
    
    ####### Hovering functionality stuff
    # Create annotations: the largest number of annotations that can ever happen is probably 12
    N_annot_max=12

    annots = [ax2.annotate("", xy=(0,0), xytext=(15,15),bbox=dict(boxstyle="round", fc="w",alpha=0.8),textcoords="offset points", arrowprops=dict(arrowstyle="->"),fontsize=9) for an in range(N_annot_max)]
    [annot.set_visible(False) for annot in annots]
    
    t1=time.perf_counter()
    print("Draw surf LEED total:",t1-t0,flush=True)
    
    # Define the hovering functions
    def update_annot(ind, mousexy):
        '''
        This function updates the textboxes that label the LEED spots on which the mouse is hovering.
        Input:  - ind, dictionary. Contains details of the mouse movement event
                - mousexy, tuple. (x,y) position of mouse in data coordinates
        Returns: - number of annotations that have been updated
        '''

        '''
        ind['ind'] contains the indices within the scatter plot array of the points on which the mouse 
        is hovering. There might be more than one point if points overlap, or if the mouse 
        is hovering in between two points that are close together. 
        I'll keep only those points that are at the smallest distance from the mouse (within 1e-8)
        to avoid overcrowding the image.
        '''
        dist=np.linalg.norm(np.subtract(LEEDSurf[ind['ind']],mousexy),axis=1)
        closest= np.where(abs(dist-dist.min())<=1e-8)
        
        pts=np.array(ind['ind'])[closest]
        '''
        Now pts contains only the indices of the points closest to the mouse.
        In the following:
        - If there are overlapping points (i.e., same fractional indices), only one textbox
          with the fractional indices (only once) + which domains overlap
        - Otherwise only one textbox with the fractional indices and the domain to which the beam belongs
        In all cases there will be textboxs for each of the beams symmetry-equivalent to those selected.
        '''
        
        sym_eq=set([x for pt in pts for x in LEEDBeamGroups[pt]])
        '''
        sym_eq contains the list of indices of all beams under the mouse, and those that 
        are symmetry-equivalent to them. A set is used to avoid duplicates.
        
        Now one needs to include also the beams that overlap with some of those that are
        symmetry-equivalent to the ones selected, as well as their symmetry equivalents.
        '''
        for spot in sym_eq.copy():
            '''
            Find whether each element in sym_eq is superimposed to another LEED spot.
            If it is, add it to sym_eq. Also store all of the beams that are symmetry equivalent to it
            '''
            superposed=LEEDSuperposed[spot]
            for sup in superposed:
                if sup<Nbeams_tot:    # since LEEDSuperposed contains also spots from the other domains, if only one is shown
                    sym_eq |= set(LEEDBeamGroups[sup]) # add each superposed spot and its symmetry-equivalent
        '''
        Now that the list of spots is complete, group them into entries of uniques,
        such that each entry contains coincident beams
        '''        
        #get only the unique ones:
        uniques=[list(s) for s in set(tuple(u) for u in LEEDSuperposed[list(sym_eq)])] # the trick with the set is to remove duplicates
        
        N_annots=len(uniques)
        eps=1e-8 # tolerance factor to decide if spots are along the kx or ky axes
        for idx in range(N_annots):
            spots=uniques[idx]
            pos=Scatt.get_offsets()[spots[0]] # always take the first element in spots, since they anyway overlap
            annots[idx].xy=pos # position of arrow tip
            '''
            Prepare text for textbox:
            - Get the fractional index of the spot
            - If there's only one domain, print the index only
            - Otherwise, check if the spot superimposes with others. If yes, print which domains overlap,
              unless the spot is '(0, 0)'
            - Also in this case, the fractional index will appear only once
            TODO: domains that give glide-extinct spots should have a '(g)' appended to the domain number
            '''
            FractInd=LEEDNames[spots[0]]     # any of the overlapping spots will do, 
                                        # as they all have the same fractional index
            if NDomains==1 or FractInd=='(0 0)':
                DomStr=''            
            else:
                if len(spots)>1:
                    DomStr=', doms. {}'.format("+".join([str(sp//Nbeams_dom+1) for sp in spots]))
                else:
                    DomStr=', dom. %d'%(spots[0]//Nbeams_dom+1)
            text='%s%s'%(FractInd, DomStr)
            annots[idx].set_text(text) # text in the box
            # Align the textbox quadrant-wise
            if pos[0]>eps:
                annots[idx].set_ha('left')
            elif abs(pos[0])<eps:
                annots[idx].set_ha('center')
            else:
                annots[idx].set_ha('right')
            if pos[1]>eps:
                annots[idx].set_va('bottom')
            elif abs(pos[1])<eps:
                annots[idx].set_va('center')
            else:
                annots[idx].set_va('top')
            if np.linalg.norm(pos)<eps:
                annots[idx].set_va('bottom')
                annots[idx].set_ha('left')
            # And set its position such that the arrows are kinda radial
            arrowlength=20
            if np.linalg.norm(pos)>eps:
                annots[idx].xyann=arrowlength*pos/np.linalg.norm(pos)
            else:
                annots[idx].xyann=arrowlength/np.sqrt(2)*np.array([1,1])
        
        return N_annots
    
    def hover(event):
        if event.inaxes == ax2:
            cont, ind = Scatt.contains(event)
            if cont:
                mouse_pos=(event.xdata,event.ydata)
                N=update_annot(ind,mouse_pos)
                [annot.set_visible(True) for annot in annots[0:N]]
                [annot.set_visible(False) for annot in annots[N:]]
                fig.canvas.draw_idle()
            elif any([annot.get_visible for annot in annots]):
                [annot.set_visible(False) for annot in annots]
                fig.canvas.draw_idle()
        return

    fig.canvas.mpl_connect("motion_notify_event", hover)
    
    return (NDomains,LEEDEnergy,Rotation)

EnableButton(EnUp,False)
(DomsPlot,EnergyPlot,RotPlot)=PlotLEED()
PlotRealSpace(Rotation=RotPlot)


### Handle domains on/off event

def toggle_domains(event):
    global DomsPlot,EnergyPlot,RotPlot
    old_lab=ToggleDoms.label.get_text()
    idx=show_doms_text.index(old_lab)
    ToggleDoms.label.set_text(show_doms_text[(idx+1)%2])
    ax2.clear()
    if idx==0: #i.e., the old label was "Show Domains"
        (DomsPlot,EnergyPlot,RotPlot)=PlotLEED(NDomains=NDoms,LEEDEnergy=EnergyPlot,Rotation=RotPlot)
    else: #i.e., the old label was "Hide Domains"
        (DomsPlot,EnergyPlot,RotPlot)=PlotLEED(NDomains=1,LEEDEnergy=EnergyPlot,Rotation=RotPlot)
    fig.canvas.draw_idle()
    return

ToggleDoms.on_clicked(toggle_domains)

#### Handle events of Energy TextBox and up/down buttons

def EnergyChanged(text):
    global DomsPlot,EnergyPlot,RotPlot
    Enew=float(text)
    if Enew>=10 and Enew<=Energy:
        ax2.clear()
        (DomsPlot,EnergyPlot,RotPlot)=PlotLEED(NDomains=DomsPlot,LEEDEnergy=Enew,Rotation=RotPlot)
        fig.canvas.draw_idle()
    else:
        EnergyBox.set_val(str(EnergyPlot)) # put same value as before
    return

def EnergyUpDown(event):
    Scaling=1.2
    if event.inaxes==ax_enUp:
        if EnergyPlot*Scaling<=Energy:
            NewEn=EnergyPlot*Scaling
            EnableButton(EnUp,True)
            EnableButton(EnDown,True)
        elif EnergyPlot<Energy:
            NewEn=Energy
            EnableButton(EnUp,False)
        else:
            EnableButton(EnUp,False)
            return
    else:
        if EnergyPlot/Scaling>=10:
            NewEn=EnergyPlot/Scaling
            EnableButton(EnUp,True)
            EnableButton(EnDown,True)
        elif EnergyPlot>10:
            NewEn=10
            EnableButton(EnDown,False)
        else:
            EnableButton(EnDown,False)
            return

    EnergyChanged(NewEn)
    EnergyBox.set_val("{:.1f}".format(NewEn))
    return

EnergyBox.on_submit(EnergyChanged)
EnUp.on_clicked(EnergyUpDown)
EnDown.on_clicked(EnergyUpDown)

#### Handle rotation angle changes events
def RotationChanged(text):
    global DomsPlot,EnergyPlot,RotPlot
    
    ThetaNew=float(text)
    if ThetaNew>180:
        RotBox.set_val("{:.2f}".format(ThetaNew-360))
    elif ThetaNew<-180:
        RotBox.set_val("{:.2f}".format(ThetaNew+360))
    else:
        pass
    
    ThetaNew=np.radians(ThetaNew)
    NewRot=np.array([[np.cos(ThetaNew), np.sin(ThetaNew)],[-np.sin(ThetaNew),np.cos(ThetaNew)]])
    
    [ax.clear() for ax in [ax1,ax2]]
    (DomsPlot,EnergyPlot,RotPlot)=PlotLEED(NDomains=DomsPlot,LEEDEnergy=EnergyPlot,Rotation=NewRot)
    PlotRealSpace(Rotation=RotPlot)
    fig.canvas.draw_idle()
    return

def RotUpDown(event):
    Step=10
    OldTheta=np.degrees(np.arccos(RotPlot[0,0]))*np.sign(RotPlot[0,1])
    if event.inaxes==ax_rotCW:
        NewTheta=OldTheta-Step
    else:
        NewTheta=OldTheta+Step
    RotBox.set_val("{:.1f}".format(NewTheta))
    RotationChanged(NewTheta)
    return

def RotLowIndex(event):
    if event.inaxes in [ax_10tox,ax_10toy]:
        Rot=RotateVtoX(bulk[0])
        if event.inaxes==ax_10toy:
            DeltaTheta=90
        else:
            DeltaTheta=0
    else:
        Rot=RotateVtoX(bulk[1])
        if event.inaxes==ax_01toy:
            DeltaTheta=90
        else:
            DeltaTheta=0
    NewTheta=np.degrees(np.arccos(Rot[0,0]))*np.sign(Rot[0,1])+DeltaTheta
    RotBox.set_val("{:.1f}".format(NewTheta))
    RotationChanged(NewTheta)
    return

RotBox.on_submit(RotationChanged)
RotCW.on_clicked(RotUpDown)
RotCCW.on_clicked(RotUpDown)
OneZeroToX.on_clicked(RotLowIndex)
OneZeroToY.on_clicked(RotLowIndex)
ZeroOneToX.on_clicked(RotLowIndex)
ZeroOneToY.on_clicked(RotLowIndex)

tend=time.perf_counter()
print("total time",tend-tstart,flush=True)

plt.show()