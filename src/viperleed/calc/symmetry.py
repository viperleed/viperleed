"""Functions for determining and setting slab symmetry."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-19'
__license__ = 'GPLv3+'

import copy
import itertools
import logging
import re

import numpy as np
import scipy.spatial as sps

from viperleed.calc.classes.atom_containers import AtomList
from viperleed.calc.classes.slab import AlreadyMinimalError
from viperleed.calc.classes.sym_entity import SymPlane
from viperleed.calc.files import parameters
from viperleed.calc.lib import leedbase
from viperleed.calc.lib.base import addUnequalPoints
from viperleed.calc.lib.base import dist_from_line
from viperleed.calc.lib.base import rotation_matrix
from viperleed.calc.lib.base import rotation_matrix_order
from viperleed.calc.lib.math_utils import angle

logger = logging.getLogger(__name__)


def getSymPosLists(sl, rp, pointlist, output=False):
    """Generates and returns a symposlist and hexsymposlist based on the
    pointlist, for example the list of cartesian in-plane atom positions
    in the lowest-occupied layer."""

    def uniqueSymPosList(sl, rp, spl, verbose, description=""):
        eps = rp.SYMMETRY_EPS
        abst = np.transpose(sl.ab_cell)
        tree = sps.KDTree(spl)
        if len(spl) > 1e6:
            logger.warning(
                "Approximate search for symmetry positions will "
                "be applied due to very large number of candidates - check "
                "result for errors! If the search fails, consider lowering "
                "the SYMMETRY_EPS z component or setting the "
                "SYMMETRY_FIND_ORI parameter to False.")
            rp.setHaltingLevel(2)
            approximate = min([np.linalg.norm(abst[0]),
                               np.linalg.norm(abst[1])])/100
            # significantly speeds up the search especially for
            #   large unit cells, but might give some errors.
        else:
            approximate = 0
        usepoint = [True]*len(spl)
        if verbose:
            t = int(len(spl)/10.0)
        for (i, p) in enumerate(spl):
            if verbose:
                if (i+1) % t == 0:
                    logger.debug(description+": "
                                 + str(int((i+1)/t)*10)+"%")
            if usepoint[i]:
                for j in tree.query_ball_point(p, eps,
                                               eps=approximate)[1:]:
                    usepoint[j] = False
        spl = list(itertools.compress(spl, usepoint))
        return(spl)

    symposlist = [np.array([0., 0.])]    # always check the origin
    hexsymposlist = []
    symposlist.extend(pointlist)
    symposlist.extend([(p1+p2)/2 for (p1, p2) in
                       itertools.combinations(pointlist, 2)])
    if sl.celltype == "hexagonal":
        hexsymposlist = [(p1+p2+p3)/3 for (p1, p2, p3) in
                         itertools.combinations(pointlist, 3)]
    # collapse to base unit cell:
    symposlist = list(np.dot(sl.ab_cell,
                             (np.dot(np.linalg.inv(sl.ab_cell),
                                     np.array(symposlist).transpose())
                              % 1.0)).T)
    if len(hexsymposlist) > 0:
        hexsymposlist = list(np.dot(sl.ab_cell,
                                    (np.dot(np.linalg.inv(sl.ab_cell),
                                            np.array(hexsymposlist).T)
                                     % 1.0)).T)
    # remove duplicates:
    verbose = (False if (len(symposlist)+len(hexsymposlist) < 1e5
                         or not output) else True)
    if verbose:
        logger.debug("Found {} candidates, removing duplicates..."
                     .format(len(symposlist)+len(hexsymposlist)))
    # symposlist:
    symposlist = uniqueSymPosList(sl, rp, symposlist, verbose,
                                  description="Atoms and pairs")
    # hexsymposlist:
    if len(hexsymposlist) > 1:
        hexsymposlist = uniqueSymPosList(sl, rp, hexsymposlist, verbose,
                                         description="Triplets")
    return(symposlist, hexsymposlist)


def findBulkSymmetry(sl, rp):
    """Checks the bulk slab for screw axes and glide planes."""
    eps = rp.SYMMETRY_EPS
    abt = sl.ab_cell.T.copy()
    rotsfound = []
    glidesfound = []
    ts = copy.deepcopy(sl)
    rp2 = copy.deepcopy(rp)
    # optimize C vector, removing vacuum below first. Notice that
    # this is the safer way, as it allows us to relax assumptions
    # on the z periodicity of sl
    ts.remove_vacuum_at_bottom(rp2)
    try:
        ts.ensure_minimal_c_vector(rp2)
    except AlreadyMinimalError:
        pass
    else:
        logger.debug('Bulk unit cell could be reduced with repeat vector '
                     '[{:.5f} {:.5f} {:.5f}]'.format(*(rp2.BULK_REPEAT)))

    # figure out what to check
    pcands = ts.get_candidate_layer_periods(eps)
    if len(pcands) == 0:
        return
    nl = ts.n_sublayers
    # check for screw axes
    checkrots = []
    if nl % 2 == 0:
        checkrots.extend([2, 4])
    if ts.celltype == "hexagonal" and (nl % 3 == 0):
        checkrots.extend([3, 6])
    for per in pcands:
        for ro in [ro for ro in checkrots if ro not in rotsfound]:
            if ts.is_bulk_screw_symmetric(ro, per, eps):
                rotsfound.append(ro)
    sl.bulk_screws = rotsfound
    if len(rotsfound) > 0:
        logger.debug("Bulk screw axes found: " +
                     ", ".join([str(v) for v in rotsfound]))
    # check for glide planes
    ori = np.array([0, 0])
    checkglides = [SymPlane(ori, abt[0], abt), SymPlane(ori, abt[1], abt),
                   SymPlane(ori, abt[0] + abt[1], abt),
                   SymPlane(ori, abt[0] - abt[1], abt)]
    if ts.celltype == "hexagonal":
        checkglides.extend([SymPlane(ori, 2*abt[0] + abt[1], abt),
                            SymPlane(ori, abt[0] + 2*abt[1], abt)])
    for per in pcands:
        for gl in [gl for gl in checkglides if gl not in glidesfound]:
            if ts.is_bulk_glide_symmetric(gl, per, eps):
                glidesfound.append(gl)
    sl.bulk_glides = glidesfound
    if len(rotsfound) > 0:
        logger.debug("Bulk glide planes found: " +
                     ", ".join([str(gl.par) for gl in glidesfound]))
    return


def findSymmetry(sl, rp, bulk=False, output=True, forceFindOri=False):
    """Reduces the unit cell if necessary and finds the plane group of the
    slab. Stores the plane group and the higher-symmetry direction of the
    unit cell, if there is one."""
    celltype = "ERROR - not recognized"
    planegroup = ""  # plane group will be stored in Hermann-Mauguin notation
    eps = rp.SYMMETRY_EPS
    # reduce surface unit cell
    abst = sl.ab_cell.T  # surface unit cell, transposed
#        usurf = np.array([[1,0],[0,1]])
    if rp.SYMMETRY_FIX != "p1":
        abst, usurf, celltype = leedbase.reduceUnitCell(abst)
    else:
        # if symmetry is switched off, don't try to change the cell.
        celltype, usurf = leedbase.checkLattice(abst)
    # usurf tracks unit cell changes
    # reduce bulk unit cell
    if not bulk:
        abbt = np.dot(np.linalg.inv(rp.SUPERLATTICE), abst)
        # bulk ab unit cell, transposed
        abbt, ubulk, _ = leedbase.reduceUnitCell(abbt)
        # ubulk tracks unit cell changes
    utr = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]])
    utr[:2, :2] = usurf
    if not np.array_equal(utr, np.identity(3, dtype=int)):
        if (np.array_equal(utr, np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]))
                and output):
            logger.info("The POSCAR unit cell was changed from an acute "
                        "to an obtuse form.")
            rp.checklist.append(
                "The POSCAR unit vector definitions have changed. Make sure "
                "to check beam labels for compatibility with new unit cell.")
        elif output:
            logger.warning("The POSCAR unit cell was not in its highest "
                           "symmetry form. The unit cell will be modified "
                           "with the transformation matrix: \n"+str(utr))
            rp.setHaltingLevel(1)
            rp.checklist.append(
                "The POSCAR unit vector definitions have changed. Make sure "
                "to check beam labels for compatibility with new unit cell.")
        # MODIFY SUPERLATTICE PARAMETER
        if not bulk and rp.superlattice_defined:
            rp.SUPERLATTICE = np.dot(usurf, np.dot(rp.SUPERLATTICE,
                                                   np.linalg.inv(ubulk)))
            parameters.modify(rp, "SUPERLATTICE")
        # MODIFY SYMMETRY_FIX PARAMETER
        if "[" in rp.SYMMETRY_FIX and not bulk:
            rgx = re.compile(r'\s*(?P<group>(pm|pg|cm|rcm|pmg))\s*\[\s*'
                             + r'(?P<i1>[-012]+)\s+(?P<i2>[-012]+)\s*\]')
            m = rgx.match(rp.SYMMETRY_FIX)
            targetsym = m.group("group")
            tspar = [int(m.group("i1")), int(m.group("i2"))]
            cartdir = np.dot(tspar, abst)
            newab = np.dot(sl.ab_cell, np.transpose(usurf))
            newdir = np.dot(np.linalg.inv(newab), cartdir)
            newdir = newdir / min(newdir)
            s = (targetsym+"[{:.0f} {:.0f}]".format(newdir[0], newdir[1]))
            parameters.modify(rp, "SYMMETRY_FIX", s)
        # MODIFY UNIT CELL
        sl.update_cartesian_from_fractional()
        sl.ucell_mod.append(('rmul', utr.T))
        sl.ucell = np.dot(sl.ucell, utr.T)
        # same as np.transpose(np.dot(utr,np.transpose(sl.ucell)))
        sl.collapse_cartesian_coordinates(update_origin=True)
        # gets fractional coordinates in the new unit cell and
        #   collapses appropriately

    # check cell type again
    abst = sl.ab_cell.T
    celltype, _ = leedbase.checkLattice(abst)
    sl.celltype = celltype
    if output:
        logger.info("Found unit cell type: "+celltype)
        logger.info("Starting symmetry search...")
    # FIND HIGHEST SYMMETRY ORIGIN
    sl.collapse_cartesian_coordinates()
    # create a testslab: C projected to Z
    ts = copy.deepcopy(sl)
    if bulk:        # check whether there are at least 2 atomic layers
        ts.create_sublayers(eps.z)
        if ts.n_sublayers < 2:
            ts = ts.with_double_thickness()
    ts.project_c_to_z()
    ts.sort_by_z()

    bigslab = copy.deepcopy(ts)
    # will have atoms duplicated and shifted to 4 unit cells
    #  ([0,0], [0,1], [1,0], [1,1])
    bigslab.atlist.strict = False
    for at in bigslab.atlist[:]:
        for i in range(0, 2):
            for j in range(0, 2):
                if not (i == 0 and j == 0):
                    tmpat = at.duplicate()
                    tmpat.pos[0] += i
                    tmpat.pos[1] += j
    bigslab.update_atom_numbers()
    bigslab.atlist.strict = True
    bigslab.update_cartesian_from_fractional(update_origin=True)
    # bigslab.full_update(rp)   can't do this - would collapse coordinates!
    bigslab.create_sublayers(eps.z)

    # find the lowest occupancy sublayer; comparing candidate
    #   axes / planes to this one will be fastest
    lowocclayer = bigslab.fewest_atoms_sublayer
    minlen = lowocclayer.n_atoms

    # find candidate positions for symmetry points / planes:
    if not rp.SYMMETRY_FIND_ORI and not forceFindOri and not bulk:
        symposlist = [np.array([0., 0.])]
        # only check origin, planes are defined explicitly for this case
        hexsymposlist = []
    else:
        if output:
            logger.debug("Generating candidate high-symmetry positions "
                         "from layer with "+str(int(minlen/4))+" atoms...")
        if minlen > 400:
            logger.warning(
                "The number of atoms in the smallest sublayer is very large. "
                "This can make the symmetry search take a very long time. To "
                "avoid this, either decrease the z component of the "
                "SYMMETRY_EPS parameter, or set the SYMMETRY_FIND_ORI "
                "parameter to False.")
        pl = [at.cartpos[0:2] for at in lowocclayer]
        symposlist, hexsymposlist = getSymPosLists(sl, rp, pl, output)

    comsymposlist = addUnequalPoints(symposlist, hexsymposlist, eps,
                                     uniqueLists=True)

    # we're done with the bigger slab, actually testing symmetry operations
    #   can be done just on the basic one.
    ts.create_sublayers(eps.z)
    lowocclayer = ts.sublayers[bigslab.sublayers.index(lowocclayer)]
    del bigslab

    # find potential rotation axes:
    if (rp.SYMMETRY_FIND_ORI or forceFindOri) and output:
        logger.debug("Checking for rotation axes: "
                     + str(len(comsymposlist))+" candidates...")
    # test potential rotation axes:
    toprotsym = 0   # keep track of the highest rotational symmetry so far

    toprotlist = []  # positions of highest rot-symmetry points
    for p in comsymposlist:
        rotsymorder = 0
        if ts.is_rotation_symmetric(p, 2, eps):
            rotsymorder = 2
            if celltype == "square":
                if ts.is_rotation_symmetric(p, 4, eps):
                    rotsymorder = 4
        if celltype == "hexagonal":
            if ts.is_rotation_symmetric(p, 3, eps):
                rotsymorder = 3
                if ts.is_rotation_symmetric(p, 6, eps):
                    rotsymorder = 6
        if rotsymorder > toprotsym:     # new best point found
            toprotsym = rotsymorder
            toprotlist = [p]
        elif rotsymorder != 0 and rotsymorder == toprotsym:
            toprotlist.append(p)
    if toprotsym > 0 and output:
        logger.debug("Highest rotation axis has order " + str(toprotsym))

    if toprotsym == 0 or toprotlist:
        # no rotation symmetry or multiple rotation centers with same order
        if output:
            logger.debug("Checking for mirror/glide planes...")

        # check for mirror/glides
        mirror = False
        glide = False
        symplanelist = []
        test_mirror_dirs = [(1, 0), (0, 1), (1, 1), (1, -1)]
        if celltype == "hexagonal":
            test_mirror_dirs.extend([(1, 2), (2, 1)])
        if not rp.SYMMETRY_FIND_ORI and not forceFindOri and not bulk:
            for (pa, pb) in [(0, 0), (0.25, 0.25), (0.25, -0.25)]:
                for (i, j) in test_mirror_dirs:
                    symplanelist.append(SymPlane(pa*abst[0]+pb*abst[1],
                                                 i*abst[0]+j*abst[1], abst,
                                                 index2=True))
        else:
            if toprotsym == 0:
                checklist = symposlist
            else:
                checklist = toprotlist
            for (k, pos) in enumerate(checklist):
                for (i, j) in test_mirror_dirs:
                    symplanelist.append(SymPlane(pos, i*abst[0]+j*abst[1],
                                                 abst, index2=True))
        if output:
            logger.debug(str(len(symplanelist))
                         + " candidates for mirror/glide planes found...")
        for spl in symplanelist:    # test the candidates
            if ts.is_mirror_symmetric(spl, eps):
                spl.type = "mirror"
                mirror = True
            elif ts.is_mirror_symmetric(spl, eps, glide=True):
                spl.type = "glide"
                glide = True
        droptypes = ["none"]
        if toprotsym > 0 and mirror:
            droptypes.append("glide")
        i = 0
        while i < len(symplanelist):
            if symplanelist[i].type in droptypes:
                symplanelist.pop(i)
            else:
                i += 1

    if toprotsym > 0:
        # identify origin
        if symplanelist:
            # if there are no symmetry planes, then keep the old list
            toprotlist = [spl.pos for spl in symplanelist]
        topsympoint = toprotlist[0]
        corners = (np.zeros(2), abst[0], abst[1], abst[0]+abst[1])
        mindist = min(np.linalg.norm(topsympoint-c) for c in corners)
        for p in toprotlist[1:]:
            if min(np.linalg.norm(p-c) for c in corners) < mindist:
                # for points with equal rotational symmetry, prioritize the
                #   one closest to the origin to avoid shifting the unit cell
                #   randomly every time
                mindist = min(np.linalg.norm(p-c) for c in corners)
                topsympoint = p
        # shift origin
        sl.translate_atoms_2d(-topsympoint)
        ts.translate_atoms_2d(-topsympoint)

    if toprotsym == 0: # should be an else
        # identify group:
        oriplane = None
        if not mirror:
            if not glide:
                planegroup = "p1"   # no origin shift required.
            else:
                planegroup = "pg"
                for spl in symplanelist:
                    if oriplane is None:
                        oriplane = spl
                    elif (spl.distanceFromOrigin(abst) <
                          oriplane.distanceFromOrigin(abst)):
                        # prioritize planes close to the origin of cell (1,1)
                        oriplane = spl
        else:
            if not glide:
                if celltype in ("square", "rectangular"):
                    planegroup = "pm"
                else:
                    planegroup = "p1"
            else:
                planegroup = "cm"
            for spl in symplanelist:
                if spl.type == "mirror":
                    if oriplane is None:
                        oriplane = spl
                    elif (spl.distanceFromOrigin(abst) <
                          oriplane.distanceFromOrigin(abst)):
                        oriplane = spl
                        # prioritize planes close to the origin of cell (1, 1)
            if planegroup == "cm":
                if celltype in ("square", "rectangular"):
                    # both mirrors and glides. if parallel to unit vectors: rcm
                    if tuple(oriplane.par) in [(1, 0), (0, 1)]:
                        planegroup = "rcm"
                elif (celltype == "hexagonal"
                      and tuple(oriplane.par) not in [(1, 1), (1, -1)]):
                    # Special case - requires rotation of unit cell to bring
                    # mirror/glide along either (11) or (1-1) direction.
                    abst, oriplane = mirror_to_diagonal(sl, rp, abst, oriplane)
        if oriplane is not None and oriplane.distanceFromOrigin(abst) > eps:
            # shift to closest point on oriplane
            shiftv = (np.array([oriplane.dir[1], -oriplane.dir[0]])
                      * oriplane.distanceFromOrigin(abst))
            if dist_from_line(
                    oriplane.pos, oriplane.pos+oriplane.dir, shiftv) > eps:
                shiftv = -1*shiftv
            # ts is not used any more in this case, otherwise those atoms
            #  would have to be shifted as well.
            sl.translate_atoms_2d(-shiftv)
        if oriplane:
            oriplane.pos = np.array([0, 0])
            sl.orisymplane = oriplane

    if not planegroup:
        # start by checking special case: in cmm, there are two inequivalent
        #  2fold axes, one of which would not have been found yet -> shift
        #  there (potentially), test
        if toprotsym == 2 and celltype in ["hexagonal", "rhombic"]:
            shiftslab = copy.deepcopy(ts)
            shiftslab.translate_atoms_2d(-abst[0]/2)
            # test diagonal mirror at shifted origin
            spl = SymPlane(np.array([0, 0]), (abst[0]+abst[1]), abst)
            if shiftslab.is_mirror_symmetric(spl, eps):
                planegroup = "cmm"
                ts = shiftslab
                sl.translate_atoms_2d(-abst[0]/2)  # correct origin

    if not planegroup:
        efftype = ""    # effective cell type
        if celltype == "hexagonal":
            if not (toprotsym == 3 or toprotsym == 6):
                efftype = "rhombic"
            else:
                # test mirror plane along unit vector at origin
                spl = SymPlane(np.array([0, 0]), abst[0], abst)
                if ts.is_mirror_symmetric(spl, eps):
                    if toprotsym == 6:
                        planegroup = "p6m"
                    else:
                        planegroup = "p31m"
                else:
                    if toprotsym == 6:
                        planegroup = "p6"
                    else:
                        # test mirror plane along both diagonals at origin
                        found = False
                        for i in [+1, -1]:
                            spl = SymPlane(np.array([0, 0]),
                                           (abst[0]+(i*abst[1])), abst)
                            if ts.is_mirror_symmetric(spl, eps):
                                found = True
                                break
                        if found:
                            planegroup = "p3m1"
                        else:
                            planegroup = "p3"
        elif celltype == "square":
            if toprotsym == 4:
                # test mirror plane along unit vector at origin
                spl = SymPlane(np.array([0, 0]), abst[0], abst)
                if ts.is_mirror_symmetric(spl, eps):
                    planegroup = "p4m"
                else:
                    # test glide plane along diagonal at origin
                    spl = SymPlane(np.array([0, 0]), (abst[0]+abst[1]),
                                   abst)
                    if ts.is_mirror_symmetric(spl, eps, glide=True):
                        planegroup = "p4g"
                    else:
                        planegroup = "p4"
            else:   # p1 p2 pm pg pmm pmg pgg
                # test mirror plane along both diagonals at origin
                found = False
                for i in [+1, -1]:
                    spl = SymPlane(np.array([0, 0]), (abst[0]+i*abst[1]),
                                   abst)
                    if ts.is_mirror_symmetric(spl, eps):
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
            for i in [+1, -1]:
                spl = SymPlane(np.array([0, 0]), (abst[0]+i*abst[1]), abst)
                if ts.is_mirror_symmetric(spl, eps):
                    found = True
                    break
            if not found:
                efftype = "oblique"
                # since we shouldn't be here if there is no 2fold rotation,
                #   this should simply be p2...
            else:
                if toprotsym == 2:
                    planegroup = "cmm"
                else:
                    planegroup = "cm"
                    logger.warning("Unexpected point encountered in "
                                   "findSymmetry routine: FS001")
                    rp.setHaltingLevel(2)
        if celltype == "rectangular" or efftype == "rectangular":
            # test mirror plane along both unit vectors at origin
            mirs = [False, False]
            for i in range(0, 2):
                spl = SymPlane(np.array([0, 0]), abst[i], abst)
                if ts.is_mirror_symmetric(spl, eps):
                    mirs[i] = True
            if toprotsym == 2:
                if mirs[0] and mirs[1]:
                    planegroup = "pmm"
                else:
                    # test glide planes at origin
                    gldplane = None
                    for i in range(0, 2):
                        spl = SymPlane(np.array([0, 0]), abst[i], abst)
                        if ts.is_mirror_symmetric(spl, eps, glide=True):
                            spl.type = "glide"
                            gldplane = spl
                    if gldplane is not None:
                        planegroup = "pmg"
                        sl.orisymplane = gldplane
                    else:
                        # test glide plane at a/4:
                        spl = SymPlane(abst[0]/4, abst[1], abst)
                        if ts.is_mirror_symmetric(spl, eps, glide=True):
                            planegroup = "pgg"
                        else:
                            planegroup = "p2"
                if planegroup in ["pmm", "pmg", "pgg"]:
                    # each of these might be a mis-identified rcmm
                    if ts.is_rotation_symmetric((abst[0]+abst[1])/4, 2, eps):
                        planegroup = "rcmm"
            else:
                logger.warning("Unexpected point encountered in "
                               "findSymmetry routine: FS002")
                rp.setHaltingLevel(2)
        if celltype == "oblique" or efftype == "oblique":
            if toprotsym == 2:
                planegroup = "p2"
            else:
                planegroup = "p1"
                logger.warning("Unexpected point encountered in "
                               "findSymmetry routine: FS003")
                rp.setHaltingLevel(2)

    sl.planegroup = planegroup
    sl.foundplanegroup = planegroup
    if planegroup in ["pm", "pg", "cm", "rcm", "pmg"]:
        sl.foundplanegroup = planegroup+str(sl.orisymplane.par)
    else:
        sl.foundplanegroup = planegroup
    if output:
        logger.info("Found plane group: "+sl.foundplanegroup)
    if planegroup == "rcm" and not rp.SYMMETRY_FIX and not bulk:
        logger.warning(
            "The given unit cell could be reduced to half the size in a "
            "centered representation. Consider reducing the unit cell size, "
            "or using SYMMETRY_FIX to set the symmetry to p1, pm or pg.")
        rp.setHaltingLevel(1)
    if planegroup == "rcmm" and not rp.SYMMETRY_FIX and not bulk:
        logger.warning(
            "The given unit cell could be reduced to half the size in a "
            "centered representation. Consider reducing the unit cell size, "
            "or using SYMMETRY_FIX to set the symmetry to p1, p2, pm, pg, "
            "pmm, pmg, or pgg.")
        rp.setHaltingLevel(1)

    # CHECK IF USER WANTS TO MANUALLY REDUCE THE SLAB SYMMETRY, AND
    #    WHETHER THE GIVEN REDUCTION IS LEGAL
    if rp.SYMMETRY_FIX and not bulk and sl.symbaseslab is None:
        planegroup = setSymmetry(sl, rp, rp.SYMMETRY_FIX)

    return planegroup

def mirror_to_diagonal(sl, rp, abst, oriplane):
    """Rotate cell to bring oriplane along a diagonal.

    The correct direction among (11) and (1-1) is the one
    forming a 60deg angle with the symmetry plane.
    """
    directions = (abst[0]+abst[1], abst[0]-abst[1])
    angles = [angle(d, oriplane.dir) for d in directions]
    dev_from_60 =  [(a+1e-5) % np.radians(60) < 1e-3
                                    for a in angles]
    if not any(dev_from_60) or all(dev_from_60):
        err = ("The POSCAR cell could not be "
                               "rotated to a higher symmetry")
        logger.error(err)
        raise RuntimeError(err)
    idx = 0 if dev_from_60[0] else 1
    m = rotation_matrix(angles[idx], dim=3)
    sl.update_cartesian_from_fractional()
    sl.ucell_mod.append(('lmul', m))
    sl.ucell = np.dot(m, sl.ucell)
    abst = sl.ab_cell.T
    direction = (abst[0]+abst[1], abst[0]-abst[1])[idx]
    sl.collapse_cartesian_coordinates(update_origin=True)
    oriplane = SymPlane(np.array([0, 0]), direction, abst)
    logger.warning("The POSCAR unit cell was changed to a higher "
                   "symmetry form. Make sure to check beam labels for "
                   "compatibility with new unit cell.")
    rp.checklist.append("Check modified unit cell")
    rp.setHaltingLevel(2)
    return abst,oriplane


def setSymmetry(sl, rp, targetsym):
    """Sets the symmetry of the slab, based on the one found by
    findSymmetry. Can be the planegroup from findSymmetry, or a reduction
    from it."""

    def invalidDirectionMessage(rp, planegroup, targetsym, setHaltingTo=2):
        logger.warning(
            "Invalid direction given for symmetry reduction "
            "from {} to {}. Input will be ignored, proceeding without "
            "symmetry reduction.".format(planegroup, targetsym))
        rp.setHaltingLevel(setHaltingTo)

    abst = sl.ab_cell.T  # surface unit cell, transposed
    # set high symmetry
    if "[" not in sl.foundplanegroup:
        sl.planegroup = sl.foundplanegroup
    else:
        ds = re.search(r'\[(.*)\]', sl.foundplanegroup).group(1)
        d = (int(ds.split()[0]), int(ds.split()[1]))
        sl.planegroup = sl.foundplanegroup.split("[")[0]
        sl.orisymplane = SymPlane(np.array([0, 0]),
                                  d[0]*abst[0]+d[1]*abst[1], abst)
        if sl.planegroup in ["pg", "pmg"]:
            sl.orisymplane.type = "glide"
    planegroup = sl.planegroup
    if targetsym in [planegroup, "found"]:
        return planegroup
    # DICTIONARY FOR ALLOWED SYMMETRY REDUCTIONS:
    pgr = {
        'p1': [], 'p2': ['p1'], 'pm': ['p1'], 'pg': ['p1'], 'cm': ['p1'],
        'rcm': ['p1', 'pm', 'pg'], 'pmm': ['p1', 'p2', 'pm'],
        'pmg': ['p1', 'p2', 'pm', 'pg'], 'pgg': ['p1', 'p2', 'pg'],
        'cmm': ['p1', 'p2', 'cm'],
        'rcmm': ['p1', 'p2', 'pm', 'pg', 'rcm', 'pmm', 'pmg', 'pgg'],
        'p4': ['p1', 'p2'],
        'p4m': ['p1', 'p2', 'pm', 'cm', 'pmm', 'cmm', 'p4'],
        'p4g': ['p1', 'p2', 'pg', 'cm', 'pgg', 'cmm', 'p4'],
        'p3': ['p1'], 'p3m1': ['p1', 'p3', 'cm'], 'p31m': ['p1', 'p3', 'cm'],
        'p6': ['p1', 'p2', 'p3'],
        'p6m': ['p1', 'p2', 'cm', 'cmm', 'p3', 'p3m1', 'p31m', 'p6']}
    if '[' not in targetsym:
        if ((targetsym in ['pm', 'pg', 'cm', 'rcm', 'pmg']
                and planegroup != 'pmg')
                or targetsym == 'cmm' and planegroup == 'p6m'):
            logger.warning(
                "Symmetry reduction from "+planegroup+" to "+targetsym
                + " requires a direction, which was not given. Input will be "
                "ignored, proceeding without symmetry reduction.")
            rp.setHaltingLevel(1)
            targetsym = planegroup
    else:
        rgx = re.compile(r'\s*(?P<group>(pm|pg|cm|rcm|pmg|cmm))\s*\[\s*'
                         + r'(?P<i1>[-012]+)\s+(?P<i2>[-012]+)\s*\]')
        m = rgx.match(targetsym)
        targetsym = m.group('group')
        tspar = [int(m.group('i1')), int(m.group('i2'))]
        # if unit cell was changed, adapt directions
        cellchange = False
        for op in sl.ucell_mod:
            if op[0] in ['lmul', 'rmul']:
                cellchange = True
                break
        if cellchange:
            logger.warning(
                "A symmetry change was requested relative to a unit cell "
                "vector, but the unit cell has been modified. Attempting to "
                "interpret direction in the old coordinate system...")
            rp.setHaltingLevel(1)
            tspar = np.dot(np.linalg.inv(sl.ab_cell),
                           np.dot(sl.ucell_ori[:2, :2], tspar))
            for i in range(0, 2):
                tspar[i] = round(tspar[i])
    if targetsym == planegroup:
        return planegroup
    if targetsym not in pgr[planegroup]:
        logger.warning(
            "Symmetry reduction to "+targetsym+" was "
            "requested. This is not a valid symmetry reduction from "
            "the detected symmetry of "+planegroup+". Symmetry "
            + planegroup+" will be used.")
        rp.setHaltingLevel(2)
    else:
        if targetsym not in ['pm', 'pg', 'cm', 'rcm', 'pmg', 'cmm']:
            sl.planegroup = targetsym
        else:
            # NOW NEED TO GO THROUGH ALL THE MORE TRICKY TRANSFORMS
            if planegroup == 'rcm':  # reducing to: pg, pm
                if targetsym == 'pg':  # shift origin to glide plane
                    shiftv = (0.25
                              * np.dot(np.array(sl.orisymplane.par[1],
                                                -sl.orisymplane.par[0]),
                                       abst))
                    sl.translate_atoms_2d(-shiftv)
                    sl.orisymplane.type = "glide"
                    # since the origin shifts and the direction stays the
                    #  same, nothing needs to be changed about the symplane
                sl.planegroup = targetsym
            elif planegroup == 'pmm':   # reducing to: pm
                if targetsym == 'pm':
                    if (tspar[0], tspar[1]) in [(1, 0), (0, 1),
                                                (-1, 0), (0, -1)]:
                        # no harm in allowing negative directions
                        sl.planegroup = targetsym
                        sl.orisymplane = SymPlane(np.array([0, 0]),
                                                  np.dot(tspar, abst), abst)
                    else:
                        invalidDirectionMessage(rp, planegroup, targetsym)
                else:
                    sl.planegroup = targetsym
                    logger.warning("Unexpected point encountered "
                                   "in setSymmetry routine: FS004")
                    rp.setHaltingLevel(1)
            elif planegroup == 'pmg':   # reducing to: pm, pg
                if targetsym == 'pm':  # needs origin shift
                    shiftv = 0.25*np.dot(sl.orisymplane.par, abst)
                    sl.translate_atoms_2d(-shiftv)
                    sl.orisymplane = SymPlane(
                        np.array([0, 0]), np.array(sl.orisymplane.dir[1],
                                                   -sl.orisymplane.dir[0]),
                        abst)
                sl.planegroup = targetsym
            elif planegroup == 'pgg':   # reducing to: pg
                if (tspar[0], tspar[1]) in [(1, 0), (0, 1), (-1, 0), (0, -1)]:
                    shiftv = 0.25*np.dot(np.array(tspar[1], -tspar[0]), abst)
                    sl.translate_atoms_2d(-shiftv)
                    sl.orisymplane = SymPlane(np.array([0, 0]),
                                              np.dot(tspar, abst), abst)
                    sl.orisymplane.type = "glide"
                    sl.planegroup = targetsym
                else:
                    invalidDirectionMessage(rp, planegroup, targetsym)
            elif planegroup == 'cmm':   # reducing to: cm
                if (tspar[0], tspar[1]) in [(1, 1), (1, -1),
                                            (-1, -1), (-1, 1)]:
                    sl.orisymplane = SymPlane(np.array([0, 0]),
                                              np.dot(tspar, abst), abst)
                    sl.planegroup = targetsym
                else:
                    invalidDirectionMessage(rp, planegroup, targetsym)
            elif planegroup == 'rcmm':
                # reducing to: pm, pg, rcm, pmg
                if (tspar[0], tspar[1]) in [(1, 0), (0, 1), (-1, 0), (0, -1)]:
                    if targetsym in ['pm', 'pg', 'rcm']:
                        sl.orisymplane = SymPlane(np.array([0, 0]),
                                                  np.dot(tspar, abst), abst)
                        if targetsym == 'pg':
                            # shift origin to glide plane
                            shiftv = 0.25*np.dot(np.array(tspar[1],
                                                          -tspar[0]), abst)
                            sl.translate_atoms_2d(-shiftv)
                            sl.orisymplane.type = "glide"
                        sl.planegroup = targetsym
                    elif targetsym == "pmg":
                        shiftv = 0.25*(abst[0]+abst[1])
                        sl.translate_atoms_2d(-shiftv)
                        sl.orisymplane = SymPlane(np.array([0, 0]),
                                                  np.dot(tspar, abst), abst)
                        sl.orisymplane.type = "glide"
                    sl.planegroup = targetsym
                else:
                    invalidDirectionMessage(rp, planegroup, targetsym)
            elif planegroup == 'p4m':   # reducing to: pm, cm, cmm
                if targetsym == 'cmm':
                    sl.planegroup = targetsym
                elif targetsym == 'pm':
                    if (tspar[0], tspar[1]) in [(1, 0), (0, 1),
                                                (-1, 0), (0, -1)]:
                        sl.orisymplane = SymPlane(np.array([0, 0]),
                                                  np.dot(tspar, abst), abst)
                        sl.planegroup = targetsym
                    else:
                        invalidDirectionMessage(rp, planegroup, targetsym)
                elif targetsym == 'cm':
                    if (tspar[0], tspar[1]) in [(1, 1), (1, -1),
                                                (-1, -1), (-1, 1)]:
                        sl.orisymplane = SymPlane(np.array([0, 0]),
                                                  np.dot(tspar, abst), abst)
                        sl.planegroup = targetsym
                    else:
                        invalidDirectionMessage(rp, planegroup, targetsym)
            elif planegroup == 'p4g':   # reducing to: pg, cm, cmm
                allowed = True
                if targetsym == 'cmm':
                    shiftv = 0.5*abst[0]
                elif targetsym == 'pg':
                    if (tspar[0], tspar[1]) in [(1, 0), (0, 1),
                                                (-1, 0), (0, -1)]:
                        shiftv = 0.25*np.dot(np.array(tspar[1], -tspar[0]),
                                             abst)
                        sl.orisymplane = SymPlane(np.array([0, 0]),
                                                  np.dot(tspar, abst), abst)
                        sl.orisymplane.type = 'glide'
                    else:
                        allowed = False
                elif targetsym == 'cm':
                    if (tspar[0], tspar[1]) in [(1, 1), (1, -1),
                                                (-1, -1), (-1, 1)]:
                        shiftv = 0.5*abst[0]
                        sl.orisymplane = SymPlane(np.array([0, 0]),
                                                  np.dot(tspar, abst), abst)
                    else:
                        allowed = False
                if allowed:
                    sl.translate_atoms_2d(-shiftv)
                    sl.planegroup = targetsym
                else:
                    invalidDirectionMessage(rp, planegroup, targetsym)
            elif planegroup == 'p3m1':  # reducing to: cm
                if (tspar[0], tspar[1]) in [(1, -1), (-1, 1), (1, 2),
                                            (-1, -2), (2, 1), (-2, -1)]:
                    chir = 1
                    if (abs(np.dot(abst[1],
                                   np.dot(np.array([[-0.5, -np.sqrt(3)/2],
                                                    [np.sqrt(3)/2, -0.5]]),
                                          abst[0]))) > 0.01):
                        chir = -1   # left-handed unit cell -> invert rotations
                    if (tspar[0], tspar[1]) in [(1, 2), (-1, -2)]:
                        sl.rotate_unit_cell(6*chir)  # rotate 60° clockwise
                    elif (tspar[0], tspar[1]) in [(2, 1), (-2, -1)]:
                        sl.rotate_unit_cell(-6*chir)  # rotate 60° countercl.
                    abst = sl.ab_cell.T
                    sl.orisymplane = SymPlane(np.array([0, 0]),
                                              abst[0]-abst[1], abst)
                    sl.planegroup = targetsym
                else:
                    invalidDirectionMessage(rp, planegroup, targetsym)
            elif planegroup == 'p31m':  # reducing to: cm
                if (tspar[0], tspar[1]) in [(1, 0), (-1, 0), (0, 1), (0, -1),
                                            (1, 1), (-1, -1)]:
                    chir = 1
                    if abs(np.dot(abst[1], np.dot(np.array(
                            [[-0.5, -np.sqrt(3)/2], [np.sqrt(3)/2, -0.5]]),
                            abst[0]))) > 0.01:
                        chir = -1   # left-handed unit cell
                    if (tspar[0], tspar[1]) in [(1, 0), (-1, 0)]:
                        sl.rotate_unit_cell(6*chir)  # rotate 60° clockwise
                    elif (tspar[0], tspar[1]) in [(0, 1), (0, -1)]:
                        sl.rotate_unit_cell(-6*chir)  # rotate 60° countercl.
                    abst = sl.ab_cell.T
                    sl.orisymplane = SymPlane(np.array([0, 0]),
                                              abst[0]+abst[1], abst)
                    sl.planegroup = targetsym
                else:
                    invalidDirectionMessage(rp, planegroup, targetsym)
            elif planegroup == 'p6m':   # reducing to: cm, cmm
                # exclusively do 60° rotations, set mirrors accordingly
                if (tspar[0], tspar[1]) in [
                        (1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1),
                        (1, -1), (-1, 1), (1, 2), (-1, -2), (2, 1), (-2, -1)]:
                    sl.planegroup = targetsym
                    chir = 1
                    if abs(np.dot(abst[1], np.dot(np.array(
                           [[-0.5, -np.sqrt(3)/2], [np.sqrt(3)/2, -0.5]]),
                            abst[0]))) > 0.01:
                        chir = -1   # left-handed unit cell
                    if (tspar[0], tspar[1]) in [(1, 1), (-1, -1),
                                                (1, -1), (-1, 1)]:
                        if targetsym == 'cm':
                            sl.orisymplane = SymPlane(
                                np.array([0, 0]), np.dot(tspar, abst), abst)
                    elif (tspar[0], tspar[1]) in [(1, 2), (-1, -2),
                                                  (2, 1), (-2, -1)]:
                        if (tspar[0], tspar[1]) in [(1, 2), (-1, -2)]:
                            sl.rotate_unit_cell(6*chir)  # rotate 60° clockwise
                        elif (tspar[0], tspar[1]) in [(2, 1), (-2, -1)]:
                            sl.rotate_unit_cell(-6*chir)  # rotate 60° countercl.
                        abst = sl.ab_cell.T
                        if targetsym == 'cm':
                            sl.orisymplane = SymPlane(
                                np.array([0, 0]), abst[0]-abst[1], abst)
                    elif (tspar[0], tspar[1]) in [(1, 0), (-1, 0),
                                                  (0, 1), (0, -1)]:
                        if (tspar[0], tspar[1]) in [(1, 0), (-1, 0)]:
                            sl.rotate_unit_cell(6*chir)  # rotate 60° clockwise
                        elif (tspar[0], tspar[1]) in [(0, 1), (0, -1)]:
                            sl.rotate_unit_cell(-6*chir)  # rotate 60° countercl.
                        abst = sl.ab_cell.T
                        if targetsym == 'cm':
                            sl.orisymplane = SymPlane(
                                np.array([0, 0]), abst[0]+abst[1], abst)
                else:
                    invalidDirectionMessage(rp, planegroup, targetsym)
            else:
                logger.warning("Unexpected point encountered in "
                               "findSymmetry routine: FS005")
                rp.setHaltingLevel(1)
    if targetsym == sl.planegroup:
        logger.info("The symmetry for the slab was reduced to "
                    + targetsym + ", as requested.")
    return planegroup


def enforceSymmetry(sl, rp, planegroup="fromslab",
                    movement='fromparams', rotcell=False):
    """Finds how atoms are linked to each other based on the planegroup.
    If the planegroup argument is not given, the planegroup assigned to
    the slab will be used. Otherwise, the given planegroup has to be a
    subgroup of the highest symmetry planegroup found for the slab. Sets
    movement = True or False to force or suppress recalculating atom
    positions to perfectly fit the symmetry; keep at default to follow
    SYMMETRIZE_INPUT parameter. Set rotcell = False to avoid rotating the
    unit cell such that a mirror plane is parallel to x, y, or 45°."""
    if planegroup == "fromslab":
        if sl.planegroup != "unknown":
            planegroup = sl.planegroup
        else:
            logger.warning("Call to enforceSymmetry before findSymmetry; "
                           "running findSymmetry now.")
            planegroup = sl.findSymmetry(rp)
    if movement == 'fromparams':
        nomove = not rp.SYMMETRIZE_INPUT
    else:
        if type(movement) == bool:
            nomove = not movement
        else:
            nomove = not rp.SYMMETRIZE_INPUT
            logger.warning("enforceSymmetry: Invalid 'movement' variable "
                           "passed. Using SYMMETRIZE_INPUT parameter instead.")
    eps = rp.SYMMETRY_EPS
    abst = sl.ab_cell.T  # surface unit cell, transposed

    # FIND ATOM LINKING - HERE WORK WITH sl INSTEAD OF ts, SINCE WE WANT
    #   TO ASSIGN PROPERTIES TO INDIVIDUAL ATOMS
    for at in sl:  # first put all atoms in a list of their own
        at.linklist = [at]
        at.symrefm = np.identity(2)
    if not planegroup == "p1":  # p1 has no symmetry to check for
        sl.create_sublayers(eps.z)
        sl.sort_original()
        sl.collapse_cartesian_coordinates()
        # TEST ROTATION AT ORIGIN - TESTING ONLY HIGHEST ROTATIONAL ORDER
        #    IS ENOUGH
        if planegroup not in ["p1", "pm", "pg", "cm", "rcm"]:
            if planegroup in ["p2", "pmm", "pmg", "pgg", "cmm", "rcmm"]:
                toprotsym = 2
            elif planegroup in ["p3", "p3m1", "p31m"]:
                toprotsym = 3
            elif planegroup in ["p4", "p4m", "p4g"]:
                toprotsym = 4
            elif planegroup in ["p6", "p6m"]:
                toprotsym = 6
            else:
                logger.warning("Unexpected point encountered in "
                               "enforceSymmetry routine: ES001")
                rp.setHaltingLevel(1)
            tmpslab = copy.deepcopy(sl)
            tmpslab.rotate_atoms(toprotsym)
            m = np.linalg.inv(rotation_matrix_order(toprotsym))
            for (sli, sl1) in enumerate(sl.sublayers):
                for (ati, at1) in enumerate(sl1):
                    for (atj, at2) in enumerate(tmpslab.sublayers[sli]):
                        if sl1.atlist[atj] in at1.linklist:
                            # don't check atoms that are already linked
                            continue
                        if not at1.is_same_xy(at2, eps):
                            continue
                        # combine the two linklists
                        at1.linklist.extend(sl1.atlist[atj].linklist)
                        for at3 in sl1.atlist[atj].linklist:
                            at3.symrefm = np.linalg.multi_dot([
                               np.linalg.inv(at2.symrefm),
                               at3.symrefm, m, at1.symrefm])
                            if not at3 == sl1.atlist[atj]:
                                at3.linklist = at1.linklist
                        sl1.atlist[atj].linklist = at1.linklist
                        break
        # TEST MIRROR AND GLIDE PLANES
        if planegroup not in ["p2", "p3", "p4", "p6"]:
            if planegroup in ["pm", "pg", "cm", "pmg", "rcm"]:
                # two possibilities, stored earlier
                testplane = sl.orisymplane
            elif planegroup in ["pmm", "p4m", "p6m", "rcmm"]:
                # mirror at 0
                testplane = SymPlane(np.array([0, 0]), abst[0], abst)
            elif planegroup in ["pgg", "p4g"]:
                # glide plane parallel to a at b/4
                testplane = SymPlane(abst[1]/4, abst[0], abst, ty="glide")
            elif planegroup in ["cmm", "p31m"]:
                # mirror a+b
                testplane = SymPlane(np.array([0, 0]), abst[0]+abst[1], abst)
            elif planegroup == "p3m1":
                # mirror a-b
                testplane = SymPlane(np.array([0, 0]), abst[0]-abst[1], abst)
            else:
                logger.warning("Unexpected point encountered in "
                               "enforceSymmetry routine: ES002")
                rp.setHaltingLevel(1)
            tmpslab = copy.deepcopy(sl)
            tmpslab.mirror_atoms(testplane)
            ang = angle(np.array([1, 0]), testplane.dir)
            rotm = rotation_matrix(ang)
            m = np.dot(rotm, np.dot(np.array([[1, 0], [0, -1]]),
                                    np.linalg.inv(rotm)))
            for (sli, sl1) in enumerate(sl.sublayers):
                for (ati, at1) in enumerate(sl1):
                    for (atj, at2) in enumerate(tmpslab.sublayers[sli]):
                        # don't check atoms that are already linked
                        if sl1.atlist[atj] in at1.linklist:
                            continue
                        if not at1.is_same_xy(at2, eps):
                            continue
                        # combine the two linklists
                        at1.linklist.extend(sl1.atlist[atj].linklist)
                        for at3 in sl1.atlist[atj].linklist:
                            at3.symrefm = np.linalg.multi_dot([
                                np.linalg.inv(at2.symrefm),
                                at3.symrefm, m, at1.symrefm])
                            if not at3 == sl1.atlist[atj]:
                                at3.linklist = at1.linklist
                        sl1.atlist[atj].linklist = at1.linklist
                        break
    sl.linklists = []     # re-create linklists
    for at in sl:
        if len(at.linklist) > 1 and at.linklist not in sl.linklists:
            # don't keep the linklists of length 1
            sl.linklists.append(at.linklist)
    # check that all linked atoms have the same sites
    for ll in sl.linklists:
        if len(set([at.site for at in ll])) > 1:
            logger.warning(
                "Symmetry-equivalent atoms are assigned to different sites: "
                + "; ".join(["{}: {}".format(at, at.site.label) for at in ll]))
            rp.setHaltingLevel(1)
    # FIND ALLOWED MOVE DIRECTIONS FOR ATOMS
    lockpoints = []     # full list of rotation points (no duplication)
    ori = np.array([0., 0.])
    if planegroup in ["p2", "pmm", "pmg", "pgg", "cmm", "p4", "p4m", "p4g",
                      "p6", "p6m", "rcmm"]:
        lockpoints = [ori, 0.5*abst[0], 0.5*abst[1], 0.5*(abst[0]+abst[1])]
        if planegroup in ["p6", "p6m"]:
            lockpoints.extend([(2*abst[0]+abst[1])/3,
                               (abst[0]+2*abst[1])/3])
        if planegroup == "rcmm":
            lockpoints.extend([0.25*(abst[0]+abst[1]),
                               0.75*(abst[0]+abst[1]),
                               0.25*abst[0]+0.75*abst[1],
                               0.75*abst[0]+0.25*abst[1]])
    elif planegroup in ["p3", "p3m1", "p31m"]:
        lockpoints = [ori, (2*abst[0]+abst[1])/3, (abst[0]+2*abst[1])/3]
    lockplanes = []
    # full list of mirror planes; no glide planes, since those don't restrict
    #   single atom movement
    if planegroup in ["pm", "rcm"]:
        # include duplication of planes at ori+a / ori+b to avoid having to
        #  check atom positions +- a/b
        if np.array_equal(sl.orisymplane.par, [1, 0]):
            lockplanes = [sl.orisymplane,
                          SymPlane(abst[1]/2, abst[0], abst),
                          SymPlane(abst[1], abst[0], abst, collapse=False)]
        else:
            lockplanes = [sl.orisymplane,
                          SymPlane(abst[0]/2, abst[1], abst),
                          SymPlane(abst[0], abst[1], abst, collapse=False)]
    if planegroup == "cm":
        if np.array_equal(sl.orisymplane.par, [1, 1]):
            lockplanes = [sl.orisymplane,
                          SymPlane((abst[0] - abst[1]) / 2,
                                   abst[0] + abst[1], abst, collapse=False),
                          SymPlane((abst[1] - abst[0]) / 2,
                                   abst[0] + abst[1], abst, collapse=False)]
        else:
            lockplanes = [sl.orisymplane,
                          SymPlane((abst[0] + abst[1]) / 2,
                                   abst[0] - abst[1], abst),
                          SymPlane((abst[0] + abst[1]),
                                   abst[0] - abst[1], abst, collapse=False)]
    if planegroup in ["pmm", "p4m", "rcmm"]:
        for i in range(0, 2):
            lockplanes.append(SymPlane(ori, abst[i], abst))
            lockplanes.append(SymPlane(abst[abs(i-1)]/2, abst[i], abst))
            lockplanes.append(SymPlane(abst[abs(i-1)], abst[i], abst,
                                       collapse=False))
    if planegroup in ["cmm", "p4m"]:
        lockplanes.append(SymPlane(ori, abst[0]+abst[1], abst))
        lockplanes.append(SymPlane((abst[0]+abst[1])/2, abst[0]-abst[1],
                                   abst))
        if planegroup == "cmm":
            lockplanes.extend([SymPlane(ori, abst[0]-abst[1], abst),
                               SymPlane(abst[0]+abst[1], abst[0]-abst[1],
                                        abst, collapse=False),
                               SymPlane((abst[0]-abst[1])/2, abst[0]+abst[1],
                                        abst, collapse=False),
                               SymPlane((abst[1]-abst[0])/2, abst[0]+abst[1],
                                        abst, collapse=False)])
    if planegroup == "pmg":
        if np.array_equal(sl.orisymplane.par, [1, 0]):
            lockplanes = [SymPlane(abst[0]/4, abst[1], abst),
                          SymPlane(abst[0]*3/4, abst[1], abst)]
        else:
            lockplanes = [SymPlane(abst[1]/4, abst[0], abst),
                          SymPlane(abst[1]*3/4, abst[0], abst)]
    if planegroup == "p4g":
        lockplanes = [SymPlane((abst[0]+abst[1])/4, abst[0]-abst[1], abst),
                      SymPlane((abst[0]+abst[1])*3/4, abst[0]-abst[1], abst),
                      SymPlane((abst[0]-abst[1])/4, abst[0]+abst[1], abst,
                               collapse=False),
                      SymPlane((abst[1]-abst[0])/4, abst[0]+abst[1], abst,
                               collapse=False)]
    if planegroup in ["p3m1", "p6m"]:
        lockplanes.extend([SymPlane((abst[0]+abst[1])/2, abst[0]-abst[1],
                                    abst),
                           SymPlane(ori, abst[0]+2*abst[1], abst,
                                    index2=True),
                           SymPlane(ori, 2*abst[0]+abst[1], abst,
                                    index2=True),
                           SymPlane(abst[0]+abst[1], abst[0]+2*abst[1], abst,
                                    index2=True, collapse=False),
                           SymPlane(abst[0]+abst[1], 2*abst[0]+abst[1], abst,
                                    index2=True, collapse=False)])
    if planegroup in ["p31m", "p6m"]:
        lockplanes.extend([SymPlane(ori, abst[0], abst),
                           SymPlane(ori, abst[1], abst),
                           SymPlane(abst[0], abst[1], abst, collapse=False),
                           SymPlane(abst[1], abst[0], abst, collapse=False),
                           SymPlane(ori, abst[0]+abst[1], abst)])
    ts = copy.deepcopy(sl)
    ts.project_c_to_z()
    ts.collapse_cartesian_coordinates()
    for at in ts:
        # first check points
        for p in lockpoints:
            if at.is_same_xy(p, eps):
                if not nomove:
                    at.cartpos[0:2] = p
                at.freedir = 0  # lock completely
                break
        # then if not locked yet, check planes
        if not at.freedir == 0:
            for pl in lockplanes:
                d = dist_from_line(pl.pos, pl.pos+pl.dir, at.cartpos[:2])
                if d < eps:
                    at.freedir = pl.par
                    if not nomove:  # shift atom onto plane
                        shiftv = np.array([pl.dir[1], -pl.dir[0]]) * d
                        if (dist_from_line(pl.pos, pl.pos+pl.dir,
                                           at.cartpos[:2] + shiftv)
                                > d * 1.1):
                            shiftv = -1 * shiftv
                        at.cartpos[:2] += shiftv
                    break
    for at in sl:
        at2 = [a for a in ts if a.num == at.num][0]
        at.freedir = at2.freedir
        at.cartpos = at2.cartpos
    # average positions for linked atoms
    if not nomove and not planegroup == "p1":
        sl.collapse_cartesian_coordinates()
        releps = [eps / np.linalg.norm(abst[j]) for j in range(0, 2)]
        mvslabs = []
        if planegroup not in ["pm", "pg", "cm", "rcm"]:
            for i in range(0, toprotsym-1):
                if len(mvslabs) == 0:
                    tmpslab = copy.deepcopy(sl)
                else:
                    tmpslab = copy.deepcopy(mvslabs[-1])
                tmpslab.rotate_atoms(toprotsym, ori)
                mvslabs.append(tmpslab)
        if planegroup not in ["p2", "p3", "p4", "p6"]:
            tmpslab = copy.deepcopy(sl)
            tmpslab.mirror_atoms(testplane)
            mvslabs.append(tmpslab)
            if planegroup not in ["pm", "pg", "cm", "rcm"]:
                for i in range(0, toprotsym-1):
                    tmpslab = copy.deepcopy(mvslabs[-1])
                    tmpslab.rotate_atoms(toprotsym, ori)
                    mvslabs.append(tmpslab)
        for (llind, ll) in enumerate(sl.linklists):
            for at in ll:
                psum = np.copy(at.cartpos)
                pn = 1
                for ms in mvslabs:
                    found = False
                    for atind in range(0, len(ll)):
                        at2 = ms.linklists[llind][atind]
                        complist = [at2.cartpos[:2]]
                        for j in range(0, 2):
                            if abs(at2.pos[j]) < releps[j]:
                                complist.append(at2.cartpos[:2]+abst[j])
                            if abs(at2.pos[j]-1) < releps[j]:
                                complist.append(at2.cartpos[:2]-abst[j])
                        if len(complist) == 3:
                            # corner - add the diagonally opposed one
                            complist.append(complist[1]+complist[2]
                                            - complist[0])
                        for p in complist:
                            if np.linalg.norm(p - at.cartpos[:2]) < eps:
                                psum += np.append(p, at2.cartpos[2])
                                pn += 1
                                found = True
                                break
                        if found:
                            break
                if pn != len(mvslabs) + 1:
                    logger.warning(
                        "During symmetrization of the slab, not all "
                        "symmetry-equivalent atoms were found for {}. Atom "
                        "will not be symmetrized.")
                else:
                    at.cartpos = psum / pn
    sl.collapse_cartesian_coordinates(update_origin=True)
    if not rotcell:
        # !!! THIS IS NOW THE DEFAULT, rotcell is NEVER true.
        #  if it stays like this, consider deleting the following lines..
        return
    # after everything else is done, rotate unit cell (in x,y without
    #   changing fractional coordinates) if necessary:
    # because TensErLEED is faster if mirror planes are along x, y, or x+y
    mirrordirs = []
    # if a or b are already along one of the "allowed" directions, we don't
    #   need to do anything
    if planegroup in ["cm", "cmm", "p3m1", "p31m", "p6", "p6m"]:
        # rotate to one of the diagonals
        mirrordirs.append((abst[0]-abst[1])
                          / np.linalg.norm(abst[0]-abst[1]))
        mirrordirs.append((abst[0]+abst[1])
                          / np.linalg.norm(abst[0]+abst[1]))
    elif planegroup in ["pm", "rcm", "pmm", "pmg", "rcmm", "p4m", "p4g"]:
        mirrordirs.append(abst[0]/np.linalg.norm(abst[0]))
        mirrordirs.append(abst[1]/np.linalg.norm(abst[1]))
        if planegroup in ["p4m", "p4g"]:
            mirrordirs.append((abst[0]+abst[1])
                              / np.linalg.norm(abst[0]+abst[1]))
            mirrordirs.append((abst[0]-abst[1])
                              / np.linalg.norm(abst[0]-abst[1]))
    if len(mirrordirs) > 0:
        found = False
        for d in mirrordirs:
            for i in range(0, 2):
                if abs(np.dot(np.array([1, 0]), d)-1) < 0.001:
                    found = True
        if not found:
            ang = angle(np.array([1, 0]), mirrordirs[0])
            rotm = rotation_matrix(ang, dim=3)
            sl.ucell = np.dot(rotm, sl.ucell)
            for i in range(0, 3):
                for j in range(0, 3):
                    if abs(sl.ucell[i, j]) < 1e-6:
                        sl.ucell[i, j] = 0
            sl.update_cartesian_from_fractional()
            # modify BEAM_INCIDENCE
            if rp.THETA != 0:
                logger.debug("Modifying BEAM_INCIDENCE parameter")
                rp.PHI += np.degrees(ang)
                parameters.modify(rp, "BEAM_INCIDENCE")
    return


def getSymBaseSymmetry(sl, rp):
    """Runs the symmetry search for the symbaseslab, then transfers atom
    linking to translationally equivalent atoms in the extended slab."""
    if sl.symbaseslab is None:
        logger.error("getSymBaseSymmetry: No symmetry base slab defined.")
        raise RuntimeError("getSymBaseSymmetry called without symmetry base "
                           "slab.")
    if sl.symbaseslab.planegroup == "unknown":
        findSymmetry(sl.symbaseslab, rp, forceFindOri=True)
        enforceSymmetry(sl.symbaseslab, rp, rotcell=False)
    sl.update_cartesian_from_fractional()
    for ll in sl.symbaseslab.linklists:
        newll = []
        for ssl_at in ll:
            at = [a for a in sl if a.num == ssl_at.num][0]
            newll.append(at)
            at.linklist = newll
            at.symrefm = np.copy(ssl_at.symrefm)
    for at in [at for at in sl if at.duplicate_of]:
        at.duplicate_of.linklist.append(at)
        at.linklist = at.duplicate_of.linklist
        at.symrefm = np.copy(at.duplicate_of.symrefm)
        if rp.SYMMETRIZE_INPUT:
            cv = at.cartpos[:2] - at.duplicate_of.cartpos[:2]
            v = np.append(np.round(np.dot(np.linalg.inv(
                                     sl.symbaseslab.ab_cell), cv)), 0.)
            at.cartpos = (at.duplicate_of.cartpos
                          + np.dot(sl.symbaseslab.ucell, v))
    sl.update_fractional_from_cartesian()
    sl.linklists = []
    for at in sl:
        if len(at.linklist) > 1 and at.linklist not in sl.linklists:
            sl.linklists.append(at.linklist)
