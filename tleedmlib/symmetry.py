# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 13:58:39 2020

@author: Florian Kraushofer

Functions for determining and setting slab symmetry
"""

import copy
import itertools
import logging
import re

import numpy as np

from viperleed.tleedmlib import leedbase
from viperleed.tleedmlib.base import (addUnequalPoints, angle, dist_from_line,
                                      rotation_matrix_order, rotation_matrix)
from viperleed.tleedmlib.classes.slab import SymPlane
from viperleed.tleedmlib.files import parameters

logger = logging.getLogger("tleedm.symmetry")


def findBulkSymmetry(sl, rp):
    """Checks the bulk slab for screw axes and glide planes."""
    eps = rp.SYMMETRY_EPS
    epsz = rp.SYMMETRY_EPS_Z
    uct = np.transpose(copy.copy(sl.ucell))
    abt = uct[:2, :2]
    rotsfound = []
    glidesfound = []
    ts = copy.deepcopy(sl)
    ts.sort_by_z()
    ts.collapseCartesianCoordinates()
    ts.createSublayers(epsz)
    # optimize C vector
    newC = ts.getMinC(rp)
    if newC is not None:
        logger.debug("Bulk unit cell could be reduced with repeat vector "
                     "[{:.5f} {:.5f} {:.5f}]".format(*(-newC)))
        # apply new unit cell
        ts.atlist = [at for at in ts.atlist
                     if at.cartpos[2] > ts.topat_ori_z - abs(newC[2])]
        ts.layers[0].atlist = ts.atlist
        ts.layers = [ts.layers[0]]
        ts.layers[0].isBulk = True
        rp2 = copy.deepcopy(rp)
        rp2.SUPERLATTICE = np.array([[1, 0], [0, 1]], dtype=float)
        rp2.BULK_REPEAT = -newC
        ts = ts.makeBulkSlab(rp2)
    # figure out what to check
    pcands = ts.getCandidateLayerPeriod(eps)
    if len(pcands) == 0:
        return
    nl = len(ts.sublayers)
    # check for screw axes
    checkrots = []
    if nl % 2 == 0:
        checkrots.extend([2, 4])
    if ts.celltype == "hexagonal" and (nl % 3 == 0):
        checkrots.extend([3, 6])
    for per in pcands:
        for ro in [ro for ro in checkrots if ro not in rotsfound]:
            if ts.isBulkScrewSymmetric(ro, per, eps):
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
            if ts.isBulkGlideSymmetric(gl, per, eps):
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
    epsz = rp.SYMMETRY_EPS_Z
    ori = np.array([0., 0.])
    # reduce surface unit cell
    abst = sl.ucell[:2, :2].T  # surface unit cell, transposed
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
            newsl = ("SUPERLATTICE M = {:.0f} {:.0f}, {:.0f} {:.0f}"
                     .format(*[x for y in rp.SUPERLATTICE for x in y]))
            parameters.modifyPARAMETERS(rp, "SUPERLATTICE", newsl,
                                        include_left=True)
        # MODIFY SYMMETRY_FIX PARAMETER
        if "[" in rp.SYMMETRY_FIX and not bulk:
            rgx = re.compile(r'\s*(?P<group>(pm|pg|cm|rcm|pmg))\s*\[\s*'
                             + r'(?P<i1>[-012]+)\s+(?P<i2>[-012]+)\s*\]')
            m = rgx.match(rp.SYMMETRY_FIX)
            targetsym = m.group("group")
            tspar = [int(m.group("i1")), int(m.group("i2"))]
            cartdir = np.dot(tspar, abst)
            newab = np.dot(sl.ucell[:2, :2], np.transpose(usurf))
            newdir = np.dot(np.linalg.inv(newab), cartdir)
            newdir = newdir / min(newdir)
            s = (targetsym+"[{:.0f} {:.0f}]".format(newdir[0], newdir[1]))
            parameters.modifyPARAMETERS(rp, "SYMMETRY_FIX", s)
        # MODIFY UNIT CELL
        sl.getCartesianCoordinates()
        sl.ucell_mod.append(('rmul', utr.T))
        sl.ucell = np.dot(sl.ucell, utr.T)
        # same as np.transpose(np.dot(utr,np.transpose(sl.ucell)))
        sl.collapseCartesianCoordinates(updateOrigin=True)
        # gets fractional coordinates in the new unit cell and
        #   collapses appropriately

    # check cell type again
    abst = sl.ucell[:2, :2].T
    celltype, _ = leedbase.checkLattice(abst)
    sl.celltype = celltype
    if output:
        logger.info("Found unit cell type: "+celltype)
        logger.info("Starting symmetry search...")

    # FIND HIGHEST SYMMETRY ORIGIN
    sl.collapseCartesianCoordinates()
    # create a testslab: C projected to Z
    ts = copy.deepcopy(sl)
    if bulk and any(abs(sl.ucell[:2, 2]) > eps):
        # off-normal bulk repeat vectors may break symmetry
        #  -> double bulk for symmetry testing
        ts = ts.doubleBulkSlab()
    ts.createSublayers(epsz)
    ts.projectCToZ()
    ts.sort_by_z()

    lowocc_layer = ts.getLowOccLayer()
    compare_to = ts.getCompareCoords(eps)

    rotorders = [2]
    if celltype == "square":
        rotorders.append(4)
    if celltype == "hexagonal":
        rotorders.extend([3, 6])
    # matrices needed to find the axis from two rotation-equivalent atoms:
    rot_scale_mat = {
        2: np.eye(2)*0.5,
        3: rotation_matrix(np.pi/6) / np.sqrt(3),
        4: rotation_matrix(np.pi/4) / np.sqrt(2),
        6: rotation_matrix(np.pi/3)}

    # check rotation symmetries, starting with highest
    toprotlist = []
    toprotsym = 0
    if output:
        logger.debug("Checking for rotation axes...")
    for order in rotorders[::-1]:
        if toprotsym != 0:
            break  # we only need the highest symmetry axis
        rotated_slab = copy.deepcopy(ts)
        rotated_slab.rotateAtoms(order)
        lowocc_rotated = rotated_slab.sublayers[
            ts.sublayers.index(lowocc_layer)]
        if rotated_slab.isTranslationSymmetric_2D(
                ori, eps, compare_to=compare_to):
            # origin is highest-order rotation axis, we're done here
            toprotlist = [ori]
            toprotsym = order
            break
        if not (rp.SYMMETRY_FIND_ORI or forceFindOri):
            continue
        for i in range(0, len(lowocc_layer.atlist)):
            translation = (lowocc_layer.atlist[i].cartpos[:2]
                           - lowocc_rotated.atlist[0].cartpos[:2])
            if rotated_slab.isTranslationSymmetric_2D(
                    translation, eps, compare_to=compare_to):
                v = (lowocc_layer.atlist[0].cartpos[:2]
                     - lowocc_layer.atlist[i].cartpos[:2])
                toprotsym = order
                toprotlist.append(lowocc_layer.atlist[0].cartpos[:2]
                                  + np.dot(rot_scale_mat[order], v))
    if toprotsym > 0:
        if output:
            logger.debug("Highest rotation axis has order " + str(toprotsym))
        # add implied other rotation axes
        add_axes = {
            2: ((0, 0.5), (0.5, 0), (0.5, 0.5)),
            3: ((1/3, 2/3), (2/3, 1/3)),
            4: ((0.5, 0.5),),
            6: tuple()}
        for p in toprotlist[:]:
            toprotlist.extend([p + abst[0]*i + abst[1]*j
                               for (i, j) in add_axes[toprotsym]])
    mirror = False
    glide = False
    symplanelist = []
    test_mirror_dirs = [(1, 0), (0, 1), (1, 1), (1, -1)]
    if celltype == "hexagonal":
        test_mirror_dirs.extend([(1, 2), (2, 1)])
    if not rp.SYMMETRY_FIND_ORI and not forceFindOri and not bulk:
        # check origin and selected points
        for (pa, pb) in [(0, 0), (0.25, 0.25), (0.25, -0.25)]:
            for (i, j) in test_mirror_dirs:
                symplanelist.append(SymPlane(pa*abst[0]+pb*abst[1],
                                             i*abst[0]+j*abst[1], abst,
                                             index2=True))
    else:
        # check all highest-order rotation axes
        for (k, pos) in enumerate(toprotlist):
            for (i, j) in test_mirror_dirs:
                symplanelist.append(SymPlane(pos, i*abst[0]+j*abst[1],
                                             abst, index2=True))
    if output and symplanelist:
        logger.debug(f"Testing {len(symplanelist)} candidates for "
                     "mirror/glide planes...")
    for spl in symplanelist:    # test the candidates
        if ts.isMirrorSymmetric(spl, eps):
            spl.type = "mirror"
            mirror = True
        elif ts.isMirrorSymmetric(spl, eps, glide=True):
            spl.type = "glide"
            glide = True
    if mirror:
        # drop the glides, they are implied in rotation+mirror
        symplanelist = [spl for spl in symplanelist if spl.type != "glide"]

    # no rotation axes: we had no candidates. generic search for mirrors/glides
    if toprotsym == 0 and (rp.SYMMETRY_FIND_ORI or forceFindOri or bulk):
        if output:
            logger.debug("Checking for mirror/glide planes...")
        for (i, j) in test_mirror_dirs:
            symplane = SymPlane(ori, i*abst[0]+j*abst[1], abst)
            glidevec = (symplane.par[0]*abst[0]+symplane.par[1]*abst[1])/2
            mirrored_slab = copy.deepcopy(ts)
            mirrored_slab.mirror(symplane)
            lowocc_mirrored = mirrored_slab.sublayers[
                ts.sublayers.index(lowocc_layer)]
            to_check = [(k, at.cartpos[:2]) for (k, at) in
                        enumerate(lowocc_layer.atlist)]
            # also check translations by whole unit vectors
            # this is necessary to detect the glide planes in e.g. cm
            to_check.extend([(0, lowocc_layer.atlist[0].cartpos[:2]
                             + va*abst[0] + vb*abst[1])
                             for (va, vb) in ((1, 0), (0, 1), (1, 1))])
            for (k, cartpos) in to_check:
                translation = (cartpos
                               - lowocc_mirrored.atlist[0].cartpos[:2])
                for is_glide in (True, False):
                    if (not is_glide and np.linalg.norm(translation) > eps
                            and abs(np.dot(translation, symplane.dir)) > 1e-3):
                        # we only want translations perpendicular to a mirror
                        continue
                    if is_glide:
                        for sign in (+1, -1):
                            v = translation + glidevec*sign
                            if (abs(np.dot(v, symplane.dir)) < 1e-3
                                    or np.linalg.norm(v) < eps):
                                break
                        else:
                            continue
                    if mirrored_slab.isTranslationSymmetric_2D(
                            translation, eps, compare_to=compare_to):
                        found_symplane = copy.deepcopy(symplane)
                        found_symplane.pos = (
                            lowocc_layer.atlist[0].cartpos[:2]
                            + cartpos) / 2
                        found_symplane.collapse(abst)
                        if not is_glide:
                            mirror = True
                            found_symplane.type = "mirror"
                        else:
                            glide = True
                            found_symplane.type = "glide"
                        symplanelist.append(found_symplane)
                        # also add implied 2nd plane (may be equivalent!)
                        symplanelist.append(found_symplane.get_implied(abst))

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
        for at in sl.atlist:
            at.cartpos[:2] -= topsympoint
        for at in ts.atlist:
            at.cartpos[:2] -= topsympoint
        sl.ucell_mod.append(('add', -topsympoint))
        sl.getFractionalCoordinates()
        ts.getFractionalCoordinates()
    else:
        # identify some groups with NO rotation here:
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
            for at in sl.atlist:
                at.cartpos[0:2] -= shiftv
            # ts is not used any more in this case, otherwise those atoms
            #  would have to be shifted as well.
            sl.ucell_mod.append(('add', -shiftv))
            sl.getFractionalCoordinates()
        if oriplane:
            oriplane.pos = np.array([0, 0])
            sl.orisymplane = oriplane

    if not planegroup:
        # !!! THIS SHOULD BE OBSOLETE WITH NEW ALGO - TEST!
        # start by checking special case: in cmm, there are two inequivalent
        #  2fold axes, one of which would not have been found yet -> shift
        #  there (potentially), test
        if toprotsym == 2 and celltype in ["hexagonal", "rhombic"]:
            shiftslab = copy.deepcopy(ts)
            for at in shiftslab.atlist:
                at.cartpos[:2] -= abst[0]/2
            shiftslab.getFractionalCoordinates()
            # test diagonal mirror at shifted origin
            spl = SymPlane(np.array([0, 0]), (abst[0]+abst[1]), abst)
            if shiftslab.isMirrorSymmetric(spl, eps):
                planegroup = "cmm"
                ts = shiftslab
                # correct origin
                for at in sl.atlist:
                    at.cartpos[0:2] -= abst[0]/2
                sl.ucell_mod.append(('add', -abst[0]/2))
                sl.getFractionalCoordinates()

    if not planegroup:
        efftype = ""    # effective cell type
        if celltype == "hexagonal":
            if not (toprotsym == 3 or toprotsym == 6):
                efftype = "rhombic"
            else:
                # test mirror plane along unit vector at origin
                spl = SymPlane(np.array([0, 0]), abst[0], abst)
                if ts.isMirrorSymmetric(spl, eps):
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
                            if ts.isMirrorSymmetric(spl, eps):
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
                if ts.isMirrorSymmetric(spl, eps):
                    planegroup = "p4m"
                else:
                    # test glide plane along diagonal at origin
                    spl = SymPlane(np.array([0, 0]), (abst[0]+abst[1]),
                                   abst)
                    if ts.isMirrorSymmetric(spl, eps, glide=True):
                        planegroup = "p4g"
                    else:
                        planegroup = "p4"
            else:   # p1 p2 pm pg pmm pmg pgg
                # test mirror plane along both diagonals at origin
                found = False
                for i in [+1, -1]:
                    spl = SymPlane(np.array([0, 0]), (abst[0]+i*abst[1]),
                                   abst)
                    if ts.isMirrorSymmetric(spl, eps):
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
                if ts.isMirrorSymmetric(spl, eps):
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
                if ts.isMirrorSymmetric(spl, eps):
                    mirs[i] = True
            if toprotsym == 2:
                if mirs[0] and mirs[1]:
                    planegroup = "pmm"
                else:
                    # test glide planes at origin
                    gldplane = None
                    for i in range(0, 2):
                        spl = SymPlane(np.array([0, 0]), abst[i], abst)
                        if ts.isMirrorSymmetric(spl, eps, glide=True):
                            spl.type = "glide"
                            gldplane = spl
                    if gldplane is not None:
                        planegroup = "pmg"
                        sl.orisymplane = gldplane
                    else:
                        # test glide plane at a/4:
                        spl = SymPlane(abst[0]/4, abst[1], abst)
                        if ts.isMirrorSymmetric(spl, eps, glide=True):
                            planegroup = "pgg"
                        else:
                            planegroup = "p2"
                if planegroup in ["pmm", "pmg", "pgg"]:
                    # each of these might be a mis-identified rcmm
                    if ts.isRotationSymmetric((abst[0]+abst[1])/4, 2, eps):
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
    dev_from_60 = [(a+1e-5) % np.radians(60) < 1e-3
                   for a in angles]
    if not any(dev_from_60) or all(dev_from_60):
        err = ("The POSCAR cell could not be "
               "rotated to a higher symmetry")
        logger.error(err)
        raise RuntimeError(err)
    idx = 0 if dev_from_60[0] else 1
    m = rotation_matrix(angles[idx], dim=3)
    sl.getCartesianCoordinates()
    sl.ucell_mod.append(('lmul', m))
    sl.ucell = np.dot(m, sl.ucell)
    abst = sl.ucell[:2, :2].T
    direction = (abst[0]+abst[1], abst[0]-abst[1])[idx]
    sl.collapseCartesianCoordinates(updateOrigin=True)
    oriplane = SymPlane(np.array([0, 0]), direction, abst)
    logger.warning("The POSCAR unit cell was changed to a higher "
                   "symmetry form. Make sure to check beam labels for "
                   "compatibility with new unit cell.")
    rp.checklist.append("Check modified unit cell")
    rp.setHaltingLevel(2)
    return abst, oriplane


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

    abst = np.transpose(sl.ucell[:2, :2])  # surface unit cell, transposed
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
            tspar = np.dot(np.linalg.inv(sl.ucell[:2, :2]),
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
                    for at in sl.atlist:
                        at.cartpos[:2] -= shiftv
                    sl.ucell_mod.append(('add', -shiftv))
                    sl.getFractionalCoordinates()
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
                    for at in sl.atlist:
                        at.cartpos[:2] -= shiftv
                    sl.ucell_mod.append(('add', -shiftv))
                    sl.getFractionalCoordinates()
                    sl.orisymplane = SymPlane(
                        np.array([0, 0]), np.array(sl.orisymplane.dir[1],
                                                   -sl.orisymplane.dir[0]),
                        abst)
                sl.planegroup = targetsym
            elif planegroup == 'pgg':   # reducing to: pg
                if (tspar[0], tspar[1]) in [(1, 0), (0, 1), (-1, 0), (0, -1)]:
                    shiftv = 0.25*np.dot(np.array(tspar[1], -tspar[0]), abst)
                    for at in sl.atlist:
                        at.cartpos[:2] -= shiftv
                    sl.ucell_mod.append(('add', -shiftv))
                    sl.getFractionalCoordinates()
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
                            for at in sl.atlist:
                                at.cartpos[0:2] -= shiftv
                            sl.ucell_mod.append(('add', -shiftv))
                            sl.getFractionalCoordinates()
                            sl.orisymplane.type = "glide"
                        sl.planegroup = targetsym
                    elif targetsym == "pmg":
                        shiftv = 0.25*(abst[0]+abst[1])
                        for at in sl.atlist:
                            at.cartpos[0:2] -= shiftv
                        sl.ucell_mod.append(('add', -shiftv))
                        sl.getFractionalCoordinates()
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
                    for at in sl.atlist:
                        at.cartpos[:2] -= shiftv
                    sl.ucell_mod.append(('add', -shiftv))
                    sl.getFractionalCoordinates()
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
                        sl.rotateUnitCell(6*chir)  # rotate 60° clockwise
                    elif (tspar[0], tspar[1]) in [(2, 1), (-2, -1)]:
                        sl.rotateUnitCell(-6*chir)  # rotate 60° countercl.
                    abst = sl.ucell[:2, :2].T
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
                        sl.rotateUnitCell(6*chir)  # rotate 60° clockwise
                    elif (tspar[0], tspar[1]) in [(0, 1), (0, -1)]:
                        sl.rotateUnitCell(-6*chir)  # rotate 60° countercl.
                    abst = sl.ucell[:2, :2].T
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
                            sl.rotateUnitCell(6*chir)  # rotate 60° clockwise
                        elif (tspar[0], tspar[1]) in [(2, 1), (-2, -1)]:
                            sl.rotateUnitCell(-6*chir)  # rotate 60° countercl.
                        abst = sl.ucell[:2, :2].T
                        if targetsym == 'cm':
                            sl.orisymplane = SymPlane(
                                np.array([0, 0]), abst[0]-abst[1], abst)
                    elif (tspar[0], tspar[1]) in [(1, 0), (-1, 0),
                                                  (0, 1), (0, -1)]:
                        if (tspar[0], tspar[1]) in [(1, 0), (-1, 0)]:
                            sl.rotateUnitCell(6*chir)  # rotate 60° clockwise
                        elif (tspar[0], tspar[1]) in [(0, 1), (0, -1)]:
                            sl.rotateUnitCell(-6*chir)  # rotate 60° countercl.
                        abst = sl.ucell[:2, :2].T
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
    epsz = rp.SYMMETRY_EPS_Z
    abst = sl.ucell[:2, :2].T  # surface unit cell, transposed

    # FIND ATOM LINKING - HERE WORK WITH sl INSTEAD OF ts, SINCE WE WANT
    #   TO ASSIGN PROPERTIES TO INDIVIDUAL ATOMS
    for at in sl.atlist:  # first put all atoms in a list of their own
        at.linklist = [at]
        at.symrefm = np.identity(2)
    if not planegroup == "p1":  # p1 has no symmetry to check for
        sl.createSublayers(epsz)
        sl.sortOriginal()
        sl.collapseCartesianCoordinates()
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
            tmpslab.rotateAtoms(toprotsym)
            tmpslab.collapseCartesianCoordinates()
            m = np.linalg.inv(rotation_matrix_order(toprotsym))
            for (sli, sl1) in enumerate(sl.sublayers):
                for (ati, at1) in enumerate(sl1.atlist):
                    for (atj, at2) in enumerate(tmpslab.sublayers[sli]
                                                .atlist):
                        if sl1.atlist[atj] in at1.linklist:
                            # don't check atoms that are already linked
                            continue
                        if not at1.isSameXY(at2.cartpos[:2], eps):
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
            g = True if testplane.type == "glide" else False
            tmpslab = copy.deepcopy(sl)
            tmpslab.mirror(testplane, glide=g)
            tmpslab.collapseCartesianCoordinates()
            ang = angle(np.array([1, 0]), testplane.dir)
            rotm = rotation_matrix(ang)
            m = np.dot(rotm, np.dot(np.array([[1, 0], [0, -1]]),
                                    np.linalg.inv(rotm)))
            for (sli, sl1) in enumerate(sl.sublayers):
                for (ati, at1) in enumerate(sl1.atlist):
                    for (atj, at2) in enumerate(tmpslab.sublayers[sli]
                                                .atlist):
                        # don't check atoms that are already linked
                        if sl1.atlist[atj] in at1.linklist:
                            continue
                        if not at1.isSameXY(at2.cartpos[:2], eps):
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
    for at in sl.atlist:
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
    ts.projectCToZ()
    ts.collapseCartesianCoordinates()
    for at in ts.atlist:
        # first check points
        for p in lockpoints:
            if at.isSameXY(p, eps):
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
    for at in sl.atlist:
        at2 = [a for a in ts.atlist if a.oriN == at.oriN][0]
        at.freedir = at2.freedir
        at.cartpos = at2.cartpos
    # average positions for linked atoms
    if not nomove and not planegroup == "p1":
        sl.collapseCartesianCoordinates()
        releps = [eps / np.linalg.norm(abst[j]) for j in range(0, 2)]
        mvslabs = []
        if planegroup not in ["pm", "pg", "cm", "rcm"]:
            for i in range(0, toprotsym-1):
                if len(mvslabs) == 0:
                    tmpslab = copy.deepcopy(sl)
                else:
                    tmpslab = copy.deepcopy(mvslabs[-1])
                tmpslab.rotateAtoms(toprotsym)
                tmpslab.collapseCartesianCoordinates()
                mvslabs.append(tmpslab)
        if planegroup not in ["p2", "p3", "p4", "p6"]:
            tmpslab = copy.deepcopy(sl)
            tmpslab.mirror(testplane, glide=g)
            tmpslab.collapseCartesianCoordinates()
            mvslabs.append(tmpslab)
            if planegroup not in ["pm", "pg", "cm", "rcm"]:
                for i in range(0, toprotsym-1):
                    tmpslab = copy.deepcopy(mvslabs[-1])
                    tmpslab.rotateAtoms(toprotsym)
                    tmpslab.collapseCartesianCoordinates()
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
    sl.collapseCartesianCoordinates(updateOrigin=True)
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
            sl.getCartesianCoordinates()
            # modify BEAM_INCIDENCE
            if rp.THETA != 0:
                logger.debug("Modifying BEAM_INCIDENCE parameter")
                rp.PHI += np.degrees(ang)
                parameters.modifyPARAMETERS(
                    rp, "BEAM_INCIDENCE",
                    "{:.3f} {:.3f}".format(rp.THETA, rp.PHI)
                    )
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
    sl.getCartesianCoordinates()
    for ll in sl.symbaseslab.linklists:
        newll = []
        for ssl_at in ll:
            at = [a for a in sl.atlist if a.oriN == ssl_at.oriN][0]
            newll.append(at)
            at.linklist = newll
            at.symrefm = np.copy(ssl_at.symrefm)
    for at in [at for at in sl.atlist if at.duplicateOf is not None]:
        at.duplicateOf.linklist.append(at)
        at.linklist = at.duplicateOf.linklist
        at.symrefm = np.copy(at.duplicateOf.symrefm)
        if rp.SYMMETRIZE_INPUT:
            cv = at.cartpos[:2] - at.duplicateOf.cartpos[:2]
            v = np.append(np.round(np.dot(np.linalg.inv(
                                     sl.symbaseslab.ucell[:2, :2]), cv)), 0.)
            at.cartpos = (at.duplicateOf.cartpos
                          + np.dot(sl.symbaseslab.ucell, v))
    sl.getFractionalCoordinates()
    sl.linklists = []
    for at in sl.atlist:
        if len(at.linklist) > 1 and at.linklist not in sl.linklists:
            sl.linklists.append(at.linklist)
