"""
=======================================
   ViPErLEED Graphical User Interface
=======================================
 *** module guilib.leedsim.classes ***

Created: 2010-01-12
Author: Michele Riva

Contains Qt-independent classes used by widgets to display the real space
lattices and the LEED pattern
"""

from fractions import Fraction
import re
import itertools

import numpy as np
from scipy import spatial as sp
from matplotlib import cm  # color maps
from matplotlib import colors as mpl_colors

import guilib as gl


degrees = '\u00b0'


class RealSpace():
    def __init__(self, params):
        gl.check_leed_params(params)  # this raises exceptions

        self.superlatticeM = params['SUPERLATTICE']
        # self.surfGroup = gl.PlaneGroup(params['surfGroup'])

        # Set up the lattices
        surfB = params['surfBasis']

        # self.fov is the real space field of view.
        # contains at least 4 unit cells
        self.fov = 4.2*max(surfB.ravel())
        self.surf = gl.Lattice(surfB,
                               space='real',
                               group=params['surfGroup'],
                               limit=self.fov)
        bulkB = self.get_BulkBasis()
        self.bulk = gl.Lattice(bulkB,
                               space='real',
                               group=params['bulkGroup'],
                               limit=self.fov)
        self.bulk.group.screws_glides = (params.get('bulk3Dsym', 'None'),
                                         self.bulk.type)

    def get_BulkBasis(self):
        return np.dot(np.linalg.inv(self.superlatticeM), self.surf.basis)

    def get_bulkToX(self, direction):
        if direction not in (0, 1):
            raise ValueError("First positional argument of get_bulkToX() "
                             f"should be 0 or 1. Got {direction!r} instead.")

        basis = self.bulk.basis[direction]
        theta = np.arctan2(basis[1]/np.linalg.norm(basis),
                           basis[0]/np.linalg.norm(basis))

        return np.degrees(theta)


class LEEDPattern():
    def __init__(self, params):
        self.maxEnergy = params['eMax']
        aperture = params.get('screenAperture', 110)
        self.maxScreenRadius = self.get_ScreenRadius(self.maxEnergy, aperture)
        self.superlatticeM = params['SUPERLATTICE']
        self.surfGroup = gl.PlaneGroup(params['surfGroup'])

        self.build_Lattices(params)
        self.build_LEEDPattern()

    def build_Lattices(self, params):
        # Create dummy surface lattice just to get the reciprocal basis
        surf = gl.Lattice(params['surfBasis'])
        surfB = surf.reciprocal_basis()

        # Now make the reciprocal lattices:
        self.surfR = gl.Lattice(surfB,
                                space='reciprocal',
                                group=params['surfGroup'],
                                limit=self.maxScreenRadius)
        bulkB = self.get_BulkBasis()
        self.bulkR = gl.Lattice(bulkB,
                                space='reciprocal',
                                group=params['bulkGroup'],
                                limit=self.maxScreenRadius)
        self.bulkR.group.screws_glides = (params.get('bulk3Dsym', 'None'),
                                          self.bulkR.type)
        self.nBeams = self.surfR.nbeams()

    def build_LEEDPattern(self):

        # 1) get the operations that generate domains
        self.domOps = self.get_DomainOperations()
        self.nDoms = len(self.domOps)
        self.domSuperlattices = np.array([np.dot(self.superlatticeM, op)
                                          for op in self.domOps])

        # 2) generate the kx ky coordinates for all the domains
        self.doms = [np.dot(np.dot(self.surfR.hk,
                                   np.linalg.inv(superlattice).transpose()),
                            self.bulkR.basis)
                     for superlattice in self.domSuperlattices]

        # 3) prepare the full LEED pattern of the surface, including all
        # domains
        surfLEEDx = np.array([self.doms[i][:, 0]
                              for i in range(self.nDoms)]).flatten()
        surfLEEDy = np.array([self.doms[i][:, 1]
                              for i in range(self.nDoms)]).flatten()
        self.surfLEED = np.array(list(zip(surfLEEDx, surfLEEDy)))

        # 4) find superposed beams
        self.superposedBeams = self.get_supSpots()

        # 5) prepare the formatted names of the fractional indices
        self.names = np.array([
            self.formatFractionalIndices(self.surfR.hk, m)
            for m in self.domSuperlattices])

        # 6) and build the subpatterns
        self.build_subpatterns()

    def get_supSpots(self):
        # Group all the beams into separate groups, such that each group
        # contains the indices of overlapping beams (can be one index if non
        # overlapping)
        tree = sp.cKDTree(self.surfLEED)
        return np.array(tree.query_ball_point(x=self.surfLEED, r=1e-8))

    def get_beamGrouping(self, domains=None):
        self.firstDomRefs = self.get_FirstDomainSymmetry()

        if domains is None:
            domains = range(self.nDoms)
            # THE LINES BELOW ARE PROBABLY WRONG FOR A LIST OF DOMAINS
            # nDoms = self.nDoms
        # else:
            # nDoms = len(domains)

        # beamsRefs contains at index [i, j] a list of the indices of beams of
        # domain domains[i] equivalent to beam j (including j)
        beamsRefs = []
        # for dom in range(nDoms):  # THIS IS PROBABLY WRONG FOR A LIST OF DOMAINS
        for dom in domains: 
            beamsRefs.append([crossRef + dom*self.nBeams
                              for crossRef in self.firstDomRefs])
        return np.array(beamsRefs, dtype=object)

    def get_equivalentSpots(self, domains=None):
        if not hasattr(domains, '__len__') and domains is not None:
            raise TypeError("Invalid type for domains. "
                            "Should be None or an array-like.")
        if domains is None:
            # output spots from all domains
            domains = range(self.nDoms)
        elif not all(isinstance(dom, int) for dom in domains):
            raise ValueError("Invalid domain index. "
                             "All indices should be integers")
        elif any(dom < 0 or dom >= self.nDoms for dom in domains):
            raise ValueError("Domain index out of range. "
                             f"Indices should be between 0 and {self.nDoms}.")

        names = self.names[domains]
        beamRefs = self.get_beamGrouping(domains)

        # get ready to prepare one dictionary per each domain:
        # key: beam name
        # value: names of beams equivalent to key
        #        Notice that for each key (=beam) in one domain, one should
        #        also include as equivalent beams those from other domains
        #        at the same positional index.
        #        !!!THIS IS A PROBLEM FOR NON-NORMAL INCIDENCE!!!
        #        --> get this list (the same for all domains) ready first
        allRefs = [np.concatenate(beamRefs[:, beam]).astype(int)
                   for beam in range(self.nBeams)]
        refNames = [set(self.names.ravel()[beamRef]) for beamRef in allRefs]

        # now pack the dictionaries
        # Notice that it is necessary to have one dictionary per domain since
        #   some domains will give rise to a subset of all the LEED spots:
        #   e.g., in a (2x1) on square bulk, one domain will give spot (1/2, 0)
        #   the other domain spot (0, 1/2)
        # Having all the dictionaries allows to determine which of the
        #   superposed spots are indeed equivalent (see next comment) as well
        #   as which domains overlap
        domsDicts = [dict(zip(domNames, refNames)) for domNames in names]

        # now figure out which of the beams are symmetry equivalent
        # taking superposition of spots into account:
        # spots originating from beams that are symmetry equivalent within
        # each domain, are symmetry equivalent only if all domains contribute
        # with symmetry-equivalent beams in all spots.
        #
        # For example (hexagonal bulk):
        # dom1: (1 0), (-1 0), (1 -1), (-1 1) equivalent
        # dom2: (1 0), (-1 0), (0 -1), (0 1)  equivalent
        # dom3: (1 0), (-1 0)                 equivalent
        # --> (1 0) and (-1 0) are equivalent
        #
        # dom1: (1 0),  (-1 0),  (1 -1), (-1 1)    equivalent
        # dom2:                  (1 -1), (-1 1)    equivalent
        # dom3: (0, 1), (0, -1), (1 -1), (-1 1)    equivalent
        # --> (1 -1) and (-1 1) are equivalent, but are
        #     not equivalent to [(1 0), (-1 0)]
        eqBeams = []
        for name in names.ravel():
            setLst = [dDict[name] for dDict in domsDicts if name in dDict]
            if setLst:
                newBeams = set.intersection(*setLst)
                # remove the keys of the beams that have
                # already been processed from each domain
                for domDict in domsDicts:
                    for beam in newBeams:
                        if beam in domDict:
                            del domDict[beam]
                eqBeams.append(newBeams)

        # sort within each equivalence group
        lstBeams = [sorted(list(beams), key=self.beamsSortCriterion)
                    for beams in eqBeams]
        # and by energy
        sortedBeams = sorted(lstBeams, key=self.sortEnergy)

        beamsWithIndices = []
        overlapDoms = []
        for (i, beams) in enumerate(sortedBeams):
            for beam in beams:
                overlapDoms = [dom + 1 for dom in domains
                               if beam in self.names[dom].ravel()]
                beamsWithIndices.append([beam, i, overlapDoms])

        # Now find which of the beams are extinct due to glide, place them
        # in a list of beam indices in the same order as in self.names[i]
        extinct = self.get_ExtinctFirstDomain()

        extinctNames = self.names[:, extinct]
        # notice that extinctNames is an empty list if there are no glide
        # extinct beams.
        for b, (beam, group, doms) in enumerate(beamsWithIndices):
            extDoms = [dom for dom in doms if beam in extinctNames[dom-1]]
            # Notice the dom-1 because of the way the doms list is
            # built in the previous for loop

            beamsWithIndices[b] = [beam, group, doms, extDoms]
            # now, only in case extDoms is of length 1, i.e., one domain
            # only gives rise to the extinct beam, the group index goes
            # negative
            # !!! Perhaps will be changed with the better version in base.py !!!
            if len(extDoms) == 1:
                beamsWithIndices[b][1] *= -1

        return beamsWithIndices

    def beamsSortCriterion(self, beamstr):
        beam = eval(beamstr)  # returns a tuple

        sortH = beam[0]
        sortK = beam[1]

        return (-sortH-sortK, -sortH)
    
    def sortEnergy(self, beamlist):
        beam = eval(beamlist[0])  # returns a tuple
        # sort by Energy, i.e., by g^2:
        return np.linalg.norm(np.dot(beam, self.bulkR.basis))

    def build_subpatterns(self):
        # Now sort the sub-patterns:
        # - There will always be a firstLEED and a domsLEED
        #   This is to handle the plotting of the first domain only.
        #   All entries in these lists are of type LEEDsubpattern
        #
        # - firstLEED, contains at most two entries: the first one
        #   includes all non-extinct beams, the second one
        #   the extinct ones, if any are present
        # - domsLEED is empty if only one domain is present.
        #   Otherwise, it contains:
        #   * one entry with superimposed spots
        #   * at most two entries for each domain.
        #     first entry: non-extinct; second entry: extinct (if any)

        doms = self.doms
        symEq = self.get_beamGrouping()
        names = self.names

        maskFirst = np.full(self.nBeams, True)
        maskFirst[self.get_ExtinctFirstDomain()] = False

        self.firstLEED = [LEEDsubpattern(doms[0][maskFirst],
                                         symEq[0][maskFirst],
                                         names[0][maskFirst])]
        if all(maskFirst):  # no glide-extinct beams
            pass
        else:
            glideAlpha = 0.2
            glideScale = 2
            maskExtinct = np.invert(maskFirst)
            self.firstLEED.append(LEEDsubpattern(doms[0][maskExtinct],
                                                 symEq[0][maskExtinct],
                                                 names[0][maskExtinct],
                                                 color='k',
                                                 alpha=glideAlpha,
                                                 marker='x',
                                                 sizeScale=glideScale
                                                 ))

        self.domsLEED = []
        self.domColors = None
        if self.nDoms > 1:
            # If there are domains, there always will be overlapping spots,
            # at least those at integer orders

            # superposedMask is True at all the beam indices of beams that
            # are superposed to others
            superposedMask = np.array([len(sup) > 1
                                       for sup in self.superposedBeams])

            # figure out which domains superpose:
            overlapDoms = [[x//self.nBeams + 1 for x in sup]
                           for sup in self.superposedBeams[superposedMask]]
            self.domsLEED = [LEEDsubpattern(self.surfLEED[superposedMask],
                                            symEq.ravel()[superposedMask],
                                            names.ravel()[superposedMask],
                                            domain=-1,
                                            color='gray',
                                            overlapDoms=overlapDoms
                                            )]

            # define colors
            self.domColors = cm.gnuplot(np.linspace(0.1, 0.9, self.nDoms))

            for dom in range(self.nDoms):
                mask = np.invert(
                    superposedMask[dom*self.nBeams:(dom+1)*self.nBeams])
                maskDom = np.logical_and(maskFirst, mask)
                self.domsLEED.append(LEEDsubpattern(doms[dom][maskDom],
                                                    symEq[dom][maskDom],
                                                    names[dom][maskDom],
                                                    domain=dom+1,
                                                    color=self.domColors[dom]
                                                    )
                                     )
                if all(maskFirst):  # no glide extinct
                    pass
                else:
                    maskDomExtinct = np.logical_and(maskExtinct, mask)
                    self.domsLEED.append(
                        LEEDsubpattern(doms[dom][maskDomExtinct],
                                       symEq[dom][maskDomExtinct],
                                       names[dom][maskDomExtinct],
                                       domain=dom+1,
                                       color=self.domColors[dom],
                                       alpha=glideAlpha,
                                       marker='x',
                                       sizeScale=glideScale))

    def get_BulkBasis(self):
        return np.dot(self.superlatticeM.transpose(), self.surfR.basis)

    def get_ScreenRadius(self, en, aperture=110.0):
        """
        aperture: float, default=110.0
                  degrees of aperture of the solid angle captured by the
                  LEED screen. The current (2020-08-17) default value is taken
                  from the dimensions of the screen of the ErLEED optics. The
                  MCP SpectaLEED from Omicron appears to have an equivalent
                  aperture of ~80°
        """
        elMass = 9.109e-31  # kg
        elQ = 1.60218e-19   # C
        hbar = 1.05457e-34  # J*s
        # The next one is the prefactor in AA^-1 eV^(-1/2)
        rt2me_hbar = np.sqrt(2*elMass*elQ) / hbar*1e-10

        return rt2me_hbar*np.sqrt(en)*np.sin(np.radians(aperture)/2)

    def get_FirstDomainSymmetry(self):
        """
        Returns
        -------
        crossRefs: list of np.arrays, length == len(self.surfR.hk)
                   Each np.array contains the positional indices of the
                   beams equivalent to the current one (including
                   self-reference)
        """

        hk = self.surfR.hk
        N = len(hk)
        symOps = self.surfR.group.group_ops()

        # at first, use sets to avoid duplicates
        crossRefs = [{i} for i in range(N)]

        for op in symOps:
            hkT = np.dot(op, hk.transpose()).transpose().tolist()
            # transform hk according to the operation of the group
            for i in range(N):
                idx = hkT.index(hk[i].tolist())
                # find which transformed hk is equal to the current beam
                # NB: there's always only one for each (beam, operation) pair
                #
                # and add it to the set of cross references for the current
                # beam
                crossRefs[i] |= {idx}

        # then convert each set to a np.array of int, so that it's iterable
        # and it's easy to process later
        return np.array([np.array(list(x), dtype=int) for x in crossRefs], 
                        dtype=object)   # TODO: check dtype=object

    def get_ExtinctFirstDomain(self):
        """
        Returns
        -------
        list of indices of beams of the first domain that are extinct
        due to glide
        """
        hk = self.surfR.hk
        group = self.surfR.group.group

        if 'g' in group:  # group has glide
            if '[1 0]' in group:  # pg or pmg
                idx = np.where([x[1] == 0 and abs(x[0] % 2) == 1
                                for x in hk])[0].tolist()
            elif '[0 1]' in group:  # pg or pmg
                idx = np.where([x[0] == 0 and abs(x[1] % 2) == 1
                                for x in hk])[0].tolist()
            else:  # pgg or p4g
                idx = np.where([(x[1] == 0 and abs(x[0] % 2) == 1)
                                or (x[0] == 0 and abs(x[1] % 2) == 1)
                                for x in hk])[0].tolist()
        else:
            idx = []

        return idx

    def get_DomainOperations(self):
        """
        Finds the symmetry operations of the bulk lattice that give distinct
        domains in LEED.
        NB: this function gives different results than what
        one would get from LEEDPat, as we're interested also in the symmetry
        relations between the intensities of the spots, while LEEDPat cares
        only about the presence or not of any spot.

        Returns
        -------
        list of operations

        -------
        """
        # New version written on 2020-06-22. The old version is kept at the end
        # in commented form.
        #
        # The current version is based on the concept of co-sets of a group.
        # Given a group G and a subgroup H, the left co-set of H with respect
        # to the group operation g of G is
        #           gH = {g*h : h in H}.
        # Additionally, we take into account that each element of G is found in
        # exactly only one co-set (e.g., H is the identity co-set).
        #
        # The operations g_i that generate distinct co-sets are those that will
        # generate distinct lattices

        def __array_2_tuple(arr):
            """
            Convenience function that is used to convert a 2x2 array to a 2x2
            tuple
            """
            return tuple(map(tuple, arr))

        # 1) project the operations from the surface group to the bulk; also
        #    rounding to integers
        transform = np.linalg.inv(self.superlatticeM)
        surf_ops = tuple(op.round().astype(int)
                         for op
                         in self.surfR.group.transform(transform))

        # 2) keep track of the operations of the bulk group that are already in
        #    a co-set, and of the group operations giving distinct co-sets
        coset_ops = set()
        _ops = []

        # 3) run through the bulk operations, adding the operations of new
        #    co-sets to the set above
        bulk_ops = self.bulkR.group.group_ops(include_3d=True)

        for b_op in bulk_ops:
            if b_op in coset_ops:
                continue
            _ops.append(b_op)
            coset = set(__array_2_tuple(np.dot(b_op, surf_op))
                        for surf_op in surf_ops)
            coset_ops.update(coset)

        return _ops

        # ********************** OLD VERSION STARTS HERE **********************
        # group = self.bulkR.group.group_ops()
        # M = self.superlatticeM

        # _ops = []
        # allOps = self.bulkR.group.allOps

        # for op in group:
            # if any(np.array_equal(np.dot(M, op), np.dot(Ei, np.dot(M, x)))
                   # for Ei in allOps for x in _ops):
                # pass
            # else:
                # _ops.append(op)

        # return _ops

    def formatFractionalIndices(self, hk, m):
        """
        Method that generates fractional indices of a superstructure m whose
        integer indices are hk.

        Parameters
        ----------
        hk, 2D-array
            Integer indices along the two directions to be formatted

        m, 2x2 np.array
            superlattice matrix

        Returns
        -------
        names, 1D array
            each element is a string with formatted names for the indices
        """

        # The beam indices [hb,kb] for the surface indices [h,k] are
        #    [hb, kb] = np.dot([h, k], G),
        # with G = np.linalg.inv(M).transpose()
        hkbulk = np.dot(np.linalg.inv(m), hk.transpose()).transpose()

        # now get the formatting done:
        # all fractional indices are of type nn/mu with nn integer, and
        mu = abs(np.round(np.linalg.det(m)).astype(int)).item()

        # notice that, to prevent Fraction().limit_denominator() from
        # overflowing I'm casting the numpy-type numbers to their native
        # python type via .item()
        names = [', '.join(str(Fraction(hh.item()).limit_denominator(mu))
                           for hh in hkb) for hkb in hkbulk]

        return np.array(names)

    def rotatePattern(self, angle, which='surf'):
        if which not in ['surf', 'surface', 'bulk']:
            raise ValueError("Only 'surf', 'surface', or 'bulk' are "
                             "acceptable values for rotatePattern")
        if hasattr(angle, '__len__'):
            raise ValueError('Angle input is not a scalar')

        if which in ['surf', 'surface']:
            angle = np.radians(angle)
            rot = np.array([[np.cos(angle), np.sin(angle)],
                            [-np.sin(angle), np.cos(angle)]])
            for pat in [*self.firstLEED, *self.domsLEED]:
                pat.transformCoords(rot)
            return None
        return self.bulkR.get_rotated_lattice(angle)


class LEEDsubpattern():
    """
    --> I originally intended to use this guy for the
         determination of beam equivalence and
         handling of hovering annotations
    """
    def __init__(self, beams, symBeams, names, domain=1, color='k',
                 alpha=None, marker='o', sizeScale=1, overlapDoms=None):
        if not all([len(beams) == len(x) for x in [symBeams, names]]):
            raise ValueError('Input shape mismatch. The length of symBeams'
                             'and/or _names differs from the one of beams.')
        if not isinstance(domain, int):
            raise ValueError("domain must be an integer.")

        self.beams = beams
        self.rotBeams = self.beams
        self.symmetricBeams = symBeams
        self.names = names
        self.domain = domain
        if alpha is None:
            self.color = np.array([mpl_colors.to_rgba(color, alpha=1)])
        else:
            self.color = np.array([mpl_colors.to_rgba(color, alpha=alpha)])
        self.marker = marker
        self.sizeScale = sizeScale
        self.overlappingDomains = overlapDoms

    def transformCoords(self, transf):
        if np.equal(np.shape(transf), (2, 2)).all():
            self.rotBeams = np.dot(self.beams, transf)
        else:
            raise


class Woods():
    common = {'p(1\u00d71)', 'p(1\u00d72)', 'p(2\u00d71)', 'p(2\u00d72)',
              'p(3\u00d71)', 'p(3\u00d72)', 'p(3\u00d73)', 'c(2\u00d72)',
              'c(4\u00d74)', 'c(6\u00d72)', 'c(8\u00d72)'}
    examples = {
        'Oblique': common,
        'Square': common | {'c(4\u00d72)',
                            'p(\u221a2\u00d7\u221a2)R45' + degrees,
                            'p(\u221a5\u00d7\u221a5)R26.6' + degrees,
                            'p(2\u221a2\u00d7\u221a2)R45' + degrees,
                            'c(3\u221a2\u00d7\u221a2)R45' + degrees,
                            'c(5\u221a2\u00d7\u221a2)R45' + degrees},
        'Rectangular': common,
        'Hexagonal': common | {'c(4\u00d72)',
                               'p(\u221a3\u00d7\u221a3)R30' + degrees,
                               'p(\u221a7\u00d7\u221a7)R19.1' + degrees,
                               'p(2\u221a3\u00d72\u221a3)R30' + degrees},
        'Rhombic': common,
        }
    
    def woodsToMatrix(self, woods, bulk):
        parsed = self.parseWoods(woods)
        if parsed is None:
            return None
        g1 = parsed[1]
        g2 = parsed[2]
        alpha = np.radians(parsed[3])
        R1 = np.linalg.norm(bulk[0])
        R2 = np.linalg.norm(bulk[1])
        q = R2/R1
        omega = np.arccos(np.dot(bulk[0], bulk[1])/(R1*R2))
        
        m = np.array([[g1*np.sin(omega-alpha), g1*np.sin(alpha)/q],
                      [-g2*q*np.sin(alpha), g2*np.sin(omega+alpha)]])
        m = m/np.sin(omega)
        
        if parsed[0] == 'c':
            m = np.dot([[1, 1], [-1, 1]], m)/2
        
        if self.isCommensurate(m):
            m = np.array([int(np.round(mij)) 
                          for mij in m.ravel()]).reshape(m.shape)
        else:
            m = None  # incommensurate
        return m
    
    def parseWoods(self, woods):
        if not isinstance(woods, str):
            return None
        
        reWoods = re.compile(
            r'''^(?P<prefix>[pc])                # * primitive or centered
            \(                                   # * open parenthesis
            (?P<g1>\d+)?                         # * direction1, integer part
            (\u221a(?P<g1rt>\d+))?               # * direction1, radical part
            \u00d7                               # * times
            (?P<g2>\d+)?                         # * direction2, integer part
            (\u221a(?P<g2rt>\d+))?               # * direction2, radical part
            \)                                   # * close parenthesis
            ((R(?P<alpha>\d+(\.\d+)?))\u00b0)?$  # * rotation angle
            ''', re.VERBOSE)
        # notice that the ^ and $ anchors make sure that the string is 
        # matched as a whole
        
        # check if it matches the full woods
        m = reWoods.match(woods)
        if m is not None:
            w = m.groupdict()
            g1 = 1
            g2 = 1
            if w['g1'] is not None:
                g1 *= float(w['g1'])
            if w['g1rt'] is not None:
                g1 *= np.sqrt(float(w['g1rt']))
            if w['g2'] is not None:
                g2 *= float(w['g2'])
            if w['g2rt'] is not None:
                g2 *= np.sqrt(float(w['g2rt']))
            if w['alpha'] is not None:
                alpha = float(w['alpha'])
            else:
                alpha = 0
            return [w['prefix'], g1, g2, alpha]
        return None
    
    def matrixToWoods(self, m, bulk):
        m = np.array(m)
        
        pc = self.primitiveOrCentered(m, bulk)
        if pc:  # matrix is woods representable as primitive or centered
            toFormat = []
            for gSq in [pc[1]**2, pc[2]**2]:
                gSq = np.round(gSq)
                (g2, grt) = self.squareToProdOfSquares(gSq)
                integer = str(int(round(np.sqrt(g2))))
                
                dirStr = '' #format the direction in here
                if grt > 1: # insert root part
                    dirStr = '\u221a' + str(int(grt))
                if dirStr == '': # if there is no root part, always insert 
                                 # the integer part
                    dirStr = integer
                else:
                    if integer != '1': # otherwise place it in only if it's 
                                       # not 1
                        dirStr = integer + dirStr
                
                toFormat.append(dirStr)
            woods = '{}({})'.format(pc[0], '\u00d7'.join(toFormat))
            cosAlpha = pc[3]
            if np.abs(cosAlpha) > 1e-3 and 1-np.abs(cosAlpha) > 1e-3:
                #angle not 0, 90 nor 180
                alpha = np.round(np.degrees(np.arccos(cosAlpha)),
                                 decimals = 1)
                woods += 'R' + str(alpha) + degrees
            return woods
        return None
    
    def isCommensurate(self, m, eps=1e-3):
        if m is None:
            return False
        
        mu = np.linalg.det(m)
        if np.round(mu) == 0.0:  # matrix is singular
            return False

        if np.abs(mu/np.round(mu) - 1) > eps:  # determinant is not integer
            return False
        
        # now check whether any element is non-integer
        for mij in m.ravel():
            if np.round(mij) == 0.0:
                if abs(mij)/np.sqrt(np.abs(mu)) > eps:
                    return False
            elif np.abs(mij/np.round(mij) - 1) > eps:
                return False
        return True
    
    def isRepresentable(self, m, basis):
        transform = np.dot(m, basis)
        nBasis = np.linalg.norm(basis, axis=1)
        nTransf = np.linalg.norm(transform, axis=1)
        g1 = nTransf[0]/nBasis[0]
        g2 = nTransf[1]/nBasis[1]
        mu = np.abs(np.linalg.det(m))
        
        return abs(mu/(g1*g2) - 1) < 1e-8
    
    def primitiveOrCentered(self, m, basis):
        invT = np.array([[1, -1], [1, 1]])
        
        primitive = self.isRepresentable(m, basis)
        ctrd = self.isRepresentable(np.dot(invT, m), basis)
        if primitive:
            prefix = 'p'
        elif ctrd:
            prefix = 'c'
            m = np.dot(invT, m)
        else:
            return False
        
        nBasis = np.linalg.norm(basis, axis=1)
        transform = np.dot(m, basis)
        nTransf = np.linalg.norm(transform, axis=1)
        g1 = nTransf[0]/nBasis[0]
        g2 = nTransf[1]/nBasis[1]
        cosAlpha = np.dot(transform[0], basis[0])/(nTransf[0]*nBasis[0])
        
        return (prefix, g1, g2, cosAlpha)
    
    def squareToProdOfSquares(self, number):
        # takes number, finds all prime factors, returns a tuple, first 
        # element is a product of all primes showing up an even number of 
        # times, the second one the rest. Useful to turn, e.g., sqrt(12) 
        # into 2sqrt(3)
        
        factors = list(self.primeFactors(number))
        if not factors:
            factors = [1]
        uniqueFact = sorted(list(set(factors)))
        countFact = [((factors.count(fac) // 2)*2, factors.count(fac) % 2)
                      for fac in uniqueFact]
        (pow2, restPow) = zip(*countFact)
        (squares, remainders) = zip(*[(fact**pow, fact**rem)
                                      for (fact, pow, rem)
                                      in zip(uniqueFact, pow2, restPow)
                                      ])
        
        return (np.prod(squares), np.prod(remainders))
    
    def primeFactors(self, n):
        f = 2
        increments = itertools.chain(
                                  [1, 2, 2],
                                  itertools.cycle([4, 2, 4, 2, 4, 6, 2, 6])
                                  )
        for incr in increments:
            if f*f > n:
                break
            while n % f == 0:
                yield f
                n //= f
            f += incr
        if n > 1:
            yield n
