# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

Class containing parameters read from the PARAMETERS file, and some parameters
defined at runtime. Most default values are defined here.
"""

import numpy as np
import logging
import os
import random
import shutil

try:
    import matplotlib.pyplot as plt
except:
    plotting = False
else:
    plotting = True

from tleedmlib.files.iodeltas import checkDelta
from tleedmlib.leedbase import getMaxTensorIndex
from tleedmlib.base import available_cpu_count

logger = logging.getLogger("tleedm.rparams")

###############################################
#                CLASSES                      #
###############################################

class SearchPar:
    """Stores properties of ONE parameter of the search, i.e. what variation
    of what atom is linked to this parameter."""
    def __init__(self, atom, mode, el, deltaname):
        self.atom = atom
        self.mode = mode
        self.el = el
        self.deltaname = deltaname
        self.steps = -1     # not used for interpretation, info only
        self.restrictTo = None  # None, Index, or other search par
        self.linkedTo = None    # other search par linked via 'atom number'
        self.parabolaFit = {"min": None, "curv": None} # parabolic fit
        d = {}
        if mode == "occ":
            self.steps = len(next(iter(atom.disp_occ.values())))
        else:
            if mode == "geo":
                d = atom.disp_geo
            elif mode == "vib":
                d = atom.disp_vib
        if len(d) > 0:
            if el in d:
                k = el
            else:
                k = "all"
            self.steps = len(d[k])

class DomainParameters:
    """Stores workdir, slab and runparams objects for each domain"""
    def __init__(self, workdir, homedir, name):
        self.workdir = workdir  # path do sub-directory for domain calculation
        self.homedir = homedir  # path to main tleedm working directory
        self.name = name        # domain name as defined by user
        self.sl = None
        self.rp = None
        
        self.refcalcRequired = False
        self.tensorDir = None

class Rparams:
    """Stores the parameters found in a PARAMETERS file (default values in
    __init__), as well as some parameters defined at runtime."""
    def __init__(self):
        self.readParams = {}    # original parameters as read from file
        
        self.ATTENUATION_EPS = 0.0001
        self.BULKDOUBLING_EPS = 0.001
        self.BULKDOUBLING_MAX = 10
        self.BULK_REPEAT = None
        self.DOMAINS = []       # list of domains (name, path)
        self.DOMAIN_STEP = 1   # area step in percent for domain search
        self.ELEMENT_MIX = {}   #if any ELEMENT_MIX is defined, it will be
                    #  added to the dictionary with the element name as the
                    #  label and the splitlist as the value
        self.ELEMENT_RENAME = {}  #defines what chemical element should be
                         #  used for a POSCAR element, if ELEMENT_MIX is not
                         #  defined for that element.
        self.FILAMENT_WF = 2.65  # work function of emitting cathode
        self.FORTRAN_COMP = ["",""] #before files, after files
        self.FORTRAN_COMP_MPI = ["",""] #before files, after files - for search
        self.GAUSSIAN_WIDTH = 0.5
        self.GAUSSIAN_WIDTH_SCALING = 0.5
        self.HALTING = 2    # 2: stop for user confirmation of intermediate
                # results. 1: stop for minor warnings. 0: always stop.
        self.IV_SHIFT_RANGE = [-3, 3, -1]  # step of -1: init from data
        self.LAYER_CUTS = ["dz(1.2)"] # strings either specifying a c or z 
                    # cutoff, or a list of c coordinates separating layers
        self.LAYER_STACK_VERTICAL = True
        self.LOG_DEBUG = False
        self.LOG_SEARCH = False
        self.N_BULK_LAYERS = 1           #number of bulk layers
        self.N_CORES = 0        # number of cores
        self.PHASESHIFT_EPS = 0  # defined in updateDerivedParams
        self.PHI = 0.0          #from BEAM_INCIDENCE
        self.PLOT_COLORS_RFACTOR = None
        self.RUN = [0,1,2,3]        #what segments should be run
        self.R_FACTOR_TYPE = 1  # which R-factor to use; 1: Pendry, 2: R2,
                                #   3: Zanazzi-Jona
        self.R_FACTOR_SMOOTH = 0
        self.SCREEN_APERTURE = 110.
        self.SEARCH_BEAMS = 0   # which beams to use for search. 0: average,
                                #   1: integer, 2: fractional
        self.SEARCH_CULL = 0.1  # kill off worst part of population during 
                                #   search, replace with others
        self.SEARCH_CULL_TYPE = "genetic"  # clone: copy survivor; genetic: 
                                # offspring through random mixing of two 
                                # survivors. Random: re-initialize random
        self.SEARCH_MAX_GEN = 50000 # maximum number of generations in search
        self.SEARCH_MAX_DGEN = {"all": 0, "best": 0, "dec": 100}
                # maximum number of generations without change before search 
                #   is stopped. All: all configs, best: only 1, dec: best 10%
                #    0: don't use parameter
        self.SEARCH_MAX_DGEN_SCALING = {"all": None, "best": None, "dec": None}
                                                        # derived
        self.SEARCH_LOOP = False    # repeat displacements blocks until no 
                                    #   further improvement?
        self.SEARCH_POPULATION = 0  # trial structures in search
        self.SEARCH_START = "crandom"
        self.SITE_DEF = {}      #labels are the element names, content is
                               #  dictionaries of format {sitename, list of
                               #  atom numbers in POSCAR}
        self.STOP = False
        self.SUPERLATTICE = np.identity(2, dtype=float)
        self.SUPPRESS_EXECUTION = False
        self.SYMMETRIZE_INPUT = True
        self.SYMMETRY_CELL_TRANSFORM = np.identity(2, dtype=float)
        self.SYMMETRY_EPS = 0.1
        self.SYMMETRY_EPS_Z = 0.1
        self.SYMMETRY_FIND_ORI = True
        self.SYMMETRY_FIX = ''
        self.TENSOR_INDEX = None  # default: pick highest in Tensors folder
        self.TENSOR_OUTPUT = []  #defines whether Tensor output is required for
                                 #  a layer. Default = 1.
        self.THEO_ENERGIES = [-1, -1, -1]  #default: [20, 800, 2], initialized
                                           #  in tleedm.py / runSection / INIT
        self.THETA = 0.0        #from BEAM_INCIDENCE
        self.TL_VERSION = 0.    # requested TensErLEED version
        self.T_EXPERIMENT = None
        self.T_DEBYE = None
        self.V0_IMAG = 4.5               # !!! CHOOSE BETTER DEFAULT?
        self.V0_REAL = "default"   # 'default' will  use rundgren format with
                                   #    values from PHASESHIFTS
        self.V0_Z_ONSET = 1.0
        self.VIBR_AMP_SCALE = []   # read as list of strings, interpret later

        # script-defined values
        self.LMAX = 0     #will be calculated based on PHASESHIFT_EPS parameter
        # self.inpot_re_set = {"total": False, "bulk": False, "surface": False,
        #                      "delta": False}
        # self.inpot_im_set = {"total": False, "bulk": False, "surface": False,
        #                      "delta": False}    #accept only two inputs, then
        #                                         #  ignore & warn

        # variable states
        self.workdir = os.getcwd() 
                                # MAIN WORK DIRECTORY; where to look for files
        self.searchConvInit = {"gaussian": None, 
                              "dgen": {"all": None, "best": None, "dec": None}}
        self.searchMaxGenInit = self.SEARCH_MAX_GEN
        # script progress tracking
        self.halt = 0
        self.systemName = ""
        self.timestamp = ""
        self.manifest = ["AUX","OUT"]
        self.fileLoaded = {"PARAMETERS": True, "POSCAR": False,
            "IVBEAMS": False, "VIBROCC": False, "PHASESHIFTS": False,
            "DISPLACEMENTS": False, "BEAMLIST": False, "EXPBEAMS": False}
        self.runHistory = []   #sections that have been executed before
        self.lastOldruns = []  # copy of runHistory when last oldruns folder
                               #  was created
        self.superlattice_defined = False
        self.ivbeams_sorted = False
        self.last_R = None
        self.stored_R = {"refcalc": None, "superpos": None}
        
        # domains
        self.domainParams = []
        self.pseudoSlab = None

        # data from files
        self.beamlist = []  # lines as strings from _BEAMLIST
        self.ivbeams = []   # uses Beam class; list of beams only
        self.expbeams = []  # uses Beam class; contains intensities
        self.theobeams = {"refcalc": [], "superpos": None} # uses Beam class; 
                                                        #  contains intensities
        self.phaseshifts = []
        self.phaseshifts_firstline = "" # contains parameters for MUFTIN
        self.refcalc_fdout = ""
        self.superpos_specout = ""
        self.disp_blocks = []    # tuples (lines, name) in DISPLACEMENTS file
        self.disp_block_read = False # current displacements block read?
        self.disp_loops = []       # list of tuples (loopStart, loopEnd)
        self.controlChemBackup = None

        # search parameters
        self.searchpars = []
        self.searchResultConfig = None
        self.search_atlist = []    # atoms that are relevant for the search
        self.search_maxfiles = 0   # maximum number of delta files for one atom
        self.search_maxconc = 1    # maximum number of concentration steps

        self.indyPars = 0        # number of independent parameters
        self.mncstep = 0    # MAX NUMBER OF VARIATIONS (geo. times therm.)
                            #   IN 1 FILE
        self.search_index = 0      # which DISPLACEMENTS block is being done
        
        # plotting data
        self.rfacscatter_all = []  # tuples (gen, r, size, color)
        self.rfacscatter = []  # same, but thinned out along gens
        self.parScatter = [[]]      # tuples (gen, mean scatter, max scatter) 
                                    #   per search
        self.searchplots = [("", [], [], [], [])]  # (name, gens, min, max, 
                                                   #    mean) for each search
        self.lastParScatterFigs = {}  # complete figures for each search, 
                                      #   with search names as keys

    def storeRfacScatter(self, x, y, s, c):
        """Adds a list of points for r-factor scatter plots to 
        self.rfacscatter"""
        spacing = x[-1]/50
        pg = x[-1]
        self.rfacscatter_all.extend(list(zip(x,y,s,c)))
        self.rfacscatter = []
        for p in self.rfacscatter_all[::-1]:
            if p[0] <= pg - spacing or p[0] == pg:
                self.rfacscatter.append(p)
                pg = p[0]
        
        
    def updateDerivedParams(self):
        """Checks which derivative parameters (which cannot be calculated at
        initialization) can be calculated now"""
        # TENSOR_INDEX:
        if self.TENSOR_INDEX is None:
            self.TENSOR_INDEX = getMaxTensorIndex()
        # SEARCH_CONVERGENCE:
        if self.searchConvInit["gaussian"] is None:
            self.searchConvInit["gaussian"] = self.GAUSSIAN_WIDTH
        for k in ["all", "best", "dec"]:
            if self.SEARCH_MAX_DGEN_SCALING[k] is None:
                self.SEARCH_MAX_DGEN_SCALING[k] = int(
                                              1 / self.GAUSSIAN_WIDTH_SCALING)
            if self.searchConvInit["dgen"][k] is None:
                self.searchConvInit["dgen"][k] = self.SEARCH_MAX_DGEN[k]
        # Phaseshifts-based:
        if self.fileLoaded["PHASESHIFTS"]:
            # get highest required energy index
            hi = len(self.phaseshifts)-1
            if self.THEO_ENERGIES[1] > 0:
                for i in range(0, len(self.phaseshifts)):
                    if self.phaseshifts[i][0]*27.2116 > self.THEO_ENERGIES[1]:
                        hi = i
                        break
            # LMAX
            if self.PHASESHIFT_EPS == 0:
                self.PHASESHIFT_EPS = 0.05
            if self.LMAX == 0:  #determine value from PHASESHIFT_EPS
                lmax = 1
                for el in self.phaseshifts[hi][1]:  # only check highest energy
                    for i, val in enumerate(el):
                        if abs(val) > self.PHASESHIFT_EPS and (i+1) > lmax:
                            lmax = i+1
                if lmax < 8:
                    logger.debug("Found small LMAX value based on "
                        "PHASESHIFT_EPS parameter (LMAX="+str(lmax)+").")
                if lmax > 15:
                    lmax = 15
                    logger.info("The LMAX found based on the PHASESHIFT_EPS "
                        "parameter is greater than 15, which is currently "
                        "not supported. LMAX was set to 15.")
                self.LMAX = lmax
            else:       # sanity check: are large values ignored?
                warn = False
                highval = 0
                for el in self.phaseshifts[hi][1]:   #highest energy
                    for i, val in enumerate(el):
                        if abs(val) > 0.1 and (i+1) > self.LMAX:
                            warn = True
                            highval = max(highval, abs(val))
                if warn:
                    logger.warning("The LMAX value defined in the PARAMETERS "
                        "file leads to large phaseshift values being ignored "
                        "(highest value: {}). Consider using a higher LMAX, "
                        "or defining LMAX indirectly via PHASESHIFT_EPS."
                        .format(highval))
                    self.setHaltingLevel(1)
            # V0_REAL
            if self.V0_REAL == "default":
                llist = self.phaseshifts_firstline.split()
                c = []
                try:
                    for i in range(0,4):
                        c.append(float(llist[i+1]))
                except:
                    logger.error("Could not read Muftin parameters from "
                                  "_PHASESHIFTS file.")
                    raise
                self.V0_REAL = ("workfn-max("+str(round(c[0],2))
                        +", (("+str(round(c[1],2))+")+("+str(round(c[2],2))
                        +")/sqrt(EEV+workfn+("+str(round(c[3],2))+"))))")
    
    def updateCores(self):
        # if N_CORES is undefined, tries to find it
        if self.N_CORES != 0:
            return 0
        try:
            self.N_CORES = available_cpu_count()
        except:
            logger.error("Failed to detect number of cores.")
        logger.info("Automatically detected number of available CPUs: {}"
                     .format(self.N_CORES))
        if self.N_CORES == 0:
            logger.error("Failed to detect number of cores.")
            raise RuntimeError("N_CORES undefined, automatic detection failed")
        return
        
    def resetSearchConv(self):
        """Sets the search convergence and tracking parameters back to their 
        initial values."""
        self.controlChemBackup = None
        self.disp_block_read = False
        self.rfacscatter_all = []
        self.searchplots.append(("", [], [], [], []))
        self.parScatter.append([])
        self.SEARCH_MAX_GEN = self.searchMaxGenInit
        if self.searchConvInit["gaussian"] is not None:
            self.GAUSSIAN_WIDTH = self.searchConvInit["gaussian"]
        for k in ["best", "all", "dec"]:
            if self.searchConvInit["dgen"][k] is not None:
                self.SEARCH_MAX_DGEN[k] = self.searchConvInit["dgen"][k]
        return
        
    def setHaltingLevel(self, set_to):
        """Sets the halting level self.halt to value set_to, if that value is 
        higher than the current halting level."""
        if set_to > self.halt:
            logger.debug("Raising halting level to {}".format(set_to))
            self.halt = set_to

    def initTheoEnergies(self):
        if not -1 in self.THEO_ENERGIES:
            return
        info = False
        if self.fileLoaded["EXPBEAMS"]:
            minen = min(self.expbeams[0].intens)
            maxen = max(self.expbeams[0].intens)
            for beam in self.expbeams:
                if min(beam.intens) < minen:
                    minen = min(beam.intens)
                if max(beam.intens) > maxen:
                    maxen = max(beam.intens)
            values = [minen - 3, maxen + 3, 3]
            if -1 in self.THEO_ENERGIES[:2]:
                info = True
        else:
            values = [20, 800, 3]
        for i in range(0,3):
            if self.THEO_ENERGIES[i] == -1:
                self.THEO_ENERGIES[i] = values[i]
        if info:
            logger.debug("Initialized energy range from experimental beams "
                  "file: {:.2f} to {:.2f} eV".format(self.THEO_ENERGIES[0], 
                                                     self.THEO_ENERGIES[1]))
        d = ((self.THEO_ENERGIES[1] - self.THEO_ENERGIES[0])
                                            % self.THEO_ENERGIES[2])
        if d != 0:
            self.THEO_ENERGIES[1] += (self.THEO_ENERGIES[2] - d)

    def getFortranComp(self, comp="auto"):
        """Checks whether ifort or gfortran are present, and sets FORTRAN_COMP
        accordingly. If the 'comp' parameter is set to 'auto', will check
        ifort first, then gfortran. If 'comp' is set to 'ifort' or 'gfortran',
        only that compiler will be checked. Returns None."""
        supported = ["ifort", "gfortran"]   #supported compilers, in order of
                                            #  priority
        found = ""
        if comp == "auto":
            check = supported
        elif comp in supported:
            check = [comp]
        else:
            logger.error("Rparams.getFortranComp: requested compiler is not "
                          "supported.")
            raise RuntimeError("Fortran compiler not supported.")
        for c in check:
            if shutil.which(c, os.X_OK) != None:
                found = c
                break
        if found == "":
            if comp == "auto":
                logger.error("Rparams.getFortranComp: No fortran compiler "
                              "found.")
            else:
                logger.error("Rparams.getFortranComp: Requested fortran "
                              "compiler not found.")
            raise RuntimeError("Fortran compiler not found.")
        if found == "ifort":
            self.FORTRAN_COMP = ["ifort -O2 -I/opt/intel/mkl/include",
                "-L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 "
                "-lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl"]
            logger.debug("Using fortran compiler: ifort")
        elif found == "gfortran":
            self.FORTRAN_COMP = ["gfortran -O2", "-llapack -lpthread"]
            logger.debug("Using fortran compiler: gfortran")
        return
    
    def getFortranMpiComp(self, comp="auto"):
        """Checks whether mpiifort or mpifort are present, and sets 
        FORTRAN_COMP mpi accordingly. If the 'comp' parameter is set to 'auto',
         will check mpifort first, then mpiifort. If 'comp' is set to 
         'mpiifort' or 'mpifort', only that compiler will be checked. Returns 
         None."""
        supported = ["mpiifort", "mpifort"]   #supported compilers, in order of
                                              #  priority
        found = ""
        if comp == "auto":
            check = supported
        elif comp in supported:
            check = [comp]
        else:
            logger.error("Rparams.getFortranMpiComp: requested compiler is "
                          "not supported.")
            raise RuntimeError("Fortran MPI compiler not supported.")
        for c in check:
            if shutil.which(c, os.X_OK) != None:
                found = c
                break
        if found == "":
            if comp == "auto":
                logger.error("Rparams.getFortranMpiComp: No fortran compiler "
                              "found.")
            else:
                logger.error("Rparams.getFortranMpiComp: Requested fortran "
                              "compiler not found.")
            raise RuntimeError("Fortran MPI compiler not found.")
        if found == "mpiifort":
            self.FORTRAN_COMP_MPI = ["mpiifort -Ofast", ""]
            logger.debug("Using fortran compiler: mpiifort")
        elif found == "mpifort":
            self.FORTRAN_COMP_MPI = ["mpifort -Ofast -no-pie", ""]
            logger.debug("Using fortran compiler: mpifort")
        return
    
    def renormalizeDomainParams(self, config):
        """Takes a list of parameter indices as produced by e.g. 
        getRandomConfig, and checks it for domain parameters. If the domain 
        parameters can be multiplied by an integer value and still be in 
        range, returns the list with multiplied values; else, returns the 
        unchanged list."""
        domain_indices = [i for (i, sp) in enumerate(self.searchpars) 
                          if sp.mode == "dom"]
        if not domain_indices or all([config[i] == 1 for i in domain_indices]):
            return config
        mult = 1
        domain_steps = self.searchpars[domain_indices[0]].steps
        while (sum([config[i]-1 for i in domain_indices]) * (mult+1) 
                                                           <= domain_steps-1):
            mult += 1
        for i in domain_indices:
            config[i] = ((config[i] - 1) * mult) + 1
        return config

    def getCenteredConfig(self):
        """Returns a list of 'centered' parameter indices, i.e. all in the 
        middle of their respective range"""
        return ([int((sp.steps + 1)/2) for sp in self.searchpars])
    
    def getPredictConfig(self, best_config=None, curv_cutoff = 1e-4):
        """Returns a list of parameter indices as determined by the parabola 
        fit, if a good fit was achieved for a given parameter; all other 
        parameters are cloned from the best configuration if passed, or 
        centered if not. curv_cutoff defines the minimal curvature of the 
        parabolas to be used."""
        l = []
        for (i,sp) in enumerate(self.searchpars):
            if (sp.parabolaFit["min"] is not None and 
                    sp.parabolaFit["curv"] is not None and
                    sp.parabolaFit["curv"] > curv_cutoff):
                l.append(int(round(sp.parabolaFit["min"])))
            else:
                if best_config is not None:
                    l.append(best_config[i])
                else:
                    l.append(int((sp.steps + 1)/2))
        return self.renormalizeDomainParams(l)
    
    def getRandomConfig(self):
        """Returns a list of 'random' parameter indices, but makes sure that 
        linked parameters have the same value"""
        l = [-1] * len(self.searchpars)
        for i, sp in enumerate(self.searchpars):
            l[i] = random.randint(1, sp.steps)
        for i, sp in enumerate(self.searchpars):
            if sp.linkedTo is not None:
                l[i] = l[self.searchpars.index(sp.linkedTo)]
            elif sp.restrictTo is not None:
                if type(sp.restrictTo) == int:
                    l[i] = sp.restrictTo
                else:
                    l[i] = l[self.searchpars.index(sp.restrictTo)]
        if -1 in l:
            logger.error("Rparams.getRandomConfig failed: {}".format(l))
            return []
        return self.renormalizeDomainParams(l)
    
    def getOffspringConfig(self, parents):
        """Returns a list of parameter indices generated as a random mix of 
        the parameters from two of the configurations passed as 'parents', 
        picked at random if there are more than two. If possible, the parents 
        will be picked such that they are not identical, but the offspring may 
        nevertheless be identical with one of the parents."""
        parents = parents[:]
        if len(parents) < 2:
            if len(parents) == 0:
                logger.warning("Rparams.getOffspringConfig: Cannot create "
                    "offspring configuration without parents. Returning "
                    "random configuration")
                return self.getRandomConfig()
            else:
                logger.warning("Rparams.getOffspringConfig: Only one parent "
                    "passed. Returning clone.")
                return parents[0]
        i = random.randint(0, len(parents)-1)
        p2 = [parents.pop(i)]
        while len(p2) == 1:
            if len(parents) == 1:
                p2.append(parents[0])
                break
            i = random.randint(0, len(parents)-1)
            p = parents.pop(i)
            if p2[0] != p:
                p2.append(p)
        l = [-1] * len(self.searchpars)
        for i, sp in enumerate(self.searchpars):
            if sp.restrictTo is None and sp.linkedTo is None:
                l[i] = p2[random.randint(0,1)][i]
        for i, sp in enumerate(self.searchpars):
            if sp.restrictTo is not None:
                if type(sp.restrictTo) == int:
                    l[i] = sp.restrictTo
                else:
                    l[i] = l[self.searchpars.index(sp.restrictTo)+1]
            elif sp.linkedTo is not None:
                l[i] = l[self.searchpars.index(sp.linkedTo)+1]
        if -1 in l:
            logger.error("Rparams.getOffspringConfig failed: {}".format(l))
            return []
        return self.renormalizeDomainParams(l)
    
    def closePdfReportFigs(self):
        global plotting
        if not plotting:
            return
        
        for searchname in self.lastParScatterFigs:
            for f in self.lastParScatterFigs[searchname]:
                try:
                    plt.close(f)
                except:
                    pass

    def generateSearchPars(self, sl, subdomain=False):
        """Initializes a list of searchpar objects, and assigns delta files to
        atoms if required."""
        if self.domainParams:
            return(self.generateSearchPars_domains(sl))
        self.searchpars = []
        self.search_maxfiles = 0   # maximum number of delta files for one atom
        self.search_maxconc = 1
        self.indyPars = 0        # number of independent parameters
        self.mncstep = 0     # MAX NUMBER OF VARIATIONS (geo. times therm.)
                             #   IN 1 FILE
        eqlist = []     # track which atoms are symmetry-linked to the ones already
                        #   done to not double-count indyPars
        # get list of atoms that appear in the search
        if (2 in self.runHistory or 42 in self.runHistory 
                                 or sl.deltasInitialized):
            # if delta has been run, information what deltas exist is stored
            atlist = [at for at in sl.atlist if (not at.layer.isBulk and
                                                 len(at.deltasGenerated) > 0)]
        else:
            logger.debug("Delta-amplitudes were not calculated in current "
                          "run; looking for delta files by name.")
            # see what delta files are present
            filelist = [filename for filename in os.listdir(".") if
                        filename.startswith('DEL_')]
            delN = []
            for filename in filelist:
                try:
                    delN.append(int(filename.split('_')[1]))
                except:
                    filelist.remove(filename)
            atlist = [at for at in sl.atlist if (not at.layer.isBulk and
                                                 at.oriN in delN)]
            for at in atlist:
                deltaCandidates = [fn for fn in filelist
                                      if int(fn.split('_')[1]) == at.oriN]
                checkEls = list(at.disp_occ.keys())
                # check for vacancy:
                occlists = []
                for k in at.disp_occ:
                    occlists.append(at.disp_occ[k])
                for i in range(0,len(occlists[0])):
                    totalocc = 0.
                    for l in occlists:
                        if len(l) <= i:
                            logger.error("Inconsistent occupancy lists for "
                                          "atom "+str(at.oriN))
                            raise ValueError("Inconsistent input")
                        else:
                            totalocc += l[i]
                    if totalocc < 1 - 1e-4:
                        checkEls.append("vac")
                        break
                # now check if all deltas are there:
                for el in checkEls:
                    found = False
                    for df in [f for f in deltaCandidates 
                                               if f.split("_")[2] == el]:
                        if checkDelta(df, at, el, self):
                            found = True
                            at.deltasGenerated.append(df)
                            break
                    if not found:
                        logger.error("No appropriate Delta file found for "
                                      "{}, element {}".format(at, el))
                        raise RuntimeError("Missing Delta file")
            # sanity check: are displacements defined but deltas missing?
            for at in sl.atlist:
                # check whether at has non-trivial displacements:
                found = False
                for d in [at.disp_occ, at.disp_geo, at.disp_vib]:
                    for el in d:
                        if len(d[el]) > 1:
                            found = True
                            break
                if not found:
                    for el in at.disp_vib:
                        if at.disp_vib[el][0] != 0.0:
                            found = True
                            break
                if not found:
                    for el in at.disp_geo:
                        if np.linalg.norm(at.disp_geo[el][0]) > 0.0:
                            found = True
                            break
                if not found:
                    occlists = []
                    for k in at.disp_occ:
                        occlists.append(at.disp_occ[k])
                    for i in range(0,len(occlists[0])):
                        totalocc = 0.
                        for l in occlists:
                            if len(l) <= i:
                                break #error - will pop up again later...
                            else:
                                totalocc += l[i]
                        if totalocc < 1 - 1e-4:
                            found = True
                            break
                if found and at not in atlist:
                    logger.error("Atom {} has displacements defined, but no "
                                  "delta file was found! Run Delta-Amplitudes."
                                  .format(at.oriN))
                    raise RuntimeError("Delta file not found")
                elif not found and at in atlist:
                    atlist.remove(at) # delta file is there, but no
                                      #   displacements
            sl.deltasInitialized = True
        # sort atlist by displists
        al2 = []
        for dl in sl.displists:
            al2.extend([a for a in atlist if a in dl])
        al2.extend([a for a in atlist if a not in al2])
        atlist = al2
        md = {"geo": 1, "vib": 2, "occ": 3}
        splToRestrict = []
        for at in atlist:
            if len(at.deltasGenerated) > self.search_maxfiles:
                self.search_maxfiles = len(at.deltasGenerated)
            for fn in at.deltasGenerated:
                el = fn.split("_")[2]
                if el == "vac":
                    self.searchpars.append(SearchPar(at, "geo", "vac", fn))
                    continue
                mult = 1
                pars = 0
                for (mode, d) in [("vib", at.disp_vib),
                                  ("geo", at.disp_geo)]:
                    if el in d:
                        k = el
                    else:
                        k = "all"
                    if len(d[k]) > 1 or (len(d[k]) == 1 and 
                                         ((mode == "geo" and 
                                             np.linalg.norm(d[k][0]) > 0.)
                                          or (mode == "vib" and 
                                                      d[k][0] != 0.))):
                        pars += 1
                        sp = SearchPar(at, mode, el, fn)
                        self.searchpars.append(sp)
                        if el in at.constraints[md[mode]]:
                            k2 = el
                        else:
                            k2 = "all"
                        if k2 in at.constraints[md[mode]]:
                            c = at.constraints[md[mode]][k2]
                            if type(c) == int:
                                sp.restrictTo = c + 1
                            else:
                                splToRestrict.append((sp, c))
                        elif len(d[k]) > 1 and not at in eqlist:
                            self.indyPars += 1
                        if len(d[k]) > 1 and at in eqlist:
                            spl = [s for s in self.searchpars 
                                   if at in s.atom.displist 
                                   and s.mode == sp.mode 
                                   and (s.el == el 
                                        or el in ["", "all"])]
                            if spl:
                                sp.linkedTo = spl[0]
                    mult *= len(d[k])
                if pars == 0:
                    self.searchpars.append(SearchPar(at, "geo", el, fn))
                if mult > self.mncstep:
                    self.mncstep = mult
            sp = SearchPar(at, "occ", "", fn)
            self.searchpars.append(sp)
            occsteps = len(next(iter(at.disp_occ.values())))
            if occsteps > 1:
                if occsteps > self.search_maxconc:
                    self.search_maxconc = occsteps
                if at.constraints[3]:
                    c = list(at.constraints[3].values())[0]
                    if type(c) == int:
                        sp.restrictTo = c + 1
                    else:
                        splToRestrict.append((sp, c))
                else:
                    if not at in eqlist:  # occupation will actually vary
                        self.indyPars += 1  
                if at in eqlist:
                    spl = [s for s in self.searchpars if at 
                           in s.atom.displist and s.mode == 3]
                    if spl:
                        sp.linkedTo = spl[0]
            eqlist.extend(at.displist)  # do not consider those for future
                                        #    indyPars
        # if self.indyPars == 0:
        #     self.indyPars = 1
        splTargets = set()
        for (sp, (at, el)) in splToRestrict:
            # ind = self.searchpars.index(sp)
            spl = [s for s in self.searchpars if s != sp and s.atom == at and 
                   s.mode == sp.mode and (s.el == el or 
                                          el in ["", "all"] or
                                          sp.mode == 3)]
            if spl:
                sp.restrictTo = spl[0]
                splTargets.add(spl[0])
        for (sp, (at, el)) in [tup for tup in splToRestrict 
                                           if tup[0].restrictTo is None 
                                           and tup[0] not in splTargets]:
            logger.warning("Restricting search parameter for atom {}, "
                "element {}, mode {} failed: Could not identify target "
                "search parameter (atom {}, element {})."
                .format(sp.atom.oriN, sp.el, sp.mode, at.oriN, el))
        for (i, sp) in enumerate(self.searchpars):
            # restrict to lowest number index, resolve conflicts
            if sp.restrictTo not in self.searchpars:
                continue
            sp2 = None
            while sp2 != sp.restrictTo:
                sp2 = sp.restrictTo
                ind = self.searchpars.index(sp2)
                if (sp2.linkedTo in self.searchpars and 
                          self.searchpars.index(sp2.linkedTo) < ind):
                    sp.restrictTo = sp2.linkedTo
                elif (sp2.restrictTo in self.searchpars and 
                          self.searchpars.index(sp2.restrictTo) < ind):
                    sp.restrictTo = sp2.restrictTo
            if self.searchpars.index(sp.restrictTo) >= i:
                if sp.restrictTo.restrictTo is None:
                    sp.restrictTo.restrictTo = sp   # invert direction
                if type(sp.restrictTo.restrictTo) == int:
                    sp.restrictTo = sp.restrictTo.restrictTo
                else:
                    sp.restrictTo = None  # remove references to higher indices
        self.search_atlist = atlist
        if not subdomain:
            self.searchpars.append(SearchPar(None, "dom", "", ""))
            self.searchpars[-1].steps = 2
        return

    def generateSearchPars_domains(self, sl):
        """Runs generateSearchPars for every domain, then collates results."""
        self.searchpars = []
        self.indyPars = len(self.domainParams) - 1
        home = os.getcwd()
        for dp in self.domainParams:
            try:
                os.chdir(dp.workdir)
                r = dp.rp.generateSearchPars(dp.sl, subdomain=True)
            except:
                logger.error("Error while creating delta input for domain {}"
                             .format(dp.name))
                raise
            finally:
                os.chdir(home)
            if r != 0:
                logger.error("Error getting search parameters for domain {}: "
                             "{}".format(dp.name, r))
                raise RuntimeError("Error getting search parameters")
            for sp in [sp for sp in dp.rp.searchpars 
                       if type(sp.restrictTo) == int]:
                sp.restrictTo += len(self.searchpars)
            self.searchpars.extend(dp.rp.searchpars)
            self.indyPars += dp.rp.indyPars
        for dp in self.domainParams:
            self.searchpars.append(SearchPar(None, "dom", "", ""))
            self.searchpars[-1].steps = int(100 / self.DOMAIN_STEP) + 1
        self.search_maxfiles = max([dp.rp.search_maxfiles 
                                    for dp in self.domainParams])
        self.search_maxconc = max([dp.rp.search_maxconc 
                                    for dp in self.domainParams])
        self.mncstep = max([dp.rp.mncstep for dp in self.domainParams])
        return