# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

Class containing parameters read from the PARAMETERS file, and some parameters
defined at runtime. Most default values are defined here.
"""

import logging
import os
from pathlib import Path
import random
import shutil
from timeit import default_timer as timer

import numpy as np

try:
    import matplotlib.pyplot as plt
except Exception:
    plotting = False
else:
    plotting = True

from viperleed.tleedmlib import leedbase
from viperleed.tleedmlib.base import available_cpu_count
from viperleed.tleedmlib.checksums import (KNOWN_TL_VERSIONS,
                                           UnknownTensErLEEDVersionError)
from viperleed.tleedmlib.files.iodeltas import checkDelta

logger = logging.getLogger("tleedm.rparams")

###############################################
#                CLASSES                      #
###############################################


# Notice that the defaults in here that may be mutated during execution
# are saved as immutable types to prevent inadvertent modification of
# this global, and are rather converted to their mutable equivalent
# in the relevant places
DEFAULTS = {
    'EXPBEAMS_INPUT_FILE' : None,
    'PHASESHIFT_EPS': {
        'r': 0.1,
        'n': 0.05,
        'f': 0.01, # this is the default if nothing is given
        'e': 0.001,
    },
    'ZIP_COMPRESSION_LEVEL': 2,
    'SEARCH_EVAL_TIME':  60, # time interval between reads of SD.TL, TODO: should be dynamic?
    'RUN': (0, 1, 2, 3),
}


# parameter limits
# either tuple of (min, max) or list of allowed values
PARAM_LIMITS = {
    'LMAX': (1, 18),
    'INTPOL_DEG': ['3', '5'],
    }


class SearchPar:
    """Stores properties of ONE parameter of the search, i.e. what variation
    of what atom is linked to this parameter."""

    def __init__(self, atom, mode, el, deltaname):
        self.atom = atom
        self.mode = mode
        self.el = el
        self.deltaname = deltaname
        self.steps = 1
        self.edges = (None, None)  # the first and last value in the range
        self.center = 1  # the index closest to "no change" (Fortran index starting at 1)
        self.non_zero = False   # whether the center is truly "unchanged"
        self.restrictTo = None  # None, Index, or other search par
        self.linkedTo = None    # other search par linked via 'atom number'
        self.parabolaFit = {"min": None,
                            "err_co": np.nan, "err_unco": np.nan}
        d = {}
        if mode == "occ":
            el = next(iter(atom.disp_occ.keys()))  # look at any element
            self.steps = len(atom.disp_occ[el])
            self.center = atom.disp_center_index[mode][el] + 1 # (Fortran index starting at 1)
            self.non_zero = (abs(atom.disp_occ[el][self.center-1]
                                 - atom.site.occ[el]) >= 1e-4)
            edges = []
            for ind in (0, -1):
                edges.append(" + ".join("{:.2f} {}".format(
                    atom.disp_occ[e][ind], e) for e in atom.disp_occ
                    if atom.disp_occ[e][ind] > 0.005))
                if edges[-1] == "":
                    edges[-1] = "vac"
            self.edges = tuple(edges)
        else:
            if mode == "geo":
                d = atom.disp_geo
            elif mode == "vib":
                d = atom.disp_vib
        if len(d) > 0 and el != "vac":  # if vac: use defaults
            if el in d:
                k = el
            else:
                k = "all"
            self.steps = len(d[k])
            self.edges = (d[k][0], d[k][-1])
            if k not in atom.disp_center_index[mode]:
                self.center = atom.disp_center_index[mode]["all"] + 1
            else:
                self.center = atom.disp_center_index[mode][k] + 1
            self.non_zero = (np.linalg.norm(d[k][self.center-1]) >= 1e-4)


class DomainParameters:
    """Stores workdir, slab and runparams objects for each domain"""

    def __init__(self, workdir, homedir, name):
        self.workdir = Path(workdir)  # path to sub-directory for domain calculation
        self.homedir = Path(homedir)  # path to main tleedm working directory
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

        # FROM PARAMETERS FILE
        self.ATTENUATION_EPS = 0.001
        self.AVERAGE_BEAMS = None
        self.BULKDOUBLING_EPS = 0.001
        self.BULKDOUBLING_MAX = 10
        self.BULK_LIKE_BELOW = 0.
        self.BULK_REPEAT = None
        self.DOMAINS = []         # list of domains (name, path)
        self.DOMAIN_STEP = 1      # area step in percent for domain search
        self.ELEMENT_MIX = {}     # {element_name: splitlist}
        self.ELEMENT_RENAME = {}  # {element_name: chemical_element}
        self.EXPBEAMS_INPUT_FILE = DEFAULTS["EXPBEAMS_INPUT_FILE"]
        self.FILAMENT_WF = 2.65   # work function of emitting cathode
        self.FORTRAN_COMP = ["", ""]      # before files, after files
        self.FORTRAN_COMP_MPI = ["", ""]  # before files, after files
        self.GAUSSIAN_WIDTH = 0.5
        self.GAUSSIAN_WIDTH_SCALING = 0.5
        self.HALTING = 2    # 2: major concerns, 1: minor warnings, 0: always
        self.INTPOL_DEG = 3 # Degree of interpolation spline used in R-factor calculation
        self.IV_SHIFT_RANGE = [-3, 3, -1]  # step of -1: init from data
        self.LAYER_CUTS = ["dz(1.2)"]  # list of either str or c coordinates
        self.LAYER_STACK_VERTICAL = True
        self.LMAX = [0, 0]    # minimum and maximum LMAX
        self.LOG_DEBUG = False
        self.LOG_SEARCH = True
        self.N_BULK_LAYERS = 1           # number of bulk layers
        self.N_CORES = 0                 # number of cores
        self.OPTIMIZE = {"which": "none", "step": 0., "minpoints": 4,
                         "maxpoints": 10, "convergence": 0., "maxstep": 0.}
                    # settings for fd optimization
        self.PARABOLA_FIT = {"type": "none", "alpha": 1e-2, "mincurv": 5e-3,
                             "localize": 0}
        self.PHASESHIFT_EPS = DEFAULTS['PHASESHIFT_EPS']['f'] # changed in updateDerivedParams
        self.PHASESHIFTS_CALC_OLD = True # use old EEASiSSS version # TODO: once established, set to False or remove
        self.PHASESHIFTS_OUT_OLD = True  # output old PHASESHIFTS file # TODO: once established, set to False or remove
        self.PHI = 0.0           # from BEAM_INCIDENCE
        self.PLOT_IV = {'plot': True, 'axes': 'all', 'colors': [],
                        'legend': 'all', 'overbar': False, 'perpage': 2}
        self.RUN = self.get_default('RUN')        # what segments should be run
        self.R_FACTOR_LEGACY = True # use old runtime-compiled R-factor calculation
        self.R_FACTOR_TYPE = 1  # 1: Pendry, 2: R2, 3: Zanazzi-Jona
        self.R_FACTOR_SMOOTH = 0
        self.S_OVL = 0.3 # Muffin tin overlap parameter after Rundgren 2021, default is 0.3 - set or optimize in FD
        self.SCREEN_APERTURE = 110.
        self.SEARCH_BEAMS = 0   # 0: average, 1: integer, 2: fractional
        self.SEARCH_CULL = 0.1  # fraction of all, or absolute N if >1
        self.SEARCH_CULL_TYPE = "genetic"  # clone, genetic, random
        self.SEARCH_MAX_GEN = 100000  # maximum number of generations in search
        self.SEARCH_MAX_DGEN = {"all": 0, "best": 0, "dec": 100}
        # maximum number of generations without change before search
        #   is stopped. All: all configs, best: only 1, dec: best 10%
        #    0: don't use parameter
        self.SEARCH_MAX_DGEN_SCALING = {"all": None, "best": None, "dec": None}
        # will be defined in updateDerivedParams
        self.SEARCH_LOOP = False
        self.SEARCH_POPULATION = 0  # trial structures in search
        self.SEARCH_START = "crandom"
        self.SITE_DEF = {}   # {element_name: {sitename, [atom.oriN]}}
        self.STOP = False
        self.SUPERLATTICE = np.identity(2, dtype=float)
        self.SUPPRESS_EXECUTION = False
        self.SYMMETRIZE_INPUT = True
        self.SYMMETRY_CELL_TRANSFORM = np.identity(2, dtype=float)
        self.SYMMETRY_EPS = 0.1
        self.SYMMETRY_EPS_Z = 0.1
        self.SYMMETRY_FIND_ORI = True
        self.SYMMETRY_FIX = ''
        self.SYMMETRY_BULK = {}   # keys: group, rotation, mirror
        self.TENSOR_INDEX = None  # default: pick highest in Tensors folder
        self.TENSOR_OUTPUT = []  # per layer: write Tensor output? (0/1)
        self.THEO_ENERGIES = [-1, -1, -1]
        # default: [20, 800, 2], initialized in section INIT
        self.THETA = 0.0        # from BEAM_INCIDENCE
        self.TL_IGNORE_CHECKSUM = True
        self.TL_VERSION = 0.    # requested TensErLEED version
        self.TL_VERSION_STR = None  # TODO: replace with Version class once available
        self.T_EXPERIMENT = None
        self.T_DEBYE = None
        self.V0_IMAG = 4.5
        self.V0_REAL = "default"   # 'default' will read from PHASESHIFTS
        self.V0_Z_ONSET = 1.0
        self.VIBR_AMP_SCALE = []   # read as list of strings, interpret later
        self.ZIP_COMPRESSION_LEVEL = DEFAULTS['ZIP_COMPRESSION_LEVEL']

        # RUN VARIABLES
        self.starttime = timer()
        self.sourcedir = os.getcwd()  # where to find 'tensorleed'
        self.workdir = Path(os.getcwd())  # MAIN WORK DIRECTORY; where to find input
        self.compile_logs_dir = None
        self.searchConvInit = {
            "gaussian": None, "dgen": {"all": None, "best": None, "dec": None}}
        self.searchEvalTime = DEFAULTS['SEARCH_EVAL_TIME']  # time interval for reading SD.TL
        self.output_interval = None # changed in updateDerivedParams
        self.searchMaxGenInit = self.SEARCH_MAX_GEN
        self.searchStartInit = None
        # script progress tracking
        self.halt = 0
        self.systemName = ""
        self.timestamp = ""
        self.manifest = ["SUPP", "OUT"]
        self.fileLoaded = {
            "PARAMETERS": True, "POSCAR": False,
            "IVBEAMS": False, "VIBROCC": False, "PHASESHIFTS": False,
            "DISPLACEMENTS": False, "BEAMLIST": False, "EXPBEAMS": False}
        self.runHistory = []   # sections that have been executed before
        self.lastOldruns = []
        # copy of runHistory when last oldruns folder was created
        self.superlattice_defined = False
        self.ivbeams_sorted = False
        self.last_R = None
        self.stored_R = {"refcalc": None, "superpos": None}
        self.checklist = []  # output strings of things to check at program end

        # domains
        self.domainParams = []
        self.pseudoSlab = None

        # data from files
        self.beamlist = []  # lines as strings from BEAMLIST
        self.ivbeams = []   # uses Beam class; list of beams only
        self.expbeams = []  # uses Beam class; contains intensities
        self.theobeams = {"refcalc": [], "superpos": None}
        # uses Beam class; contains intensities
        self.phaseshifts = []
        self.phaseshifts_firstline = ""  # contains parameters for MUFTIN
        self.refcalc_fdout = ""
        self.superpos_specout = ""
        self.best_v0r = None     # best value for v0r from previous R-factor
        self.disp_blocks = []    # tuples (lines, name) in DISPLACEMENTS file
        self.disp_block_read = False  # current displacements block read?
        self.disp_loops = []          # list of tuples (loopStart, loopEnd)
        self.controlChemBackup = None

        # search parameters
        self.searchpars = []
        self.searchResultConfig = None
        self.search_atlist = []    # atoms that are relevant for the search
        self.search_maxfiles = 0   # maximum number of delta files for one atom
        self.search_maxconc = 1    # maximum number of concentration steps
        self.indyPars = 0        # number of independent parameters
        self.mncstep = 0    # max. steps (geo. times therm.) for one atom
        self.search_index = 0      # which DISPLACEMENTS block is being done

        # plotting data
        self.rfacscatter_all = []   # tuples (gen, r, size, color)
        self.rfacscatter = []       # same, but thinned out along gens
        self.parScatter = [[]]
        # tuples (gen, mean scatter, max scatter) per search
        self.searchplots = [("", [], [], [], [])]
        # (name, gens, min, max, mean) for each search
        self.lastParScatterFigs = {}
        # complete figures for each search, with search names as keys

    def get_default(self, param):
        """Return the default value of param."""
        value = DEFAULTS[param]
        if isinstance(value, tuple):
            value = list(value)
        return value

    def get_limits(self, param):
        """Return the smallest and largest acceptable values of param."""
        return PARAM_LIMITS[param]

    def total_energy_range(self):
        """Return the total overlapping energy range of experiment and
        theory. Note that this may change if experimental beams are dropped."""
        if not self.expbeams:
            return 0.
        expEnergies = []
        totalrange = 0.
        for b in self.expbeams:
            expEnergies.extend([k for k in b.intens if k not in expEnergies])
            totalrange += (min(max(b.intens), self.THEO_ENERGIES[1])
                           - max(min(b.intens), self.THEO_ENERGIES[0]))
        return totalrange

    def storeRfacScatter(self, x, y, s, c):
        """
        Adds a list of points for r-factor scatter plots to
        self.rfacscatter

        Parameters
        ----------
        x, y, s, c : lists of floats
            coordinates, size and color of points.

        Returns
        -------
        None.

        """
        spacing = x[-1]/50
        pg = x[-1]
        self.rfacscatter_all.extend(list(zip(x, y, s, c)))
        self.rfacscatter = []
        for p in self.rfacscatter_all[::-1]:
            if p[0] <= pg - spacing or p[0] == pg:
                self.rfacscatter.append(p)
                pg = p[0]

    def updateDerivedParams(self):
        """
        Checks which derivative parameters (which cannot be calculated at
        initialization) can be calculated now

        Returns
        -------
        None.

        """
        # TENSOR_INDEX:
        if self.TENSOR_INDEX is None:
            self.TENSOR_INDEX = leedbase.getMaxTensorIndex()
        # TL_VERSION:
        if self.TL_VERSION == 0.:
            path = os.path.join(self.sourcedir, "tensorleed")
            ls = [dn for dn in os.listdir(path)
                  if (os.path.isdir(os.path.join(path, dn))
                      and dn.startswith("TensErLEED"))]
            highest = 0.0
            namestr = ""
            for dn in ls:
                try:
                    s = dn.split('v')[-1]
                    f = float(s)
                    if f > highest:
                        highest = f
                        namestr = s
                except Exception:
                    pass
            self.TL_VERSION = highest
            if highest > 0.:
                logger.debug("Detected TensErLEED version " + namestr)

        # TL_VERSION_STR
        # try simple conversion to string
        self.TL_VERSION_STR = f"{self.TL_VERSION:.2f}"
        if self.TL_VERSION_STR not in KNOWN_TL_VERSIONS:
            # try again without trailing zero
            if self.TL_VERSION_STR.endswith('0'):
                self.TL_VERSION_STR = self.TL_VERSION_STR[:-1]
        if (self.TL_VERSION_STR not in KNOWN_TL_VERSIONS 
                and not self.TL_IGNORE_CHECKSUM):
            raise UnknownTensErLEEDVersionError(
                f"Unrecognized TensErLEED version: {self.TL_VERSION_STR}. "
                "Consider editing KNOWN_TL_VERSIONS global in checksums.py "
                "or setting TL_IGNORE_CHECKSUM = True."
            )

        # SEARCH_CONVERGENCE:
        if self.searchConvInit["gaussian"] is None:
            self.searchConvInit["gaussian"] = self.GAUSSIAN_WIDTH
        for k in ["all", "best", "dec"]:
            if self.SEARCH_MAX_DGEN_SCALING[k] is None:
                self.SEARCH_MAX_DGEN_SCALING[k] = int(
                                              1 / self.GAUSSIAN_WIDTH_SCALING)
            if self.searchConvInit["dgen"][k] is None:
                self.searchConvInit["dgen"][k] = self.SEARCH_MAX_DGEN[k]
        if self.output_interval is None:
            # set output interval to SEARCH_CONVERGENCE dgen, but cap at 1000
            use_dgen = min(dgen for dgen in self.searchConvInit["dgen"].values() if dgen > 0) or 1000
            self.output_interval = int(min(use_dgen or 1000, 1000))  # default to 1000 if all dgen are 0 (default)
        if self.searchStartInit is None:
            self.searchStartInit = self.SEARCH_START
        # Phaseshifts-based:
        if self.fileLoaded["PHASESHIFTS"]:
            # get highest required energy index
            hi = len(self.phaseshifts)-1
            if self.THEO_ENERGIES[1] > 0:
                for i in range(0, len(self.phaseshifts)):
                    if (self.phaseshifts[i][0]*27.211396
                            > self.THEO_ENERGIES[1]):
                        hi = i
                        break
            # LMAX
            min_set = True
            if self.PHASESHIFT_EPS == 0:
                self.PHASESHIFT_EPS = DEFAULTS['PHASESHIFT_EPS']['f']
            if self.LMAX[0] <= 0:
                self.LMAX[0] = 6
                min_set = False
            if self.LMAX[1] == 0:  # determine value from PHASESHIFT_EPS
                lmax = 1
                for el in self.phaseshifts[hi][1]:  # only check highest energy
                    for i, val in enumerate(el):
                        if abs(val) > self.PHASESHIFT_EPS and (i+1) > lmax:
                            lmax = i+1
                if lmax < 8 and not min_set:
                    logger.debug(
                        "Found small LMAX value based on "
                        "PHASESHIFT_EPS parameter (LMAX="+str(lmax)+").")
                if lmax > 18:
                    lmax = 18
                    logger.info(
                        "The LMAX found based on the PHASESHIFT_EPS "
                        "parameter is greater than 18, which is currently "
                        "not supported. LMAX was set to 18.")
                self.LMAX[1] = lmax
            else:       # sanity check: are large values ignored?
                warn = False
                highval = 0
                for el in self.phaseshifts[hi][1]:   # highest energy
                    for i, val in enumerate(el):
                        if abs(val) > 0.1 and (i+1) > self.LMAX[1]:
                            warn = True
                            highval = max(highval, abs(val))
                if warn:
                    logger.warning(
                        "The LMAX value defined in the PARAMETERS "
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
                    for i in range(0, 4):
                        c.append(float(llist[i + 1]))
                except ValueError:
                    logger.error("Could not read Muftin parameters from "
                                 "PHASESHIFTS file.")
                    raise
                self.V0_REAL = c

    def updateCores(self):
        """
        If self.N_CORES is undefined, tries to find it

        Raises
        ------
        RuntimeError
            If unable to detect number of cores.

        Returns
        -------
        None

        """
        if self.N_CORES != 0:
            return
        try:
            self.N_CORES = available_cpu_count()
        except Exception:
            logger.error("Failed to detect number of cores.")
            raise
        else:
            logger.info("Automatically detected number of available CPUs: {}"
                        .format(self.N_CORES))
        if self.N_CORES == 0:
            logger.error("Failed to detect number of cores.")
            raise RuntimeError("N_CORES undefined, automatic detection failed")
        return

    def resetSearchConv(self):
        """
        Resets the search convergence and tracking parameters to their
        initial values.

        Returns
        -------
        None.

        """
        self.controlChemBackup = None
        self.disp_block_read = False
        self.rfacscatter_all = []
        self.searchplots.append(("", [], [], [], []))
        self.parScatter.append([])
        self.SEARCH_MAX_GEN = self.searchMaxGenInit
        if self.searchConvInit["gaussian"] is not None:
            self.GAUSSIAN_WIDTH = self.searchConvInit["gaussian"]
        if self.searchStartInit is not None:
            self.SEARCH_START = self.searchStartInit
        if self.SEARCH_START == "control":
            self.SEARCH_START = "crandom"
        for k in ["best", "all", "dec"]:
            if self.searchConvInit["dgen"][k] is not None:
                self.SEARCH_MAX_DGEN[k] = self.searchConvInit["dgen"][k]
        return

    def setHaltingLevel(self, set_to):
        """
        Sets the halting level self.halt to value set_to, if that value is
        higher than the current halting level. Outputs debug message.

        Parameters
        ----------
        set_to : int
            New value for self.halt

        Returns
        -------
        None.

        """
        if set_to > self.halt:
            logger.debug("Raising halting level to {}".format(set_to))
            self.halt = set_to
        return

    def initTheoEnergies(self):
        """
        Initializes values for the THEO_ENERGIES parameter either to default,
        or if experimental beams were loaded to values corresponding to the
        experimental energy range.

        Returns
        -------
        None.

        """
        if -1 not in self.THEO_ENERGIES:
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
        for i in range(0, 3):
            if self.THEO_ENERGIES[i] == -1:
                self.THEO_ENERGIES[i] = values[i]
        if info:
            logger.debug(
                "Initialized energy range from experimental beams file: "
                "{:.2f} to {:.2f} eV".format(self.THEO_ENERGIES[0],
                                             self.THEO_ENERGIES[1]))
        d = ((self.THEO_ENERGIES[1] - self.THEO_ENERGIES[0])
             % self.THEO_ENERGIES[2])
        if d != 0:
            self.THEO_ENERGIES[1] += (self.THEO_ENERGIES[2] - d)
        return

    # TODO: eventually, these default values should be moved to some constant or other file
    def getFortranComp(self, comp="auto", skip_check=False):
        """
        Checks whether ifort or gfortran are present, and sets FORTRAN_COMP
        accordingly.

        Parameters
        ----------
        comp : str, optional
            'auto' (default): will check ifort first, then gfortran.
            'ifort' or 'gfortran': only that compiler will be checked.
        skip_check : bool, optional
            If True, will not check if the compiler is present. Default False.
        Raises
        ------
        RuntimeError
            If the given compiler is not supported, or if no compiler is found.

        Returns
        -------
        None.

        """
        # supported compilers, in order of priority
        supported = ["ifort", "gfortran"]
        found = ""
        if comp == "auto":
            check = supported
        elif comp in supported:
            check = [comp]
        else:
            logger.error("Rparams.getFortranComp: requested compiler is not "
                         "supported.")
            raise RuntimeError("Fortran compiler not supported.")
        if not skip_check:
            for c in check:
                if shutil.which(c, os.X_OK) is not None:
                    found = c
                    break
        else:
            found = check[0]
        if found == "":
            if comp == "auto":
                logger.error("Rparams.getFortranComp: No fortran compiler "
                             "found.")
            else:
                logger.error("Rparams.getFortranComp: Requested fortran "
                             "compiler not found.")
            raise RuntimeError("Fortran compiler not found.")
        if found == "ifort":
            self.FORTRAN_COMP = [
                "ifort -O2 -I/opt/intel/mkl/include",
                "-L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 "
                "-lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl"]
            logger.debug("Using fortran compiler: ifort")
        elif found == "gfortran":
            self.FORTRAN_COMP = ["gfortran -O2", "-llapack -lpthread -lblas"]
            logger.debug("Using fortran compiler: gfortran")
        return

    def getFortranMpiComp(self, comp="auto", skip_check=False):
        """
        Checks whether mpiifort or mpifort are present, and sets FORTRAN_COMP
        mpi accordingly.

        Parameters
        ----------
        comp : str, optional
            'auto' (default): will check mpiifort first, then mpifort.
            'mpiifort' or 'mpifort': only that compiler will be checked.
        skip_check : bool, optional
            If True, will not check if the compiler is present. Default False.

        Raises
        ------
        RuntimeError
            If the given compiler is not supported, or if no compiler is found.

        Returns
        -------
        None.

        """
        # supported compilers, in order of priority
        supported = ["mpiifort", "mpifort"]
        found = ""
        if comp == "auto":
            check = supported
        elif comp in supported:
            check = [comp]
        else:
            logger.error("Rparams.getFortranMpiComp: requested compiler is "
                         "not supported.")
            raise RuntimeError("Fortran MPI compiler not supported.")
        if not skip_check:
            for c in check:
                if shutil.which(c, os.X_OK) is not None:
                    found = c
                    break
        else:
            found = check[0]
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
            self.FORTRAN_COMP_MPI = ["mpifort -Ofast -no-pie -fallow-argument-mismatch", ""]
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
        """Returns a list of 'centered' parameter indices, i.e. all as close
        to 'no displacement' as possible."""
        return ([sp.center for sp in self.searchpars])

    def getPredictConfig(self, best_config=None, curv_cutoff=1e-4):
        """
        Outputs a parameter configuration as list of integers that is as close
        to the result of the parabola fit as possible. For parameters with no
        valid parabola fit result, or curvatures below the cutoff, the value
        from the best known configuration will be used instead. If the best
        known configuration is not passed, these parameters will instead
        assume the value in the center of their displacement range.

        Parameters
        ----------
        best_config : list of int, optional
            The current best configuration, to be used for parameters that
            do not have usable parabola fit results. If not passed, centered
            values will be used instead.
        curv_cutoff : float, optional
            The minimum curvature that a parameter parabola needs to have in
            order to use the parabola minimum. The default is 1e-4.

        Returns
        -------
        list of int
            Parameter values for a new configuration.

        """
        out = []
        for (i, sp) in enumerate(self.searchpars):
            if (sp.parabolaFit["min"] is not None and
                    not np.isnan(sp.parabolaFit["err_co"]) and
                    not np.isnan(sp.parabolaFit["err_unco"])):
                # sp.parabolaFit["curv"] is not None and
                # sp.parabolaFit["curv"] > curv_cutoff):
                out.append(int(round(sp.parabolaFit["min"])))
            else:
                if best_config is not None:
                    out.append(best_config[i])
                else:
                    out.append(int((sp.steps + 1)/2))
        return self.renormalizeDomainParams(out)

    def getRandomConfig(self):
        """
        Outputs a new random parameter configuration. Linked parameters are
        guaranteed to have the same value.

        Returns
        -------
        list of int
            Parameter values for a new configuration.

        """
        out = [-1] * len(self.searchpars)
        for i, sp in enumerate(self.searchpars):
            out[i] = random.randint(1, sp.steps)
        for i, sp in enumerate(self.searchpars):
            if sp.linkedTo is not None:
                out[i] = out[self.searchpars.index(sp.linkedTo)]
            elif sp.restrictTo is not None:
                if type(sp.restrictTo) == int:
                    out[i] = sp.restrictTo
                else:
                    out[i] = out[self.searchpars.index(sp.restrictTo)]
        if -1 in out:
            logger.error("Rparams.getRandomConfig failed: {}".format(out))
            return []
        return self.renormalizeDomainParams(out)

    def getOffspringConfig(self, parents):
        """
        Returns a list of parameter indices generated as a random mix of
        the parameters from two of the configurations passed as 'parents',
        picked at random if there are more than two. If possible, the parents
        will be picked such that they are not identical, but the offspring may
        nevertheless be identical with one of the parents.

        Parameters
        ----------
        parents : list of list of int
            List of surviving configurations.

        Returns
        -------
        list of int
            Parameter values for a new configuration.

        """
        parents = parents[:]
        if len(parents) < 2:
            if len(parents) == 0:
                logger.warning("Rparams.getOffspringConfig: Cannot create "
                               "offspring configuration without parents. "
                               "Returning random configuration")
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
        out = [-1] * len(self.searchpars)
        for i, sp in enumerate(self.searchpars):
            if sp.restrictTo is None and sp.linkedTo is None:
                out[i] = p2[random.randint(0, 1)][i]
        for i, sp in enumerate(self.searchpars):
            if sp.restrictTo is not None:
                if type(sp.restrictTo) == int:
                    out[i] = sp.restrictTo
                else:
                    out[i] = out[self.searchpars.index(sp.restrictTo)+1]
            elif sp.linkedTo is not None:
                out[i] = out[self.searchpars.index(sp.linkedTo)+1]
        if -1 in out:
            logger.error("Rparams.getOffspringConfig failed: {}".format(out))
            return []
        return self.renormalizeDomainParams(out)

    def closePdfReportFigs(self):
        """
        Closes the pdf figures from the Search-progress pdf files, which are
        kept in memory for the Search-report file during a run.

        Returns
        -------
        None.

        """
        global plotting
        if not plotting:
            return

        for searchname in self.lastParScatterFigs:
            for f in self.lastParScatterFigs[searchname]:
                try:
                    plt.close(f)
                except Exception:
                    pass
        return

    def generateSearchPars(self, sl, subdomain=False):
        """
        Initializes a list of Searchpar objects, and assigns delta files to
        atoms if required. Also sets the self.indyPars parameter.

        Parameters
        ----------
        sl : Slab
            The Slab object containing atom information.
        subdomain : bool, optional
            Set True if executing for multiple domains, and this is not the
            highest-level Rprarams object.

        Raises
        ------
        ValueError
            If inconsistent values are found.
        RuntimeError
            If required files are missing.

        Returns
        -------
        None.

        """
        if self.domainParams:
            return(self.generateSearchPars_domains())
        self.searchpars = []
        self.search_maxfiles = 0   # maximum number of delta files for one atom
        self.search_maxconc = 1
        self.indyPars = 0        # number of independent parameters
        self.mncstep = 0     # max. steps (geo. times therm.) for one atom
        # track which atoms are symmetry-linked to the ones already done to
        #   not double-count indyPars
        eqlist = []
        # get list of atoms that appear in the search
        if (2 in self.runHistory or 42 in self.runHistory
                or sl.deltas_initialized):
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
                except ValueError:
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
                for i in range(0, len(occlists[0])):
                    totalocc = 0.
                    for ol in occlists:
                        if len(ol) <= i:
                            logger.error("Inconsistent occupancy lists for {}"
                                         .format(at))
                            raise ValueError("Inconsistent input")
                        else:
                            totalocc += ol[i]
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
                        if np.linalg.norm(at.disp_geo[el][0]) >= 1e-4:
                            found = True
                            break
                if not found:
                    occlists = []
                    for k in at.disp_occ:
                        occlists.append(at.disp_occ[k])
                    for i in range(0, len(occlists[0])):
                        totalocc = 0.
                        for ol in occlists:
                            if len(ol) <= i:
                                break  # error - will pop up again later...
                            else:
                                totalocc += ol[i]
                        if totalocc < 1 - 1e-4:
                            found = True
                            break
                if found and at not in atlist:
                    logger.error("Atom {} has displacements defined, but no "
                                 "delta file was found! Run Delta-Amplitudes."
                                 .format(at.oriN))
                    raise RuntimeError("Delta file not found")
                elif not found and at in atlist:
                    # delta file is there, but no displacements
                    atlist.remove(at)
            sl.deltas_initialized = True
        # sort atlist by displists
        al2 = []
        for dl in sl.displists:
            al2.extend([a for a in atlist if a in dl])
        al2.extend([a for a in atlist if a not in al2])
        atlist = al2
        md = {"geo": 1, "vib": 2, "occ": 3}
        splToRestrict = []
        indep = []
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
                                         ((mode == "geo"
                                          # and np.linalg.norm(d[k][0]) >= 1e-4
                                           )
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
                        elif len(d[k]) > 1 and at not in eqlist:
                            self.indyPars += 1
                            indep.append(sp)
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
                    if at not in eqlist:  # occupation will actually vary
                        self.indyPars += 1
                        indep.append(sp)
                if at in eqlist:
                    spl = [s for s in self.searchpars if at
                           in s.atom.displist and s.mode == "occ"]
                    if spl:
                        sp.linkedTo = spl[0]
            eqlist.extend(at.displist)  # do not consider for future indyPars
        splTargets = set()
        for (sp, (at, el)) in splToRestrict:
            # ind = self.searchpars.index(sp)
            spl = [s for s in self.searchpars if s != sp and s.atom == at and
                   s.mode == sp.mode and (s.el == el or
                                          el in ["", "all"] or
                                          sp.mode == "occ")]
            if spl:
                sp.restrictTo = spl[0]
                splTargets.add(spl[0])
                if spl[0] not in indep and spl[0].steps > 1:
                    self.indyPars += 1
                    indep.append(spl[0])
        for (sp, (at, el)) in [tup for tup in splToRestrict
                               if tup[0].restrictTo is None
                               and tup[0] not in splTargets]:
            logger.warning(
                "Restricting search parameter for atom {}, "
                "element {}, mode {} failed: Could not identify target "
                "search parameter (atom {}, element {})."
                .format(sp.atom.oriN, sp.el, sp.mode, at.oriN, el))
            if sp not in indep and sp.steps > 1:
                self.indyPars += 1
                indep.append(sp)
        for (i, sp) in enumerate(self.searchpars):
            # restrict to lowest number index, resolve conflicts
            if sp.restrictTo not in self.searchpars:
                continue
            sp2 = None
            while sp2 != sp.restrictTo:
                sp2 = sp.restrictTo
                ind = self.searchpars.index(sp2)
                if (sp2.linkedTo in self.searchpars
                        and self.searchpars.index(sp2.linkedTo) < ind):
                    sp.restrictTo = sp2.linkedTo
                elif (sp2.restrictTo in self.searchpars
                      and self.searchpars.index(sp2.restrictTo) < ind):
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

    def generateSearchPars_domains(self):
        """
        Runs generateSearchPars for every domain, then collates results.

        Returns
        -------
        None.

        """
        self.searchpars = []
        self.indyPars = len(self.domainParams) - 1
        home = os.getcwd()
        for dp in self.domainParams:
            try:
                os.chdir(dp.workdir)
                dp.rp.generateSearchPars(dp.sl, subdomain=True)
            except Exception:
                logger.error("Error while creating delta input for domain {}"
                             .format(dp.name))
                raise
            finally:
                os.chdir(home)
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
