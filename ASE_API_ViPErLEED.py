# -*- coding: utf-8 -*-
"""A collection of functions that run ViPErLEED from a provided ASE object.

@Author: Alexander M. Imre
based on tleedm_from_ase.py by Florian Kraushofer

Requires viperleed to be in the system path.
"""

import sys
import os
import numpy as np
import shutil
import ase.db
import lzma, base64
from io import StringIO

import viperleed

# for refcalc_from_ase_structure
from viperleed.tleedm import run_tleedm
from viperleed.tleedmlib import Slab, Rparams
from viperleed.tleedmlib.files.poscar import writePOSCAR

# for run_rfactor_from_csv
from viperleed.tleedmlib.files.beams import readOUTBEAMS
from viperleed.tleedmlib.files.iorfactor import beamlist_to_array
from viperleed.tleedmlib.wrapped.rfactor import r_factor_new as rf
from viperleed.tleedmlib.wrapped.error_codes import error_codes


def refcalc_for_ase_structure(
    exec_path,
    ase_object,
    inputs_path=None,
    cut_symmetric_cell_height_fraction=0.4,
    use_which_part_of_cell="keep_above",
    uc_transformation_matrix=None,
    uc_scaling=None,
):
    """Performs a ViPErLEED reference calulation starting from a provided ASE atoms object.
    
    This function is a core part of the API for ViPErLEED. It alows to run a reference calulation from a ASE type object.
    The ASE object which is expected to conatin a symmetric simulation cell will be cut at a specified height fraction, transformed into a ViPErLEED Slab and written to a POSCAR file. The ViPErLEED reference calculation will then be executed on the given structure using the input files found either in inputs_path or exec_path. This ultimatetly yiels reference electron scattering amplitudes and intensities. These are collected in CSV files which are returned as strings for further processing.
    
    Parameters
    ----------
    exec_path: string
        Path where to execute the reference calculation. The directory needs to exist. Input files will be copied from inputs_path if provided - otherwise it is assumed that all input files are already in exec_path. Temporary files will be stored in a subdirectory called work.
    ase_object : ASE atoms object
        ASE type object that contains the atoms for the slab to be used with ViPErLEED. This is expected to contain a symmetric simulation cell which will be cut as described below.
    inputs_path: string, optional
        Path to input files (PARAMETERS, VIBROCC, etc.) for ViPErLEED calculation. If provided, files will be copied from inputs_path to exec_path (overwriting existing files there!).
    cut_symmetric_cell_height_fraction : float, default = 0.4
        Where to cut the ase object cell (ViPErLEED can not use symmetric cells). Cutting will be performed along the vector c AFTER all transformations (see uc_transformation_matrix) are applied.
    use_which_part_of_cell: string, default="keep_above"
        Wether to keep the part above or below cut_symmetric_cell_height_fraction (use "keep_above" or "keep_below").
    uc_transformation_matrix : np.ndarray, optional
        3x3 np array, can contain an arbitray orthogonal transforamtion that will be applied to the cell. Default is None which implies no transformation. The transformation will be applied BEFORE scaling, i.e. using the original unit cell.
        The transformation is described by an orthogonal transfomation matrix (O) which is applied to both the unit cell and all atomic positions. This transformation is essentially equivalent to a change of basis.
        The transformation of the unit cell (U) is of the form O^(-1)UO. Atom positions (v) are transformed as vOS.
        This parameter can be used to switch indices, rotate the slab, etc..
    uc_scaling : float or list of float, default = None
        Scaling to be applied to the cell. Default is None which implies no scaling. Scaling will be applied AFTER the transformation, i.e. along the new unit cell vectors.
        This can be used to apply an isotropic scaling in all directions or to strech/compress along unit cell vectors in order to change lattice constants in some direction.
        If the scaling is given as a float or a list containing one item, an isotropic scaling will be applied to the unit cell and atom positions. If a list with 3 entries is provided, the scaling will be applied along the unit cell vector in the given order.
        
    Returns
    -------
    theobeams_file_str : string
        String containg the contents of the file "THEOBEAMS.csv", i.e. the scattering intensities as calculated in the reference calculation.
    amp_real_file_str : string
        String containg the contents of the file "Complex_amplitudes_real.csv", i.e. the real part of the scattering amplitudes as calculated in the reference calculation.
    amp_imag_file_str : string
        String containg the contents of the file "Complex_amplitudes_imag.csv", i.e. the imaginary part of the scattering amplitudes as calculated in the reference calculation.
        
    Notes
    -------
    The parameter uc_transformation matrix can be used to apply some useful transformations to the Slab object. Examples include swapping vectors b and c, rotations around an axis and flipping the cell along c (in case the lower half of the symmetric cell is kept).
    For these standard applications, have a look at the matrices we provide as part of this file:
        switch_b_c_mat : switches vectors b & c
        switch_a_b_mat : switches vectors a & b
        flip_c_mat : flips cell along c
        rot_mat_c(theta) : function that generates a rotation matrix around c (by angle theta)
        rot_mat_a(theta) : function that generates a rotaiton matrix around a (by angle theta)
    """

    if not os.path.isdir(exec_path):
        # Invalid path given
        raise RuntimeError("Provided exec_path is invalid")

    input_files = [
        "PARAMETERS",
        "VIBROCC",
        "IVBEAMS",
        "PHASESHIFTS",
        "DISPLACEMENTS",
        "EXPBEAMS",
        "EXPBEAMS.csv",
    ]

    # Copy all files in the input Path
    if inputs_path is not None:
        for file in input_files:
            try:
                shutil.copy2(
                    os.path.join(inputs_path, file), os.path.join(exec_path, file)
                )
            except FileNotFoundError:
                pass

    # check for PARAMETERS file - without that we can't proceede
    parameters_name = "PARAMETERS"
    if not os.path.isfile(os.path.join(exec_path, parameters_name)):
        # No PARAMETERS file – Error
        raise RuntimeError("No PARAMETERS file found – this is required")
    # If present, we are good to go.
    # Files IVBEAMS, PHASESHIFTS and VIBROCC are not required and can be autogenerated if PARAMETERS are set correctly

    # Get temporary parameters object: dummy but needed for getSurfaceAtoms
    rp = Rparams()

    # Transfer ASE object into slab object for ViPErLEED
    slab = Slab(ase_atoms=ase_object)

    # Transformation of slab object: Rotation or isotropic streching/shrinking
    if uc_transformation_matrix is not None:
        slab.apply_matrix_transformation(uc_transformation_matrix)
    elif uc_isotropic_scaling is not None:
        slab.apply_isotropic_scaling(uc_isotropic_scaling)

    # Cut the symmetric slab in half so we can work with just the surface
    # Remove everything below cut_fraction
    if use_which_part_of_cell.lower() == "keep_above":
        slab.atlist = [
            at for at in slab.atlist if at.pos[2] > cut_symmetric_cell_height_fraction
        ]
    elif use_which_part_of_cell.lower() == "keep_below":
        slab.atlist = [
            at for at in slab.atlist if at.pos[2] < cut_symmetric_cell_height_fraction
        ]
    else:
        raise RuntimeError(
            "use_which_part_of_cell must be either 'keep_above' or 'keep_below'"
        )

    slab.updateAtomNumbers()
    slab.updateElementCount()

    # Figure out surface sites
    site_def = {}
    surface_atoms = slab.getSurfaceAtoms(rp)
    # surface species of each element
    for el in slab.elements:
        atn = [at.oriN for at in surface_atoms if at.el == el]
        if atn:
            site_def[el] = {"surf": atn}

    poscar_name = "POSCAR"
    if os.path.isfile(os.path.join(exec_path, poscar_name)):
        raise RuntimeError(
            "A 'POSCAR' file is already present in exec_path - calculation aborted. Check if a ViPErLEED calculation was already run in this directroy."
        )
    # If no POSCAR - good :)
    writePOSCAR(slab, poscar_name)

    # Take care of input files and work directory
    work_path = os.path.abspath(os.path.join(exec_path, "work"))
    os.makedirs(work_path, exist_ok=True)
    # copy input files to work directory
    # NOTE: PARAMETERS should NOT contain SITE_DEF flags.
    # VIBROCC should contain *_surf sites for all elements.
    for file in input_files:
        try:
            shutil.copy2(os.path.join(exec_path, file), os.path.join(work_path, file))
        except FileNotFoundError:
            pass

    home = os.path.abspath(".")
    os.chdir(work_path)

    # We are ready to run ViPErLEED! Have fun!
    try:
        run_tleedm(
            slab=slab,
            preset_params={"SITE_DEF": site_def},
            source=os.path.dirname(viperleed.__file__),
        )
    except:
        # If ViPErLEED fails, move back to home directory
        os.chdir(home)
        raise RuntimeError("ViPErLEED referce calculation failed.")

    # ViPErLEED should have suceeded if you arrive here!

    # read out the THEOBEAMS.csv file and complex ampliudes:
    theobeams_name = "THEOBEAMS.csv"
    amp_real_name = "Complex_amplitudes_real.csv"
    amp_imag_name = "Complex_amplitudes_imag.csv"

    with open(theobeams_name) as f:
        theobeams_file_str = f.read()
    with open(amp_real_name) as f:
        amp_real_file_str = f.read()
    with open(amp_imag_name) as f:
        amp_imag_file_str = f.read()

    # Move back home
    os.chdir(home)

    return (theobeams_file_str, amp_real_file_str, amp_imag_file_str)


def run_rfactor_from_csv(
    beams_files,
    V0i,
    beams_file_is_content = (False, False),
    V0r_shift_range=[-3, 3],
    intpol_deg=5,
    intpol_step=0.5,
):
    """A function that compute the Pendry R-factor between two CSV files as produced and used by ViPErLEED.
    
    Reads in two CSV files containing LEED-I(V) spectra and computes the mutual R-factor. The format of the files must be the same as accepted/generated by ViPErLEED.
    
    Parameters
    ----------
    beams_files : touple of two string
        The two CSV files to be read in for R-factor calculation. These can contain either theoretical or experimental beams.
        If beams_file_is_content is set as True, the string is expected to contain the contents of the CSV file, otherwise it should be a path to the file.
    V0i : float
        Imaginary part of the inner potential (in eV).
    beams_file_is_content : touple of two bool, default = (False, False)
        If either value is set to true, it is assumed that the corresponding element of beams_file contains a string witht the CSV contents rather than a path to the file. For reading, a StringIO will be used internally.
    V0r_shift_range : list of two float, default = [-3, 3]
        During R-factor calculation, the real part of the inner potential is varied within these bounds (given in eV). It is recommended to choose multiples of intpol_step as bounds.
    intpol_deg : int, default = 5
        Either 3 or 5; degree of the interpolating spline polynomial
    intpol_step : float, default = 0.5
        Step size of the energy grid the beam data will be interpolated on (in eV)
        
    Returns
    -------
    best_R : float
        The R-factor between the two files. Returns the best R factor found after variation in the range or V0r_shift.
    best_V0r : float
        The V0r value for which the best R-factor was found.
    """
    # It makes sense for the V0r_shift_range to be a multiple of intpol_step, but this neither strictly required nor enforced

    # Use Pendry R-factor - TODO: dicuss if this should be a user parameter in the API (if so, it code below will need a few ajustments)
    which_r = 1

    err_msg = ""

    # Check if file is supposed to be read from disk or from memory and pass it to readOUTBEMS accordingly.
    if (beams_file_is_content[0]):
        io_file = SringIO(beams_files[0])
        beams1 = readOUTBEAMS(filename = None, sep = ",", file_StringIO = io_file)
    else:
        beams1 = readOUTBEAMS(filename = beams_files[0], sep = ",")
        
    if (beams_file_is_content[1]):
        io_file = SringIO(beams_files[1])
        beams2 = readOUTBEAMS(filename = None, sep = ",", file_StringIO = io_file)
    else:
        beams2 = readOUTBEAMS(filename = beams_files[1], sep = ",")

    # Check for common beams
    beam1_has_correspondece = [False] * len(beams1)
    beam2_correspondence = [None] * len(beams2)

    for i, beam1 in enumerate(beams1):
        for j, beam2 in enumerate(beams2):
            if _comp_beam_index(beam1.hkfrac, beam2.hkfrac):
                beam1_has_correspondece[i] = True
                beam2_correspondence[j] = i

    for beam_id in beam2_correspondence:
        if (beam2_correspondence.count(beam_id) != 1) and (beam_id is not None):
            raise RuntimeError("Multiple beams sharing the same index encountered!")

    corr_beams1 = []
    corr_beams2 = []
    # any beams without correspondence?
    for i, beam1 in enumerate(beams1):
        if beam1_has_correspondece[i]:
            corr_beams1.append(beam1)
            corr_beams2.append(beams2[beam2_correspondence.index(i)])
        else:
            err_msg += f"Beam {beam1.label} from file 1 has no correspondence\n"
    for j, beam2 in enumerate(beams2):
        if beam2_correspondence[j] is None:
            err_msg += f"Beam {beam2.label} from file 2 has no correspondence\n"

    if err_msg:
        print(err_msg)

    # Done with checking beam indices
    n_beams = len(corr_beams1)

    # Now get beams into arrays
    beams1_en, beam1_id_start, beam1_n_E_beams, beam1_arr = beamlist_to_array(
        corr_beams1
    )
    beams2_en, beam2_id_start, beam2_n_E_beams, beam2_arr = beamlist_to_array(
        corr_beams2
    )

    # extend range by V0r_shift_rang to accomodate maximum possible shift
    minen = min(np.min(beams1_en), np.min(beams2_en)) + V0r_shift_range[0]
    maxen = max(np.max(beams1_en), np.min(beams2_en)) + V0r_shift_range[1]

    averaging_scheme = np.int32(np.arange(n_beams) + 1)

    grid = np.arange(minen, maxen + intpol_step, intpol_step)

    skip_stages = np.int32([0, 0, 0, 0, 0]) # usused - we don't skip anything
    averaging_scheme = np.int32(np.arange(n_beams) + 1) # don't average

    # interpolate and prepare beam arrays
    (
        beams1_e_start_beams_out,
        beams1_n_e_beams_out,
        beams1_intpol_intensity,
        beams1_yfunc,
        ierr,
    ) = rf.prepare_beams(
        beams1_en,
        beam1_arr,
        beam1_id_start + 1,
        beam1_n_E_beams,
        skip_stages,
        n_beams,
        averaging_scheme,
        which_r,
        intpol_deg,
        grid,
        V0i,
    )

    _check_ierr(ierr)

    (
        beams2_e_start_beams_out,
        beams2_n_e_beams_out,
        beams2_intpol_intensity,
        beams2_yfunc,
        ierr,
    ) = rf.prepare_beams(
        beams2_en,
        beam2_arr,
        beam2_id_start + 1,
        beam2_n_E_beams,
        skip_stages,
        n_beams,
        averaging_scheme,
        which_r,
        intpol_deg,
        grid,
        V0i,
    )

    _check_ierr(ierr)

    # Ready to calculate R-factor
    V0r_shift_range_int = np.int32(
        [int(V0r_shift_range[0] / intpol_step), int(V0r_shift_range[1] / intpol_step)]
    )
    V0r_center = int((V0r_shift_range_int[0] + V0r_shift_range_int[1]) / 2)
    start_guess = np.int32([V0r_shift_range_int[0], V0r_center, V0r_shift_range_int[1]])

    # Same settings as used in R-factor for ViPErLEED
    fast_search = False  # TODO change to True as default if we are sure it is better
    tol_r = (1 - 5e-2,)
    tol_r_2 = (1 - 5e-2,)
    max_fit_range = 6

    (
        best_V0r_step,
        best_V0r,
        best_R,
        n_V0r_evaluated,
        R_beams,
        numerators,
        denominators,
        n_overlapp_beams,
        ierr,
    ) = rf.r_beamset_v0r_opt_on_grid(
        which_r,
        V0r_shift_range_int,
        start_guess,
        fast_search,
        intpol_step,
        beams1_yfunc,
        beams2_yfunc,
        beams1_e_start_beams_out,
        beams2_e_start_beams_out,
        beams1_n_e_beams_out,
        beams2_n_e_beams_out,
        tol_r=tol_r,
        tol_r_2=tol_r_2,
        max_fit_range=max_fit_range,
    )
    # If encountered an error
    _check_ierr(ierr)

    return best_R, best_v0r


def _comp_beam_index(id1, id2, eps=1e-5):
    # Compares diffraction beam indices and checks if they refer to the same beam
    check = np.abs(float(id1[0] - id2[0])) + np.abs(float(id1[1] - id2[1])) < eps
    return check


def _check_ierr(ierr):
    if ierr != 0:
        raise RuntimeError(f"ViPErLEED Fortran error code {ierr}: {error_codes[ierr]}")


def rot_mat_a(theta):
    """Generates a rotation matrix around the a axis.
    
    Parameters
    ----------
    theta : flaot
        Angle of rotation in radians (use e.g. np.radians(deg) to convert from degrees).
    """

    rot_mat = [
        [1, 0, 0],
        [0, np.cos(theta), -np.sin(theta)],
        [0, np.sin(theta), np.cos(theta)],
    ]
    rot_mat = np.array(rot_mat)
    return rot_mat


def rot_mat_c(theta):
    """Generates a rotation matrix around the c axis.
    
    Parameters
    ----------
    theta : flaot
        Angle of rotation in radians (use e.g. np.radians(deg) to convert from degrees).
    """

    rot_mat = [
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1],
    ]
    rot_mat = np.array(rot_mat)
    return rot_mat


# some other usefull matrices
switch_b_c_mat = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]])  # switches vectors b & c
switch_a_b_mat = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])  # switches vectors a & b

flip_c_mat = np.array(
    [[1, 0, 0], [0, 1, 0], [0, 0, -1]]
)  # flips the cell along c (usefull in case the lower half of the symmetric cell is kept)
