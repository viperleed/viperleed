# -*- coding: utf-8 -*-
"""A collection of functions that run ViPErLEED from a provided ASE object.

@Author: Alexander M. Imre, Michele Riva
based on tleedm_from_ase.py by Florian Kraushofer

Requires viperleed to be in the system path (or on PYTHONPATH).
"""

import os
import warnings
from pathlib import Path
import numpy as np
import shutil
from io import StringIO

import viperleed

# for run_from_ase
from viperleed.tleedm import run_tleedm
from viperleed.tleedmlib import Slab, Rparams
from viperleed.tleedmlib.files.poscar import writePOSCAR
from viperleed.tleedmlib.files.parameters import readPARAMETERS, interpretPARAMETERS

# for rfactor_from_csv
from viperleed.tleedmlib.files.beams import readOUTBEAMS
from viperleed.tleedmlib.files.iorfactor import beamlist_to_array
from viperleed.tleedmlib.wrapped.rfactor import r_factor_new as rf
from viperleed.tleedmlib.wrapped.error_codes import error_codes

# for plot_iv_from_csv
from viperleed.tleedmlib.files.ivplot import plot_iv
from typing import Sequence


def run_from_ase(
    exec_path,
    ase_object,
    inputs_path=None,
    cut_cell_c_fraction=0.4,
    uc_transformation_matrix=None,
    uc_scaling=None,
):
    """Perform a ViPErLEED calulation starting from an ase.Atoms object.
    
    This function is a core part of the API for ViPErLEED. It allows to run
    a calulation from an ASE-type object.  The ASE object which is expected
    to contain a simulation cell will be cut at a specified height fraction,
    transformed into a ViPErLEED Slab and written to a POSCAR file. The
    ViPErLEED calculation will then be executed on the given structure using
    the input files found either in inputs_path or exec_path. The type of
    calculation performed depends on the value of the RUN parameter in the
    PARAMETERS file. If RUN = 1, this ultimately yields reference electron
    scattering amplitudes and intensities. These are collected in CSV files
    which are returned as strings for further processing.
    
    Parameters
    ----------
    exec_path: string
        Path where to execute the reference calculation. The directory needs
        to exist. Input files will be copied from inputs_path, if provided.
        Otherwise it is assumed that all input files are already in exec_path.
        Temporary files will be stored in a subdirectory called work.
    ase_object : ase.Atoms
        ASE-type object that contains the atoms for the slab to be used with
        ViPErLEED. This is expected to contain a simulation cell which will
        be cut as described below.
    inputs_path: string, optional
        Path to directory containing input files (PARAMETERS, VIBROCC, etc.)
        for ViPErLEED calculation. If provided, files will be copied from
        inputs_path to exec_path (overwriting existing files there!).
    cut_cell_c_fraction : float, default=0.4
        Where to cut the ase object cell (ViPErLEED expects the cell to have
        the surface at the "top" -- i.e., at high z coordinates, and bulk-like
        layers at the bottom). Cutting will be performed along the vector c
        AFTER all transformations (see uc_transformation_matrix) are applied.
        Only the part of the cell above (i.e., >=) this value will be
        kept and everything else will be discarded.
    uc_transformation_matrix : np.ndarray, optional
        3x3 array, can contain an arbitray orthogonal transforamtion (i.e.,
        a combination of mirroring and rotations) that will be applied to the
        cell. Default is None, i.e., no transformation. The transformation is
        applied BEFORE scaling, i.e., using the original unit cell.
        The transformation is described by an orthogonal transfomation matrix
        (O) which is applied to both the unit cell and all atomic positions.
        This transformation is essentially equivalent to a change of basis.
        The transformation of the unit cell (U) is of the form OUO^T. Atom
        positions (v, both fractional and cartesian) are transformed as Ov.
        This parameter can be used to switch indices, rotate the slab, etc..
        Make sure that vectors a & b do not have any components in z
        direction after transformation as this is not allowed and will raise
        an error.
    uc_scaling : float or list of floats, optional
        Scaling to be applied to the cell. Default is None, i.e., no scaling.
        Scaling will be applied AFTER the transformation, i.e. along the new
        unit cell vectors. This can be used to apply an isotropic scaling in
        all directions or to strech/compress along unit cell vectors in order
        to change lattice constants in some direction. If the scaling is given
        as a number or a 1-element sequence, an isotropic scaling will be
        applied to the unit cell and atom positions. If a sequence with 3
        entries is provided, the scaling will be applied along the unit cell
        vectors in the given order. Given uc_scaling == (s1, s2, s3), the
        new unit cell vectors will be a' = a * s1, b' = b * s2, c' = c * s3.

    Returns
    -------
    theobeams_file_str : string
        String containg the contents of the file "THEOBEAMS.csv", i.e. the
        scattering intensities as calculated in the reference calculation.
        An empty string is returned if no reference calculation was run.
    amp_real_file_str : string
        String containg the contents of the file "Complex_amplitudes_real.csv",
        i.e., the real part of the scattering amplitudes as calculated in the
        reference calculation. An empty string is returned if no reference
        calculation was run.
    amp_imag_file_str : string
        String containg the contents of the file "Complex_amplitudes_imag.csv",
        i.e., the imaginary part of the scattering amplitudes as calculated in
        the reference calculation. An empty string is returned if no reference
        calculation was run.
    v0i : float
        Imaginary part of the inner potential as read from "PARAMETERS". This is
        returned here since it is an input parameter for the R-factor calculation.
        The value is read from PARAMETERS and left unchanged - it is not a result
        of the calculation and is returned here for convenience only.
    
    Raises
    ------
    ValueError
        If exec_path is non-existent (or not a directory)
    RuntimeError
        If no PARAMETERS file was found in exec_path (after potentially
        copying over from inputs_path)
    RuntimeError
        If either of the first two unit cell vectors have a component
        perpendicular to the surface (i.e., along the z coordinate)
        after transformation of the unit cell
    
    Notes
    -------
    The parameter uc_transformation_matrix can be used to apply some useful
    transformations to the Slab object. Examples include swapping vectors b
    and c, rotations around an axis and flipping the cell along c (in case
    the lower part has to be kept, rather than the top one). For these standard
    applications, have a look at the matrices we provide as part of this module:
        switch_b_c_mat : Switches vectors b & c
        switch_a_b_mat : Switches vectors a & b
        flip_c_mat : Flips (i.e., mirrors) cell along c
        rot_mat_c(theta) : Returns a rotation matrix around c (by theta rad.)
        rot_mat_a(theta) : Returns a rotation matrix around a (by theta rad.)
    For applying multiple operations, the order does matter. You can combine
    operations via matrix multiplication (@ symbol or np.dot); the left-most
    matrix is applied last.
    """

    exec_path = Path(exec_path)
    if not exec_path.is_dir():
        # Invalid path given
        raise ValueError("Invalid exec_path: not existent, or not a directory")

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
        inputs_path = Path(inputs_path)
        for file in input_files:
            try:
                shutil.copy2(inputs_path / file, exec_path / file)
            except FileNotFoundError:
                pass

    # check for PARAMETERS file - without that we can't proceede
    parameters_name = "PARAMETERS"
    if not (exec_path / parameters_name).exists():
        # No PARAMETERS file – Error
        raise RuntimeError("No PARAMETERS file found – this is required")
    # If present, we are good to go.
    # Files PHASESHIFTS and VIBROCC are not required and can be
    # autogenerated if PARAMETERS are set correctly. IVBEAMS can
    # be generated if EXPBEAMS exists.

    # Get temporary parameters object
    rp = readPARAMETERS(exec_path / parameters_name)
    interpretPARAMETERS(rp)
    v0i = rp.V0_IMAG

    # Transfer ASE object into slab object for ViPErLEED
    slab = Slab(ase_atoms=ase_object)

    # Transformation of slab object: Rotation or isotropic streching/shrinking
    if uc_transformation_matrix is not None:
        slab.apply_matrix_transformation(uc_transformation_matrix)
    if uc_scaling is not None:
        slab.apply_scaling(uc_scaling)

    # check if the now transformed slab has any z components in vectors a & b
    # raise an error because this would mess up parts of the TensErleed
    # calculations (as we pass on only the first two components).

    a, b, _ = slab.ucell.T
    if abs(a[2]) > 1e-5 or abs(b[2]) > 1e-5:
        raise RuntimeError(
            "z component found in unit cell vector a or b. This is not "
            "allowed in the TensErLEED calculation. Check eventual "
            "transformations applied to the unit cell and make sure a "
            "and b are parallel to the surface."
        )

    # Cut the slab as given in cut_cell_c_fraction - very important
    slab.atlist = [at for at in slab.atlist if at.pos[2] >= cut_cell_c_fraction]
    # Update Slab
    slab.updateAtomNumbers()
    slab.updateElementCount()

    # Figure out surface sites
    site_def = {}
    surface_atoms = slab.getSurfaceAtoms(rp)
    # surface species of each element
    for element in slab.elements:
        atn = [at.oriN for at in surface_atoms if at.el == element]
        if atn:
            site_def[element] = {"surf": atn}

    poscar_name = "POSCAR"
    if (exec_path / poscar_name).exists():
        warnings.warn(
            "A 'POSCAR' file is already present in exec_path - calculation "
            "aborted. Check if a ViPErLEED calculation was already run in "
            "this directory.",
            category=RuntimeWarning,
        )

    writePOSCAR(slab, poscar_name)

    # Take care of input files and work directory
    work_path = exec_path / "work"
    work_path.mkdir(exist_ok=True)
    # copy input files to work directory
    # NOTE: PARAMETERS should NOT contain SITE_DEF flags.
    # VIBROCC should contain *_surf sites for all elements.
    for file in input_files:
        try:
            shutil.copy2(exec_path / file, work_path / file)
        except FileNotFoundError:
            pass

    home = Path(".").resolve()
    os.chdir(work_path)

    # We are ready to run ViPErLEED! Have fun!
    try:
        run_tleedm(
            slab=slab,
            preset_params={"SITE_DEF": site_def},
            source=os.path.dirname(viperleed.__file__),
        )
    except Exception as err:
        # If ViPErLEED fails, move back to home directory
        os.chdir(home)
        raise RuntimeError("ViPErLEED calculation failed.") from err

    # ViPErLEED should have suceeded if you arrive here. However, we may not
    # have run a refcalc (id == 1). In that case, return empty strings.
    if 1 not in rp.RUN:
        return "", "", "", v0i

    # read out the THEOBEAMS.csv file and complex amplitudes:
    theobeams_name = "THEOBEAMS.csv"
    amp_real_name = "Complex_amplitudes_real.csv"
    amp_imag_name = "Complex_amplitudes_imag.csv"

    with open(theobeams_name, "r", encoding="utf-8") as fproxy:
        theobeams_file_str = fproxy.read()
    with open(amp_real_name, "r", encoding="utf-8") as fproxy:
        amp_real_file_str = fproxy.read()
    with open(amp_imag_name, "r", encoding="utf-8") as fproxy:
        amp_imag_file_str = fproxy.read()

    # Move back home
    os.chdir(home)

    return theobeams_file_str, amp_real_file_str, amp_imag_file_str, v0i


def rfactor_from_csv(
    beams_files,
    v0i,
    beams_file_is_content=(False, False),
    v0r_shift_range=(-3, 3),
    intpol_deg=5,
    intpol_step=0.5,
):  ## TODO: add kwarg for mapping for averaging
    """Compute the Pendry R-factor between two CSV files (in ViPErLEED format).
    
    Read in two CSV files containing LEED-I(V) spectra and compute the mutual
    R-factor. The format of the files must be the same as accepted/generated by
    ViPErLEED.
    
    Parameters
    ----------
    beams_files : tuple of two string
        The two CSV files to be read in for R-factor calculation.
        These can contain either theoretical or experimental beams.
        If beams_file_is_content is set as True, the string is expected
        to contain the contents of the CSV file, otherwise it should be
        a path to the file.
    v0i : float
        Imaginary part of the inner potential (in eV).
    beams_file_is_content : tuple of two bool, default=(False, False)
        If either value is set to True, it is assumed that the
        corresponding element of beams_file contains a string with
        the CSV contents rather than a path to the file. For reading,
        a StringIO will be used internally.
    v0r_shift_range : Sequence of two float, default=(-3, 3)
        During R-factor calculation, the two beam sets can be shifted
        in energy relative to one another. v0r_shift_range determines
        the minimum and maximum shifts (given in eV). The closest
        integer multiple of intpol_step will be used for the bounds.
    intpol_deg : int, default=5
        Either 3 or 5; degree of the interpolating spline polynomial
    intpol_step : float, default=0.5
        Step size of the energy grid the beam data will be interpolated
        on (in eV).
        
    Returns
    -------
    best_R : float
        The R-factor between the two files. Returns the best R factor
        found after variation in v0r_shift_range.
    best_v0r : float
        The v0r value for which the best R-factor was found.
    """
    # Use Pendry R-factor - TODO: dicuss if this should be a user parameter
    #         in the API (if so, the code below will need a few ajustments)
    which_r = 1

    err_msg = ""

    # Check if file is supposed to be read from disk or from memory and pass
    # it to readOUTBEMS accordingly.
    beams_tmp = [None, None]
    for i, (is_content, file) in enumerate(zip(beams_file_is_content, beams_files)):
        if is_content:
            file = StringIO(file)
        beams_tmp[i] = readOUTBEAMS(filename=file, sep=",")

    beams1, beams2 = beams_tmp

    # MRiva: The following checks are potentially time consuming since they
    # require looping through beams. Perhaps we can add an extra keyword
    # argument that skips this part altogether. Then, of course, we have
    # to assume that the beam files contain exactly the same beams in
    # exactly the same order. Not sure the order is consistent between
    # consecutive runs, though (it will be if the order of the output is
    # the same as in IVBEAMS).
    # AMI: Could do, but compared to the refcalc required for the
    # theoretical specetra, it's still very fast.

    # Check if multiple beams have the same index. This should never be
    # the case and would throw off the calculation.

    for i, beams in enumerate((beams1, beams2)):
        unique_lables = set()
        for beam in beams:
            unique_lables.add(beam.hkfrac)
        if len(unique_lables) != len(beams):
            raise RuntimeError(
                f"Multiple beams sharing the same index encountered"
                f" in input beamset {i}!"
            )

    # Check for common beams
    beam_corr_1_2 = {}

    for i, beam1 in enumerate(beams1):
        for j, beam2 in enumerate(beams2):
            if beam1.hkfrac == beam2.hkfrac:
                beam_corr_1_2[beam1] = beam2

    for beam in beams1:
        if beam not in list(beam_corr_1_2.keys()):
            err_msg += f"Beam {beam.label} from file 1 has no correspondence.\n"
    for beam in beams2:
        if beam not in list(beam_corr_1_2.values()):
            err_msg += f"Beam {beam.label} from file 2 has no correspondence.\n"

    corr_beams1 = []
    corr_beams2 = []

    # Warn but don't raise
    if err_msg:
        warnings.warn(err_msg, category=RuntimeWarning)

    # Done with checking beam indices
    n_beams = len(beam_corr_1_2)

    for key, value in beam_corr_1_2.items():
        corr_beams1.append(key)
        corr_beams2.append(value)

    # Now get beams into arrays
    beams1_en, beams1_id_start, beams1_n_E_beams, beams1_arr = beamlist_to_array(
        corr_beams1
    )
    beams2_en, beams2_id_start, beams2_n_E_beams, beams2_arr = beamlist_to_array(
        corr_beams2
    )

    minen1, minen2 = beams1_en.min(), beams2_en.min()
    minen = max(minen1, minen2)
    if abs(minen1 - minen2) < abs(v0r_shift_range[0]):
        minen -= v0r_shift_range[0]

    maxen1, maxen2 = beams1_en.max(), beams2_en.max()
    maxen = min(maxen1, maxen2)
    if abs(maxen1 - maxen2) < abs(v0r_shift_range[1]):
        maxen += v0r_shift_range[1]

    grid = np.arange(minen, maxen + intpol_step, intpol_step)

    skip_stages = np.int32([0, 0, 0, 0, 0])  # usused - we don't skip anything
    averaging_scheme = np.int32(np.arange(n_beams) + 1)  # don't average

    # interpolate and prepare beam arrays
    (
        beams1_e_start_beams_out,
        beams1_n_e_beams_out,
        beams1_intpol_intensity,  # may be needed for R2
        beams1_yfunc,
        ierr,
    ) = rf.prepare_beams(
        beams1_en,
        beams1_arr,
        beams1_id_start + 1,
        beams1_n_E_beams,
        skip_stages,
        n_beams,
        averaging_scheme,
        which_r,
        intpol_deg,
        grid,
        v0i,
    )
    _check_ierr(ierr)

    (
        beams2_e_start_beams_out,
        beams2_n_e_beams_out,
        beams2_intpol_intensity,  # may be needed for R2
        beams2_yfunc,
        ierr,
    ) = rf.prepare_beams(
        beams2_en,
        beams2_arr,
        beams2_id_start + 1,
        beams2_n_E_beams,
        skip_stages,
        n_beams,
        averaging_scheme,
        which_r,
        intpol_deg,
        grid,
        v0i,
    )
    _check_ierr(ierr)

    # Ready to calculate R-factor
    v0r_shift_range_int = np.int32([round(v / intpol_step) for v in v0r_shift_range])
    v0r_center = sum(v0r_shift_range_int) // 2
    start_guess = np.int32([v0r_shift_range_int[0], v0r_center, v0r_shift_range_int[1]])

    # Same settings as used in R-factor for ViPErLEED
    fast_search = False  # TODO change to True as default if we are sure it is better
    tol_r = (1 - 5e-2,)
    tol_r_2 = (1 - 5e-2,)
    max_fit_range = 6

    (_, best_v0r, best_r, *_, ierr) = rf.r_beamset_v0r_opt_on_grid(
        which_r,
        v0r_shift_range_int,
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

    return best_r, best_v0r


def _check_ierr(ierr):
    if ierr:
        raise RuntimeError(f"ViPErLEED Fortran error code {ierr}: {error_codes[ierr]}")


def plot_iv_from_csv(
    beam_file, output_file, beam_file_is_content=False, which_beams=[]
):
    """Plots LEED-I(V) spectra directly from a CSV file.

    This function can be used to plot LEED-I(V) spectra directly from
    a CSV file without running any other calculations.
    It is assumed that the CSV file uses the format used by ViPErLEED.
    Returns nothing but writes the pdf file.


    Parameters
    ----------
    beam_file : string
        If beam_file_is_content is set as True, the string is expected
        to contain the contents of the CSV file, otherwise it should be
        a path to the file.
    beam_file_is_content : bool, default = False
        If set to True, it is assumed that the beam_file
        contains a string with the CSV contents rather than a path to
        the file. For reading, a StringIO will be used internally.
    output_pdf output_pdf : string
        Path to the PDF file where the plot will be saved. (File is
        created.)
    which_beams : range or Sequence of int or strings
        Indices specifying which beams to plot. Order is the same as in the csv.

    Raises
    ------
    RuntimeError
        If which_beams is not Sequence containing just integers or
        strings.
    """

    if beam_file_is_content:
        beam_file = StringIO(beam_file)
    beam_list = readOUTBEAMS(filename=beam_file, sep=",")

    if which_beams == [] or which_beams == "all" or which_beams is None:
        pass
    elif isinstance(which_beams, range) or all(
        [type(label) is int for label in which_beams]
    ):
        beam_list = [beam_list[i] for i in which_beams]
    elif all([type(label) is str for label in which_beams]):
        new_beam_list = []
        for beam in beam_list:
            if beam.label in which_beams:
                new_beam_list.append(beam)
        beam_list = new_beam_list
    else:
        raise RuntimeError("Invalid argument plot_beams in plot_iv_from_csv()")

    # normalize each beam by the max intensity
    for beam in beam_list:
        maxint = max(beam.intens.values())
        beam.intens = {en: i / maxint for en, i in beam.intens.items()}

    plot_iv(beam_list, output_file, labels=[beam.label for beam in beam_list])
    return


def rot_mat_a(theta):
    """Generates a rotation matrix around the a axis.
    
    The rotation is positive, i.e. clockwise when looking along a.
    
    Parameters
    ----------
    theta : float
        Angle of rotation in degrees. (Use np.degrees to convert)
    """
    theta = np.radians(theta)
    rot_mat = [
        [1, 0, 0],
        [0, np.cos(theta), -np.sin(theta)],
        [0, np.sin(theta), np.cos(theta)],
    ]
    return np.array(rot_mat)


def rot_mat_c(theta):
    """Generates a rotation matrix around the c axis.
    
    The rotation is positive, i.e. clockwise when looking along c.
    
    Parameters
    ----------
    theta : float
        Angle of rotation in degrees. (Use np.degrees to convert)
    """
    theta = np.radians(theta)
    rot_mat = [
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1],
    ]
    return np.array(rot_mat)


# some other useful matrices

# Be careful: These WILL change chirality (handedness) of the cell -
# use rotation around a or c if you want to avoid this. (Chirality
# changes can not be expressed via rotations.)
switch_b_c_mat = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]])  # swap b and c
switch_a_b_mat = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])  # swap a and b

# flip the cell along c (useful if the lower half of the cell is to be kept)
flip_c_mat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
