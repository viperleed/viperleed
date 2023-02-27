# -*- coding: utf-8 -*-
"""A collection of functions that run ViPErLEED from a provided ASE object.

@Author: Alexander M. Imre, Michele Riva
based on tleedm_from_ase.py by Florian Kraushofer

Requires viperleed to be in the system path (or on PYTHONPATH).
"""

import copy
from collections import defaultdict
from io import StringIO
import logging
import os
from pathlib import Path
import shutil
from typing import Sequence
import warnings

import numpy as np

# for run_from_ase
import viperleed
from viperleed.tleedm import run_tleedm
from viperleed.tleedmlib.classes.slab import Slab
from viperleed.tleedmlib.files import poscar
from viperleed.tleedmlib.files.parameters import (readPARAMETERS,
                                                  interpretPARAMETERS)
# for rfactor_from_csv
from viperleed.tleedmlib.files.beams import readOUTBEAMS
from viperleed.tleedmlib.files.iorfactor import beamlist_to_array
from viperleed.tleedmlib.files.ivplot import plot_iv  # for plot_iv_from_csv
from viperleed.tleedmlib.wrapped.error_codes import check_ierr

try:
    from viperleed.tleedmlib.wrapped.rfactor import r_factor_new as rf          # TODO: pylint complains if not compiled
except ImportError:
    _HAS_NEW_RFACTOR = False
else:
    _HAS_NEW_RFACTOR = True

LOGGER = logging.getLogger()
_INPUT_FILES = (
    "PARAMETERS",
    "VIBROCC",
    "IVBEAMS",
    "PHASESHIFTS",
    "DISPLACEMENTS",
    "EXPBEAMS",
    "EXPBEAMS.csv",
    )

def run_from_ase(
    exec_path,
    ase_object,
    inputs_path=None,
    cut_cell_c_fraction=0.4,
    uc_transformation_matrix=None,
    uc_scaling=None,
    cleanup_work=False,
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
    cleanup_work : bool, optional
        Whether the work directory created during execution of the
        calculation is to be removed at the end. Default is False.

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
    FileNotFoundError
        If `exec_path` is non-existent or not a directory.
    FileNotFoundError
        If no PARAMETERS file was found in `exec_path`, after
        potentially copying over from `inputs_path`
    ValueError
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
    exec_path = Path(exec_path).resolve()
    if not exec_path.is_dir():
        raise FileNotFoundError(
            f"Invalid {exec_path=}: not "
            + "existent" if not exec_path.exists() else "a directory"
            )

    _copy_inputs_to_exec_path(inputs_path, exec_path)

    # Check for PARAMETERS file - without that we can't proceed
    # Files PHASESHIFTS and VIBROCC are not required and can be
    # generated if PARAMETERS are set correctly. IVBEAMS can be
    # generated if EXPBEAMS exists.
    if not (exec_path / "PARAMETERS").is_file():
        # No PARAMETERS file – Error
        raise FileNotFoundError("No PARAMETERS file found – this is required")

    # Transfer ASE object into slab object for ViPErLEED
    slab = Slab(ase_atoms=ase_object)

    # Get temporary parameters object
    rparams = readPARAMETERS(exec_path / parameters_name)
    interpretPARAMETERS(rparams, slab=slab)

    # Transformation of slab object: Rotation or isotropic streching/shrinking
    if uc_transformation_matrix is not None:
        slab.apply_matrix_transformation(uc_transformation_matrix)
    if uc_scaling is not None:
        slab.apply_scaling(*uc_scaling)

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

    _write_poscar(slab, exec_path)

    # Take care of input files and work directory
    # NOTE: PARAMETERS should NOT contain SITE_DEF flags.
    # VIBROCC should contain *_surf sites for all elements.
    work_path = _make_work_dir(exec_path)

    home = Path.cwd()
    os.chdir(work_path)

    # We are ready to run ViPErLEED! Have fun!
    try:
        run_tleedm(slab=slab,
                   preset_params=_make_preset_params(rparams, slab),
                   source=Path(viperleed.__file__).parent,)
    except Exception as err:
        # If ViPErLEED fails, move back to home directory
        os.chdir(home)
        raise RuntimeError("ViPErLEED calculation failed") from err

    # ViPErLEED should have suceeded if you arrive here. However, we may not
    # have run a refcalc (id == 1). In that case, return empty strings.
    if 1 not in rparams.RUN:
        return "", "", "", rparams.V0_IMAG

    # read out the THEOBEAMS.csv file and complex amplitudes
    content_list = _read_refcalc_output(rparams)

    # Move back home
    os.chdir(home)

    if cleanup_work:
        try:
            shutil.rmtree(work_path)
        except OSError as err:
            LOGGER.warning(
                f"Failed to remove work directory {work_path}. Info: {err}"
                )
    return *content_list, rparams.V0_IMAG


def _copy_inputs_to_exec_path(inputs_path, exec_path):
    """Copy all tleedm input files from inputs_path to exec_path."""
    if inputs_path is None:
        return
    inputs_path = Path(inputs_path)
    for file in _INPUT_FILES:
        try:
            shutil.copy2(inputs_path / file, exec_path / file)
        except FileNotFoundError:
            pass


def _make_preset_params(rparams, slab):
    """Return preset parameters to use when running tleedm."""
    if rparams.SITE_DEF:
        return {}

    # When no sites are explicitly defined, give all
    # atoms visible from vacuum a '_surf' site label
    LOGGER.warning("Parameter SITE_DEF not specified. Double check "
                   "PARAMETERS and POSCAR to make sure all sites are "
                   "recognized correctly. Default *_surf sites definitions "
                   "will be added for atoms visible from vacuum.")
    preset_params = {}
    site_def = defaultdict(list)
    for atom in slab.getSurfaceAtoms():
        site_def[atom.el].append(atom.oriN)
    preset_params["SITE_DEF"] = site_def

    return preset_params


def _make_work_dir(exec_path):
    """Return exec_path/'work' after copying there all input files."""
    work_path = exec_path / "work"
    work_path.mkdir(exist_ok=True)
    # copy input files to work directory
    for file in _INPUT_FILES:
        try:
            shutil.copy2(exec_path / file, work_path / file)
        except FileNotFoundError:
            pass
    return work_path


def _read_refcalc_output(rparams):
    """Return the contents of the output files produced by a refcalc."""
    # List of file name and earliest version in which it appeared
    output_files = (
        ("THEOBEAMS.csv", 0),
        ("Complex_amplitudes_real.csv", 1.73),
        ("Complex_amplitudes_imag.csv", 1.73),
        )

    content_list = []
    for filename, min_version in output_files:
        _path = Path(filename).resolve()
        content_str = ""
        if _path.is_file():
            with _path.open("r", encoding="utf-8") as fproxy:
                content_str = fproxy.read()
        elif rparams.TL_VERSION >= min_version:
            LOGGER.error(f"Could not find file {filename}")
        content_list.append(content_str)
    return content_list


def _write_poscar(slab, exec_path):                                             # TODO: could be done with a flag on writePOSCAR
    """Save slab as a POSCAR file, but warn if one is already there."""
    poscar_path = exec_path / "POSCAR"
    if poscar_path.exists():
        LOGGER.warning("A 'POSCAR' file is already present in "
                       f"{exec_path} and will be overwritten.")
    poscar.writePOSCAR(slab, poscar_path)


def rfactor_from_csv(
    beams_files,
    v0i,
    beams_file_is_content=(False, False),
    v0r_shift_range=(-3, 3),
    intpol_deg=5,
    intpol_step=0.5,
    return_beam_arrays = False,
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
    return_beam_arrays : bool, default = False
        For debugging, you can optionally return the interpolated beams,
        y-functions, number of overlapping points & per-beam R-factors.
        
    Returns
    -------
    best_R : float
        The R-factor between the two files. Returns the best R factor
        found after variation in v0r_shift_range.
    best_v0r : float
        The v0r value for which the best R-factor was found.
    """
    if not _HAS_NEW_RFACTOR:
        raise ModuleNotFoundError(
            "Missing R-factor compiled Fortran extension module. "
            "Run make in viperleed/tleedmlib/wrapped, then try again",
            name='viperleed.tleedmlib.wrapped.rfactor'
            )

    # Use Pendry R-factor - TODO: discuss if this should be a user
    # parameter in the API (if so, the code below will need a few
    # adjustments)
    which_r = 1

    # Check compilation of wrapped R-factor code
    check_ierr(rf.test_compilation())

    err_msg = ""

    # Check if file is supposed to be read from disk or from memory and pass
    # it to readOUTBEMS accordingly.
    beams_tmp = [None, None]
    for i, (is_content, file) in enumerate(zip(beams_file_is_content, beams_files)):
        if is_content:
            with StringIO(file) as fproxy:
                beams_tmp[i] = readOUTBEAMS(filename=fproxy, sep=",")
        else:
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

    averaging_scheme = np.int32(np.arange(n_beams) + 1)  # don't average
    n_derivs = 1

    # interpolate and prepare beam arrays

    # initialize beam arrays
    (
        beams1_e_start_beams_out,
        beams1_n_e_beams_out,
        beams1_intpol_intensity,
        beams1_deriv_y,
        ierr,
    ) = rf.alloc_beams_arrays(n_beams, len(grid))
    check_ierr(ierr)

    ierr = rf.prepare_beams(
        beams1_en,
        beams1_arr,
        beams1_id_start + 1,
        beams1_n_E_beams,
        intpol_deg,
        n_derivs,
        grid,
        beams1_e_start_beams_out,
        beams1_n_e_beams_out,
        beams1_intpol_intensity,
        beams1_deriv_y,
    )
    check_ierr(ierr)

    beams1_der = beams1_deriv_y.copy()

    # writes y into deriv_y
    rf.pendry_y_beamset(
        beams1_intpol_intensity,
        beams1_deriv_y,
        beams1_e_start_beams_out,
        beams1_n_e_beams_out,
        v0i
    )

    # Beams 2
    (
        beams2_e_start_beams_out,
        beams2_n_e_beams_out,
        beams2_intpol_intensity,
        beams2_deriv_y,
        ierr,
    ) = rf.alloc_beams_arrays(n_beams, len(grid))
    check_ierr(ierr)

    ierr = rf.prepare_beams(
        beams2_en,
        beams2_arr,
        beams2_id_start + 1,
        beams2_n_E_beams,
        intpol_deg,
        n_derivs,
        grid,
        beams2_e_start_beams_out,
        beams2_n_e_beams_out,
        beams2_intpol_intensity,
        beams2_deriv_y,

    )
    check_ierr(ierr)

    beams2_der = beams2_deriv_y.copy()

    # writes y into deriv_y
    rf.pendry_y_beamset(
        beams2_intpol_intensity,
        beams2_deriv_y,
        beams2_e_start_beams_out,
        beams2_n_e_beams_out,
        v0i
    )
    
    # Ready to calculate R-factor
    v0r_shift_range_int = np.int32([round(v / intpol_step) for v in v0r_shift_range])
    v0r_center = sum(v0r_shift_range_int) // 2
    start_guess = np.int32([v0r_shift_range_int[0], v0r_center, v0r_shift_range_int[1]])

    # Same settings as used in R-factor for ViPErLEED
    fast_search = False  # TODO change to True as default if we are sure it is better
    tol_r = (1 - 5e-2,)
    tol_r_2 = (1 - 5e-2,)
    max_fit_range = 6

    (_, best_v0r, best_r, n_evaluated, r_beams, numerators, denominators, n_overlapping_points, ierr) = rf.r_beamset_v0r_opt_on_grid(
        which_r,
        v0r_shift_range_int,
        start_guess,
        fast_search,
        intpol_step,
        beams1_deriv_y,
        beams2_deriv_y,
        beams1_e_start_beams_out,
        beams2_e_start_beams_out,
        beams1_n_e_beams_out,
        beams2_n_e_beams_out,
        tol_r=tol_r,
        tol_r_2=tol_r_2,
        max_fit_range=max_fit_range,
    )
    # If encountered an error
    check_ierr(ierr)

    #  for debugging, we have the option of returning the beam_arrays, y-Functions etc.
    if (return_beam_arrays):
        return best_r, best_v0r, (beams1_intpol_intensity, beams2_intpol_intensity), (beams1_deriv_y, beams2_deriv_y), r_beams, n_overlapping_points, (beams1, beams2)

    return best_r, best_v0r


def plot_iv_from_csv(
    beam_file,
    output_file=None,
    beam_file_is_content=False,
    which_beams=[],
    legends=None
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
    output_pdf : pathlib Path, string or None; default = None
        Path to the PDF file where the plot will be saved, if Path or string.
        (File is created.) If None, returns list of matplotlib figures.
    which_beams : range or Sequence of int or strings
        Indices specifying which beams to plot. Order is the same as in the csv.
    legends : None or list of strings, default = None
        Legends to be used for the various files.

    Raises
    ------
    RuntimeError
        If which_beams is not Sequence containing just integers or
        strings.
    """
    if isinstance(beam_file, Sequence):
        n_beams = len(beam_file)
        assert len(beam_file_is_content) == len(beam_file)
    else:
        n_beams = 1
        assert type(beam_file_is_content) is bool
    
    all_beam_data = []
    labels = []
    for file, is_content in zip(beam_file, beam_file_is_content):
        tmp_file = StringIO(file) if is_content else file
        beam_list = readOUTBEAMS(filename=tmp_file, sep=",")

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
            raise RuntimeError("Invalid argument which_beams in plot_iv_from_csv()")

        # normalize each beam by the max intensity
        for beam in beam_list:
            maxint = max(beam.intens.values())
            beam.intens = {en: i / maxint for en, i in beam.intens.items()}
        beam_labels=[beam.label for beam in beam_list]
        all_beam_data.append(beam_list)
        if not labels:
            labels = beam_labels
        else:
            for jj in range(len(beam_labels)):
                labels[jj] += ", " + beam_labels[jj]

    return plot_iv(all_beam_data, output_file, labels=labels, legends=legends)


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
