"""A collection of functions that run ViPErLEED from a provided ASE object."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexandra Mia Imre (@alexmiame)',
    'Michele Riva (@michele-riva)',
)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-06-04'
__license__ = 'GPLv3+'

from collections import defaultdict
from dataclasses import FrozenInstanceError, dataclass
from io import StringIO
import logging
from numbers import Real
from pathlib import Path
import shutil
from typing import Any, Sequence

import numpy as np

from viperleed.calc.classes.atom_containers import AtomList
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.rparams.special.energy_range import IVShiftRange
from viperleed.calc.classes.rparams.special.energy_range import TheoEnergies
from viperleed.calc.classes.slab import Slab
from viperleed.calc.constants import DEFAULT_WORK
from viperleed.calc.files import iorfactor as rf_io
from viperleed.calc.files import parameters, poscar
from viperleed.calc.files.beams import readOUTBEAMS
from viperleed.calc.files.ivplot import plot_iv
from viperleed.calc.files.iorfactor import beamlist_to_array
from viperleed.calc.files.iorfactor import write_Rfactorpdf_from_splines
from viperleed.calc.lib import rfactor as rfactor_lib
from viperleed.calc.lib import spline_interpolation
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.dataclass_utils import set_frozen_attr
from viperleed.calc.lib.matrix import rotation_matrix
from viperleed.calc.lib.rfactor.utils import average_beam_array
from viperleed.calc.run import run_calc


_LOGGER = logging.getLogger(__name__)
_INPUT_FILES = (
    "PARAMETERS",
    "VIBROCC",
    "IVBEAMS",
    "PHASESHIFTS",
    "DISPLACEMENTS",
    "EXPBEAMS",
    "EXPBEAMS.csv",
    )


@dataclass
class SlabTransform:
    """A container for unit-cell transformations.

    Attributes
    ----------
    orthogonal_matrix : numpy.ndarray, optional
        Shape (3, 3). Can contain an arbitrary orthogonal transform
        (i.e., a combination of mirroring and rotations) that will
        be applied to the cell. The transformation is applied BEFORE
        scaling, i.e., using the original unit cell. The matrix (O)
        is applied to both the unit cell and all Cartesian atomic
        positions. The unit cell U (unit vectors as columns) is
        transformed to O @ U. Atom Cartesian coordinates (v) are
        transformed to O @ v. This attribute can be used, e.g., to
        apply a rotation or a mirroring of the slab.
    uc_scaling : float or Sequence of float, optional
        Scaling to be applied to the cell. Scaling is applied AFTER
        the transformation, i.e., along the new unit cell vectors. This
        can be used to apply an isotropic scaling in all directions or
        to stretch/compress along unit cell vectors in order to change
        lattice constants in some direction. If the scaling is given as
        a number or a 1-element sequence, an isotropic scaling will be
        applied to the unit cell and atom positions. If a sequence with
        three entries is provided, the scaling will be applied along
        the unit cell vectors in the given order. Specifically, given
        uc_scaling == (s1, s2, s3), the new unit cell vectors will be
        a' = a * s1, b' = b * s2, c' = c * s3. Default is None, i.e.,
        no scaling.
    cut_cell_c_fraction : float, optional
        Fractional position along the c vector to cut the slab.
        (ViPErLEED expects the cell to have the surface at the
        "top" -- i.e., at high z coordinates, and bulk-like layers
        at the bottom). Cutting will be performed along the vector
        c AFTER all transformations (see uc_transformation_matrix)
        are applied. Only the part of the cell above (i.e., >=)
        this value will be kept. Everything else is discarded.
        Default is zero.
    """

    orthogonal_matrix: Any = None
    uc_scaling: Any = None
    cut_cell_c_fraction: float = 0.


# We would really like to use frozen=True here, but since python 3.8.10
# a frozen cannot inherit from a non-frozen. We're forced to emulate
# the frozen-ness by overriding __setattr__ and __delattr__, and using
# the __init__ implementation for frozen dataclasses. It is NOT A GOOD
# IDEA TO INHERIT FROM THIS WORKAROUND CLASS.
@dataclass
class _ImmutableSlabTransform(SlabTransform):
    """A container for immutable slab transformations.

    This class is meant exclusively as the type for Slab
    transformations that cannot be represented via a left-
    multiplication of the unit cell with an orthogonal
    matrix. These transforms cannot perform any scaling
    nor cutting of the unit cell. They should be combined
    with SlabTransform instances for that purpose.
    """

    def __init__(self):
        """Initialize instance attributes."""
        set_frozen_attr(self, 'orthogonal_matrix', None)
        set_frozen_attr(self, 'uc_scaling', None)
        set_frozen_attr(self, 'cut_cell_c_fraction', 0.)

    def __setattr__(self, name, value):
        """Raise FrozenInstanceError."""
        if name == '__doc__':
            # One exception to allow nicer 'help' for instances
            super().__setattr__(name, value)
            return
        raise FrozenInstanceError(f"cannot assign to field {name!r}")

    def __delattr__(self, name):
        """Raise FrozenInstanceError."""
        raise FrozenInstanceError(f"cannot assign to field {name!r}")


def run_from_ase(exec_path, ase_object, inputs_path=None,
                 slab_transforms=(SlabTransform(),), cleanup_work=False):
    """Perform a ViPErLEED calculation starting from an ase.Atoms object.

    This function is a core part of the API for ViPErLEED. It allows
    to run a calculation from an ASE (Atomic Simulation Environment)
    Atoms object. The ase.Atoms object, which is expected to contain
    a simulation cell, will be cut at a specified height fraction,
    transformed into a ViPErLEED Slab and written to a POSCAR file. The
    ViPErLEED calculation is then executed on the structure using the
    input files found either in `inputs_path` or `exec_path`. The type
    of calculation performed depends on the value of the RUN parameter
    in the PARAMETERS file. If RUN = 1, this yields reference electron
    scattering amplitudes and intensities. These are collected in CSV
    files, which are returned as strings for further processing.

    Parameters
    ----------
    exec_path: string or Path
        Path where to execute the ViPErLEED calculation. The directory
        needs to exist. Input files will be copied from `inputs_path`,
        if provided. Otherwise it is assumed that all input files are
        already in `exec_path`. Temporary files will be stored in a
        subdirectory called work.
    ase_object : ase.Atoms
        ASE-type object that contains the atoms for the slab to be used
        with ViPErLEED. This is expected to contain a simulation cell,
        which will be transformed using `slab_transforms`.
    inputs_path: str or Path or None, optional
        Path to directory containing input files (PARAMETERS, VIBROCC,
        etc.) for ViPErLEED calculation. If provided, files will be
        copied from `inputs_path` to `exec_path` (silently overwriting
        existing files there!). If None or not given, no copying
        occurs. Default is None.
    slab_transforms : SlabTransform or Sequence thereof, optional
        Sequence of transformation operations to be applied to the Slab
        constructed from ase_object. These include optional "change of
        basis" operations, optional `uc_scaling` for scaling the unit
        vector lengths, and a `cut_cell_c_fraction` for selecting only
        a portion of the `ase_object`. In addition, some special
        SlabTransform objects that swap unit cell axes can also be used
        (swap_a_b, swap_b_c, and swap_c_a). See help(SlabTransform) and
        help(swap_*_*) for further details. The transforms are applied
        in the order given. The `cut_cell_c_fraction` is taken only
        from the LAST SlabTransform given. Make sure that vectors a & b
        do not have any components in z after transformation as this is
        not allowed and will raise a ValueError. Default is a single
        SlabTransform(), corresponding to no transformation and no cut.
    cleanup_work : bool, optional
        Whether the work directory created during execution of the
        calculation is to be removed at the end. Default is False.

    Returns
    -------
    theobeams_file_str : str
        Contents of the file "THEOBEAMS.csv", i.e. the scattering
        intensities as calculated in the reference calculation. An
        empty string is returned if no reference calculation was run.
    amp_real_file_str : str
        Contents of the file "Complex_amplitudes_real.csv", i.e.,
        the real part of the scattering amplitudes as calculated
        in the reference calculation. An empty string is returned
        if no reference calculation was run or if the file was not
        found.
    amp_imag_file_str : str
        Contents of the file "Complex_amplitudes_imag.csv", i.e., the
        imaginary part of the scattering amplitudes as calculated in
        the reference calculation. An empty string is returned if no
        reference calculation was run or if the file was not found.
    v0i : float
        Imaginary part of the inner potential as read from the
        PARAMETERS file. This is returned here since it is an input
        parameter for the R-factor calculation. The value is read
        from PARAMETERS and left unchanged - it is not a result of
        the calculation and is returned here for convenience only.

    Raises
    ------
    FileNotFoundError
        If `exec_path` is non-existent or not a directory.
    FileNotFoundError
        If no PARAMETERS file was found in `exec_path`, after
        potentially copying over from `inputs_path`
    ValueError
        If either of the first two unit cell vectors have a
        component perpendicular to the surface (i.e., along
        the z coordinate) after transformation of the unit cell
    RuntimeError
        If the ViPErLEED calculation fails.

    Notes
    -------
    The parameter slab_transform and its `orthogonal_matrix` attribute
    can be used to apply some useful transformations to the Slab.
    Examples include rotations around an axis and flipping the cell
    along c (in case the lower part has to be kept, rather than the
    upper one). A few special SlabTransform objects are also provided
    for swapping unit-cell vectors. Have a look at the matrices
    provided as part of this module (use as the `orthogonal_matrix`
    attribute of a SlabTransform):
    -  rot_mat_x(theta) : Rotation matrix around x by theta (deg)
    -  rot_mat_z(theta) : Rotation matrix around z by theta (deg)
    -  rot_mat_axis(axis, theta): Rotation around axis by theta (deg)
    -  flip_c_mat : Mirror matrix that flips the cell along c
    and the special SlabTransform objects:
    -  swap_a_b : Swap vectors a & b (changes handedness)
    -  swap_b_c : Swap vectors b & c (changes handedness)
    -  swap_c_a : Swap vectors c & a (changes handedness)
    For applying multiple operations, the order does matter. You
    can: (i) Combine multiple orthogonal transformation matrices
    via matrix multiplication (@ symbol or np.dot) into a single
    SlabTransform object (into its `orthogonal_matrix` attribute):
    in this case the leftmost matrix is applied last. (ii) Provide
    a sequence of SlabTransform objects in `slab_transforms`: here
    transformations are applied in the order given (i.e., the
    leftmost is the applied first). (iii) A combination of the
    previous two.
    """
    exec_path = Path(exec_path).resolve()
    if not exec_path.is_dir():
        raise FileNotFoundError(
            f"Invalid exec_path={exec_path}: not "
            + "existent" if not exec_path.exists() else "a directory"
            )

    _copy_inputs_to_exec_path(inputs_path, exec_path)

    # Check for PARAMETERS file - without that we can't proceed
    # Files PHASESHIFTS and VIBROCC are not required and can be
    # generated if PARAMETERS are set correctly. IVBEAMS can be
    # generated if EXPBEAMS exists.
    parameters_file = exec_path / "PARAMETERS"
    if not parameters_file.is_file():
        # No PARAMETERS file – Error
        raise FileNotFoundError("No PARAMETERS file found – this is required")

    # Transfer ASE object into slab object for ViPErLEED,
    # and save it as a POSCAR file for debug purposes.
    slab = _make_and_check_slab(ase_object, slab_transforms)
    _write_poscar(slab, exec_path)

    # Take care of input files and work directory
    # NOTE: PARAMETERS should NOT contain SITE_DEF flags.
    # VIBROCC should contain *_surf sites for all elements.
    work_path = _make_work_dir(exec_path)

    home = Path.cwd()
    with execute_in_dir(work_path):
        # Get temporary parameters object
        rpars = parameters.read(parameters_file)
        parameters.interpret(rpars, slab=slab, silent=False)

        # We are ready to run ViPErLEED! Have fun!
        try:
            exit_code, _ = run_calc(
                slab=slab,
                preset_params=_make_preset_params(rpars, slab),
                home=home,
                )
        except Exception as exc:
            raise RuntimeError('ViPErLEED calculation failed') from exc

        if exit_code:
            raise RuntimeError('ViPErLEED calculation failed. '
                               'See log file for details.')

        # ViPErLEED should have succeeded if you arrive here. However,
        # we may not have run a refcalc (id == 1). In that case, return
        # empty strings.
        if 1 not in rpars.RUN:
            return '', '', '', rpars.V0_IMAG

        # Read out the THEOBEAMS.csv file and complex amplitudes
        content_list = _read_refcalc_output(rpars)

    if cleanup_work:
        try:
            shutil.rmtree(work_path)
        except OSError as exc:
            _LOGGER.warning('Failed to remove work directory '
                            f'{work_path}. Info: {exc}')
    return (*content_list, rpars.V0_IMAG)


def _copy_inputs_to_exec_path(inputs_path, exec_path):
    """Copy all calc input files from inputs_path to exec_path."""
    if inputs_path is None:
        return
    inputs_path = Path(inputs_path)
    for file in _INPUT_FILES:
        try:
            shutil.copy2(inputs_path / file, exec_path / file)
        except FileNotFoundError:
            pass


def _apply_transform(slab, transform, apply_cut=False):
    """Apply a single SlabTransform `transform` to `slab`."""
    indices = None       # Unit-cell and position indices to be swapped
    if transform is swap_a_b:
        indices = (0, 1), (1, 0)
    elif transform is swap_b_c:
        indices = (2, 1), (1, 2)
    elif transform is swap_c_a:
        indices = (0, 2), (2, 0)

    if indices:
        # Swapping axes will require to make a new bulk slab for sure
        if not slab.is_bulk:
            slab.bulkslab = None
        new_axes, old_axes = ((ind,) for ind in indices)
        slab.ucell.T[new_axes] = slab.ucell.T[old_axes]
        for atom in slab:
            atom.pos[new_axes] = atom.pos[old_axes]
        return  # No scaling nor cutting

    if transform.orthogonal_matrix is not None:
        slab.apply_matrix_transformation(transform.orthogonal_matrix)
    if transform.uc_scaling is not None:
        if isinstance(transform.uc_scaling, Real):
            transform.uc_scaling = (transform.uc_scaling,)
        slab.apply_scaling(*transform.uc_scaling)

    if not apply_cut or not transform.cut_cell_c_fraction:
        return

    # Cut the slab as given in cut_cell_c_fraction - very important
    if transform.cut_cell_c_fraction:
        slab.atlist = AtomList(at for at in slab
                               if at.pos[2] >= transform.cut_cell_c_fraction)
    slab.update_atom_numbers()
    slab.update_element_count()


def _make_and_check_slab(ase_object, transforms):
    """Return a Slab from ase.Atoms with transformations applied."""
    slab = Slab.from_ase(ase_object)

    # Transformations of slab:
    # Mirror/rotation/swapping and/or stretching/shrinking
    if isinstance(transforms, SlabTransform):
        transforms = (transforms,)
    last_transform = transforms[-1]
    for transform in transforms:
        _apply_transform(slab, transform,
                         apply_cut=transform is last_transform)

    # Check if the now transformed slab has any z components in vectors
    # a & b. Error out as this would mess up parts of the TensErLEED
    # calculations, as we always pass on only the first two components.
    ab_vecs = slab.ucell.T[:2]
    if any(ab_vecs[:, 2] > 1e-5):
        raise ValueError(
            f"z component found in unit cell vector a ({ab_vecs[0]}) or b "
            f"({ab_vecs[1]}). This is not allowed in a TensErLEED calculation."
            " Check eventual transformations applied to the unit cell and make"
            " sure a and b are in the (x,y) plane."
            )
    return slab


def _make_preset_params(rparams, slab):
    """Return preset parameters to use when running viperleed.calc."""
    if rparams.SITE_DEF:
        return {}

    # When no sites are explicitly defined, give all
    # atoms visible from vacuum a '_surf' site label
    _LOGGER.warning("Parameter SITE_DEF not specified. Double check "
                    "PARAMETERS and POSCAR to make sure all sites are "
                    "recognized correctly. Default *_surf sites definitions "
                    "will be added for atoms visible from vacuum.")
    preset_params = {}
    site_def = defaultdict(list)
    for atom in slab.getSurfaceAtoms(rparams):
        site_def[atom.el].append(atom.num)
    preset_params["SITE_DEF"] = {
        element: {"surf": atom_numbers}
        for element, atom_numbers in site_def.items()
        }

    return preset_params


def _make_work_dir(exec_path):
    """Return exec_path/'work' after copying there all input files."""
    work_path = exec_path / DEFAULT_WORK
    work_path.mkdir(exist_ok=True)
    # copy input files to work directory
    for file in _INPUT_FILES:
        try:
            shutil.copy2(exec_path / file, work_path / file)
        except FileNotFoundError:
            pass
    return work_path


def _read_refcalc_output(rpars):
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
        elif rpars.TL_VERSION >= min_version:
            _LOGGER.error(f"Could not find file {filename}")
        content_list.append(content_str)
    return content_list


def _write_poscar(slab, exec_path):                                             # TODO: could be done with a flag on poscar.write
    """Save slab as a POSCAR file, but warn if one is already there."""
    poscar_path = exec_path / "POSCAR"
    if poscar_path.exists():
        _LOGGER.warning("A 'POSCAR' file is already present in "
                        f"{exec_path} and will be overwritten.")
    poscar.write(slab, poscar_path)

def _match_beams_by_label(expbeams, theobeams):
    """
    Return matched exp/theory beam lists using direct hkfrac matching.
    """
    theo_by_label = {beam.hkfrac: beam for beam in theobeams}
    matched_exp = []
    matched_theo = []
    missing_exp = []
    for beam in expbeams:
        mate = theo_by_label.get(beam.hkfrac)
        if mate is None:
            missing_exp.append(beam.label)
            continue
        matched_exp.append(beam)
        matched_theo.append(mate)
    exp_labels = {beam.hkfrac for beam in expbeams}
    missing_theo = [beam.label for beam in theobeams
                    if beam.hkfrac not in exp_labels]
    if missing_exp or missing_theo:
        parts = []
        if missing_exp:
            parts.append(
                "Experimental beams without theoretical match: "
                + ", ".join(missing_exp)
            )
        if missing_theo:
            parts.append(
                "Theoretical beams without experimental match: "
                + ", ".join(missing_theo)
            )
        _LOGGER.warning("\n".join(parts))
    if not matched_exp:
        raise ValueError("No common beams found between experimental and"
                         "theoretical data.")
    return matched_exp, matched_theo


def rfactor_from_csv(
    beams_files,
    v0i,
    beams_file_is_content=(False, False),
    v0r_shift_range=(-3.0, 3.0),
    intpol_step=0.5,
    which_r='pendry',
    beam_correspondence=None,
    plot_iv=None,
    intpol_deg=None,
    return_beam_arrays=None,
):
    """
    Compute the Pendry R-factor between two CSV files (in ViPErLEED format).

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
        a path to the file. If theoretical beams are to be averaged, then the
        first file should be experimental, the second theoretical data.
    v0i : float
        Imaginary part of the inner potential in eV.
    beams_file_is_content : tuple of two bool, default=(False, False)
        If either value is set to True, it is assumed that the
        corresponding element of beams_file contains a string with
        the CSV contents rather than a path to the file. For reading,
        a StringIO will be used internally.
    v0r_shift_range : tuple of two float, default=(-3, 3)
        During R-factor calculation, the two beam sets can be shifted
        in energy relative to one another. v0r_shift_range determines
        the minimum and maximum shifts (given in eV).
    intpol_step : float, default: 0.5
        Energy step for the interpolation/output grid.
    which_r : str, default: 'pendry'
        Selects which R-factor to use (e.g. 'pendry'/'smooth').
    beam_correspondence : sequence, optional
        If given, use this correspondence directly in `average_beam_array`
        format. In that case, all experimental beams are kept as read from the
        first file, and the theoretical beams from the second file are averaged/
        reordered according to this correspondence.
        If omitted, beams are matched only by direct label (`hkfrac`) equality.
    plot_iv : dict, optional
        If given, create a dummy Rparams object, update its PLOT_IV dict with
        contents of plot_iv, and call `write_Rfactorpdf_from_splines(...)`.
        To use defaults, pass: plot_iv={'plot': True}
    intpol_deg : ignored, deprecated
        Deprecated compatibility shim. Only cubic spline interpolation is used.
    return_beam_arrays : ignored, deprecated
        Deprecated compatibility shim. The old Fortran-array debug return path
        is not supported here.

    Returns
    -------
    best_r : float
        Best overall R-factor.
    best_v0r : float
        The V0r shift in eV corresponding to best_r.
    r_beams : np.ndarray
        Per-beam R-factors at best_v0r.
    """
    if intpol_deg is not None:
        _LOGGER.warning(
            "intpol_deg is deprecated in rfactor_from_ase(); "
            "the spline-based implementation always uses cubic splines."
        )

    if return_beam_arrays not in (None, False):
        _LOGGER.warning(
            "return_beam_arrays is deprecated in rfactor_from_ase(); "
            "the spline-based implementation no longer returns the old "
            "Fortran interpolation arrays."
        )

    # Select R-factor function
    r_func = rfactor_lib.select_rfactor(which_r)

    # Check if file is supposed to be read from disk or from memory and pass
    # it to readOUTBEMS accordingly.
    beams_tmp = [None, None]
    for i, (is_content, file) in enumerate(zip(beams_file_is_content,
                                               beams_files)):
        if is_content:
            with StringIO(file) as fproxy:
                beams_tmp[i] = readOUTBEAMS(filename=fproxy, sep=",")
        else:
            beams_tmp[i] = readOUTBEAMS(filename=file, sep=",")

    # beams are here named as expbeams and theobeams for convenience, as this
    #  is assumed if they are to be averaged. However, they could also be two
    #  sets of theoretical or experimental beams.
    expbeams, theobeams = beams_tmp

    for i, beams in enumerate((expbeams, theobeams)):
        unique_lables = set()
        for beam in beams:
            unique_lables.add(beam.hkfrac)
        if len(unique_lables) != len(beams):
            raise RuntimeError(
                f"Multiple beams sharing the same index encountered"
                f" in input beamset {i}!"
            )

    if beam_correspondence is None:
        expbeams, theobeams = _match_beams_by_label(expbeams, theobeams)
        beam_correspondence = np.arange(len(theobeams), dtype=int)
    else:
        beam_correspondence = np.asarray(beam_correspondence, dtype=int)
        if beam_correspondence.ndim != 1:
            raise ValueError("beam_correspondence must be a 1D sequence of "
                             "integers.")
        if len(beam_correspondence) != len(theobeams):
            raise ValueError(
                "beam_correspondence must have one entry per theoretical beam "
                f"({len(theobeams)}), got {len(beam_correspondence)}."
            )

    # Prepare an Rparams object to get the energy ranges right.
    # Again treat first beam as 'experiment' and second as 'theory'.
    rpars = Rparams()
    rpars.expbeams = expbeams
    rpars.V0_IMAG = v0i

    theo_grid, _, _, theo_intensities = beamlist_to_array(theobeams)
    rpars.THEO_ENERGIES = TheoEnergies.from_sorted_grid(theo_grid)
    rpars.IV_SHIFT_RANGE = IVShiftRange(*v0r_shift_range, intpol_step)

    _, theo_range, _, intpol_step = rf_io.prepare_rfactor_energy_ranges(
        rpars,
        theobeams,
        for_error=False,
        n_expand=1,   # cubic only
    )
    out_grid = np.arange(
        theo_range.min,
        theo_range.max + 0.1 * intpol_step,
        intpol_step,
    )

    # Create splines
    exp_grid, _, _, exp_intensities = beamlist_to_array(expbeams)
    exp_spline = spline_interpolation.interpolate_ragged_array(
        x=exp_grid,
        y=exp_intensities,
        bc_type="not-a-knot",
    )
    exp_spline = spline_interpolation.CachedSpline(exp_spline)

    theo_intensities = average_beam_array(theo_intensities, beam_correspondence)
    theo_spline = spline_interpolation.interpolate_ragged_array(
        x=theo_grid,
        y=theo_intensities,
        bc_type="not-a-knot",
    )
    theo_spline = spline_interpolation.CachedSpline(theo_spline)

    # sweep V0r shifts - TODO: continuous shifts?
    shifts = np.arange(
        rpars.IV_SHIFT_RANGE.start,
        rpars.IV_SHIFT_RANGE.stop + rpars.IV_SHIFT_RANGE.step,
        rpars.IV_SHIFT_RANGE.step,
    )

    r_values = np.array(
        [
            r_func(
                v0i,
                intpol_step,
                out_grid,
                data_spline_1=exp_spline,
                data_spline_2=theo_spline,
                shift_2nd_spline=shift,
                groups=None,
            )
            for shift in shifts
        ]
    )

    best_index = np.argmin(r_values)
    best_v0r = shifts[best_index]
    best_r = r_values[best_index]

    # re-build theoretical splines one more time on the shifted grid
    if best_v0r != 0:
        theo_spline = spline_interpolation.interpolate_ragged_array(
            x=theo_grid + best_v0r,
            y=theo_intensities,
            bc_type="not-a-knot",
        )
        theo_spline = spline_interpolation.CachedSpline(theo_spline)

    r_beams = r_func(
        v0i,
        intpol_step,
        out_grid,
        data_spline_1=exp_spline,
        data_spline_2=theo_spline,
        groups="beam",
    )

    if plot_iv is not None:
        rpars.PLOT_IV.update(plot_iv)
        rpars.expbeams = expbeams

        write_Rfactorpdf_from_splines(
            rpars,
            theo_spline,
            exp_spline,
            out_grid,
            r_beams,
            outname="Rfactor_plots_from_ase.pdf",
            analysis_file="Rfactor_analysis_from_ase.pdf",
        )

    return best_r, best_v0r, r_beams

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
        assert len(beam_file_is_content) == len(beam_file)
    else:
        assert type(beam_file_is_content) is bool

    all_beam_data = []
    labels = []
    for file, is_content in zip(beam_file, beam_file_is_content):
        tmp_file = StringIO(file) if is_content else file
        beam_list = readOUTBEAMS(filename=tmp_file, sep=",")
        if not beam_list:
            continue

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

    if not all_beam_data:
        # Nothing to plot
        raise ValueError("No beams to plot in beam_file")

    return plot_iv(all_beam_data, output_file, labels=labels, legends=legends)


def rot_mat_x(theta):
    """Return a rotation matrix around the x axis.

    The rotation is positive, i.e., clockwise when looking along x,
    when applied from the left to column vectors. The same rotation
    for row vectors can be obtained by multiplying on the right with
    the transpose of the return value.

    Parameters
    ----------
    theta : float
        Angle of rotation in degrees.

    Returns
    -------
    rot_mat : numpy.ndarray
        Shape (3, 3). Rotation matrix around the x axis.
    """
    theta = np.radians(theta)
    rot_mat = np.identity(3)
    rot_mat[1:, 1:] = rotation_matrix(theta, dim=2)
    return rot_mat


def rot_mat_z(theta):
    """Return a rotation matrix around the z axis.

    The rotation is positive, i.e. clockwise when looking along z,
    when applied from the left to column vectors. The same rotation
    for row vectors can be obtained by multiplying on the right with
    the transpose of the return value.

    Parameters
    ----------
    theta : float
        Angle of rotation in degrees.

    Returns
    -------
    rot_mat : numpy.ndarray
        Shape (3, 3). Rotation matrix around the z axis.
    """
    theta = np.radians(theta)
    return rotation_matrix(theta, dim=3)


def rot_mat_axis(axis, theta):
    """Return a 3D rotation matrix by `theta` around `axis`.

    Parameters
    ----------
    axis : numpy.ndarray
        A vector parallel to the rotation axis. Shape (3,).
    theta : float
        Rotation angle in degrees.

    Returns
    -------
    rotation_matrix : numpy.ndarray
        Rotation matrix around `axis` by `theta`. The rotation is
        counter-clockwise when applied to column vectors from
        the left.
    """
    # We use the 'concise form' of the axis--angle formulation from
    # en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    axis_dir = axis / np.linalg.norm(axis)

    # Cross-product matrix from stackoverflow.com/questions/66707295
    u_cross = np.cross(axis_dir, -np.identity(3))
    u_outer = np.outer(axis_dir, axis_dir)
    theta = np.radians(theta)
    return (np.cos(theta) * np.identity(3)
            + np.sin(theta) * u_cross
            + (1 - np.cos(theta)) * u_outer)


# Some other useful transformations
# Swapping unit cell vectors
swap_a_b = swap_b_a = _ImmutableSlabTransform()
swap_b_c = swap_c_b = _ImmutableSlabTransform()
swap_c_a = swap_a_c = _ImmutableSlabTransform()


# Notice: the __doc__ edit will work fine in python >=3.10, but earlier
# versions print out type(obj).__doc__ when help is called. There's no
# good way around this, other than making each a separate class.
_FMT_DOC = """
    This SlabTransform swaps the '{}' and '{}' unit cell vectors,
    but does not change the Cartesian coordinates of the atoms.

    Be careful, as this transform WILL CHANGE THE HANDEDNESS of the
    unit cell. To avoid this, you can instead use a rotation around an
    appropriate axis. For orthogonal cells, you can use a 90 degrees
    rotation around '{}'.
"""

# pylint: disable=no-member
# Perhaps a bug? __doc__ exists and is set properly
swap_a_b.__doc__ += _FMT_DOC.format('a', 'b', 'c')
swap_b_c.__doc__ += _FMT_DOC.format('b', 'c', 'a')
swap_c_a.__doc__ += _FMT_DOC.format('c', 'a', 'b')
# pylint: enable=no-member


# flip the cell along c (useful if the lower half of the cell is to be kept)
flip_c_mat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
