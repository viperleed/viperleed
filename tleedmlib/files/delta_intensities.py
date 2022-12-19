# -*- coding: utf-8 -*-
"""

@author: Tobias Hable, Alexander M. Imre

Reads in delta files.
"""
import sys
import numpy as np
from numpy import sin, cos, sqrt
import fortranformat as ff
import matplotlib.pyplot as plt
import scipy
import os
import cmath # just used for sqrt
import warnings

from tqdm import tqdm # progress bar; optional
from numba import njit, prange, config, threading_layer


def read_delta_file(filename, n_energies, read_header_only=False):
    """Read and return the contents of a 
    TensErLEED delta-amplitude file.
    
    This function reads in one file of data and stores the data 
    in arrays, which can be used in later functions 
    (ex.: calc_delta_intensities).
    
    Parameters
    ----------
    filename : str
        The filename describes which file you want to read in 
        (path from the function location to the file)  
    n_energies : int
        Number of different energies. This should be a 
        known factor for your files.
        Also possible to just loop the file until the end.
    read_header_only : bool
        If True does reads Header only and stops afterwards
    
    Returns
    -------
    (phi, theta): tuple of float
        Angles of how the beam hits the sample
    (trar1, trar2) : tuple of ndarray
        reciprocal lattice vectors
    n_beams : int
        Number of beams for which delta-amplitudes were present
        in the file 
    nc_steps: numpy.ndarray
        Number of permutations between direction deltas and
        vibration deltas  
    e_kin_array : numpy.ndarray
        Array that contains the kinetic energies of the
        elastically scattered electrons inside the crystal.
        shape=(n_energies)
        (Not the incidence energy!)
    v_imag_array : numpy.ndarray
        Imaginary part of the inner potential of the surface
        shape=(n_energies)
    VV_array : numpy.ndarray
        Real part of the inner potential of the surface
        shape=(n_energies)
    beam_indices : numpy.ndarray
        Array of beam indices, with beam_indices[i] == [h_i, k_i];
        shape=(n_beams,2)
    Cundisp : numpy.ndarray
        Position of the undisplaced atoms (always 0)
        shape=(3)
    geo_delta : numpy.ndarray
        Geometric displacement of given delta
        shape=(n_geo_vib_grid, 3);
        n_geo_vib_grid read from header_block_2
    amplitudes_ref : numpy.ndarray
        Array that contains all values of the reference amplitudes
        shape=(n_energies,n_beamsf)
    amplitudes_del : numpy.ndarray
        Array that contains all values of the delta amplitudes
        shape=(n_energies,n_vib,n_geo,n_beams)
        n_vib and n_geo read in header_block_6;
        they are the number of geometric and vibrational displacements
    """

    header_block_1 = []
    header_block_2 = []
    e_kin_array = np.full(n_energies, fill_value=np.nan)
    v_inner_array = np.full(n_energies, fill_value=np.nan, 
                            dtype=np.complex128)

    # we need three fortran format readers
    ff_reader_6E13_7 = ff.FortranRecordReader("6E13.7")
    ff_reader_10F10_5 = ff.FortranRecordReader("10F10.5")
    ff_reader_10F7_4 = ff.FortranRecordReader("10F7.4")

    # Reading in the data of a file
    try:
        with open(filename, mode="r") as file:
            content = file.readlines()
    except Exception as err:
        warnings.warn(f"Unable to read Delta file: {filename}")
        raise err

    if not content:
        raise ValueError(f"File {filename} is empty.")

    if len(content) < 2:
        raise ValueError(f"Invalid delta file {filename}. "
                         "Not enough lines. "
                         f"Found {len(contents)}, "
                         "expected at least 2.")
    # make into an iterator
    file_lines = iter(content)

    # 1.Block of Header - only 1 line - theta, phi, trar1, 
    #trar2 variables
    line = next(file_lines)
    header_block_1 = ff_reader_6E13_7.read(line)
    
    # surface unit cell vectors - what used to be 
    # trar1 is now trar[:,0].
    # Similarly trar2 -> trar[:,1]
    # [check if you also need a .T to keep the convention 
    # that trar[:,0] =  trar1
    # Also, it makes much more sense to store unit-cell vectors 
    # the other way around, such that trar1 = trar[0] (== trar[0, :]).
    # Usually makes many of the calculations easier.]
    theta, phi, *trar = header_block_1
    trar = np.array(trar).reshape(2,2)

    # 2.Block of Header - also only 1 line - n_beams, n_atoms,
    #n_geo_vib_grid variables
    line = next(file_lines)
    header_block_2 = [int(p) for p in line.split()]

    if len(header_block_2) == 2:
        n_beams, n_geo_vib_grid = header_block_2
    elif len(header_block_2) == 3:
        n_beams, n_atoms, n_geo_vib_grid = header_block_2
        if (n_atoms != 1):
            raise NotImplementedError(
                    f"Unsupported delta-amplitude file {filename}. "
                    f"Found NATOMS={n_atoms}, but only NATOMS=1 "
                    "is supported."
            )           

    else:
        raise ValueError(f"Invalid header in file {filename}: "
                         "second line should contain 2 or 3 elements. "
                         f"Found (len{header_block_2}")
    

    # 3.Block of Header - (h,k) indices of the beams
    beam_indices = read_block(reader=ff_reader_10F10_5, 
                              lines=file_lines, shape=(n_beams, 2))

    # TODO: if we decide to throw CUNDISP out of TensErLEED entirely, 
    # this block needs to become optional
    # 4.Block of Header - position of undisplaced atom 
    # (Coordinates UNDISPlaced)
    # Unused quantity - only check if it is zero (as it should be)
    pos_undisplaced = read_block(reader=ff_reader_10F7_4, 
                                 lines=file_lines, shape=(3,))
    if (np.linalg.norm(pos_undisplaced) > 1e-6):
        raise NotImplementedError(
                "A non-zero value of CUNDISP (undisplaced atom "
                f"positions) was read from Delta file {filename}. "
                "This quantity is currently unused and should always "
                "be zero. Rerun Refcalc and Delta calulation with a "
                "newer TensErLEED version."
        )

    # 5.Block of Header - geometric displacements 
    #(Coordinates DISPlaced)
    # For now, this contains, along the first axis, 
    #n_vib repetitions of the same
    # displacements. We will figure out n_vib further below, 
    #then reshape this
    geo_delta = read_block(reader=ff_reader_10F7_4, 
                           lines=file_lines, shape=(n_geo_vib_grid, 3))

    # 6.Block of Header - list of (vib, 0,0,0,0,..., vib, 0,0,0,0,...)
    vib_delta = read_block(reader=ff_reader_10F7_4, 
                           lines=file_lines, shape=(n_geo_vib_grid,))
    n_vib = sum(abs(v)>1e-4 for v in vib_delta)
    n_geo = n_geo_vib_grid // n_vib
    assert n_geo_vib_grid % n_vib == 0
    geo_delta = geo_delta.reshape(n_vib, n_geo, 3)[0, :, :].reshape(n_geo,3)
    # throw out the zeros from array vib_delta
    vib_delta = vib_delta[::n_geo]

    if read_header_only:
        return (
        (phi, theta),
        trar,
        (n_beams, n_geo, n_vib),
        beam_indices,
        geo_delta,
        vib_delta,
        e_kin_array,
        v_inner_array,
        None,
        None
        )

    # Initialize arrays for reference and delta amplitudes
    amplitudes_ref = np.full(shape=(n_energies, n_beams),
                             fill_value=np.nan, dtype=np.complex128)
    amplitudes_del = np.full(
        shape=(n_energies, n_vib, n_geo, n_beams), fill_value=np.nan, dtype=np.complex128
    )

    # maybe working arrays for transfer into amplitude arrays ?

    # End of the Header - Start of Reading in the Delta Data
    for e_index, line in enumerate(file_lines):  # Energy loop

        # Energy, VPI and VV header
        e_kin, v_imag, v_real = (v for v in ff_reader_6E13_7.read(line) if v is not None)
        # Do NOT translate energy to eV!
        if e_index < n_energies:
            e_kin_array[e_index] = e_kin
            v_inner_array[e_index] = v_real + 1j*v_imag

        # Reference amplitudes
        as_real = read_block(reader=ff_reader_6E13_7, lines=file_lines, shape=(n_beams, 2))
        if e_index < n_energies:
            amplitudes_ref[e_index, :] = as_real.view(dtype=np.complex128)[..., 0]

        # Delta amplitudes
        as_real = read_block(reader=ff_reader_6E13_7, lines=file_lines, shape=(n_geo_vib_grid*n_beams, 2))
        as_complex = as_real.view(dtype=np.complex128)
        if e_index < n_energies:
            amplitudes_del[e_index, ...] = as_complex.reshape(n_vib, n_geo, n_beams)

    if e_index > n_energies:
        raise ValueError("Number of energies does not match number of blocks in file: "
                        f"Found {e_index} blocks")

    return (
        (phi, theta),
        trar,
        (n_beams, n_geo, n_vib), #available from the shapes of matrices
        beam_indices,
        geo_delta,
        vib_delta,
        e_kin_array,
        v_inner_array,
        amplitudes_ref,
        amplitudes_del,
    )


def read_block(reader, lines, shape, dtype=np.float64):
    """ This function reads in the individual blocks of a 
    TensErLEED delta-amplitude file.
    
    Parameters
    ----------
    reader: FortranReader object
        The Fortran reader that is used on this block.
    lines: iterator
        Lines of the whole data file.
    shape: numpy.ndarray
        Shape of the array that gets filled with information.
    dtype: np.float64
        Type of numbers that can be stored in the returned array.
        
    Returns
    -------
    np.array(llist,dtype): numpy.ndarray
        Array with the contents of an individual block with the
        dimensions defined in "shape".
    """
    llist = []
    len_lim = np.prod(shape)
    for line in lines:
        llist.extend((v for v in reader.read(line) if v is not None))
        if len(llist) >= len_lim:
            break
    return np.array(llist, dtype=dtype).reshape(shape)


#@njit(fastmath=True, parallel=True, nogil=True)
def calc_delta_intensities(
    phi,  
    theta, # TODO: phi and theta can be grouped in a beam_incidence tuple
    trar1,
    trar2, # TODO: trar1 and trar2 can be grouped into one
    n_beams, # TODO: not needed, can take from the shape of beam_indices
    beam_indices,
    geo_delta, # formerly CDISP in TensErLEED
    e_kin_array,
    v_inner_array,
    amplitudes_ref,
    amplitudes_del,
    n_files, # can take from shape of geo_delta/n_geo/is_surface_atom
    n_geo,   # not needed, can take from delta_steps
    is_surface_atom,
    delta_steps,
):
    """This function reads in the values of the function Transform and 
    uses them to get the intensity_matrix

    Parameters
    ----------
    phi : float
        Polar incidence angle.
    theta : float
        Azimuthal incidence angle.
    trar1 : numpy.ndarray
        Reciprocal lattice vectors; shape=(2, ).
    trar2 : numpy.ndarray
        Reciprocal lattice vectors; shape=(2, ).
    n_beams : int
        Number of beams for which delta-amplitudes were present
        in the file   
    beam_indices : numpy.ndarray
        Array of beam indices, with beam_indices[i] == [h_i, k_i];
        shape=(n_beams,2) 
    geo_delta: numpy.ndarray
        Geometric displacements of the atoms
        TODO: formerly CDisp; can be taken from delta_steps
    e_kin_array : numpy.ndarray
        Array that contains the kinetic energies of the
        elastically scattered electrons inside the crystal.
        shape=(n_energies)
        (Not the incidence energy!) 
    v_inner_array: numpy.ndarray
        Real and imaginary part of the inner potential
        shape=(n_energies); filled with complex values
    amplitudes_ref: numpy.ndarray
        Array that contains all values of the reference amplitudes 
        shape=(n_energies,n_beams)
    amplitudes_del: numpy.ndarray
        Array that contains all values of the delta amplitudes
        for all atoms/files
        shape=(n_files,n_energies,n_vib,n_geo,n_beams)
    n_files: int
        Number of files for this surface
    n_geo: numpy.ndarray
        Number of geometric displacements
        shape=(n_files)
        TODO: can be taken from "delta_steps"
    is_surface_atom: numpy.ndarray
        Array filled with Bool that indicates if the atom is
        on the surface
        shape=(n_files)
    delta_steps: numpy.ndarray
        Filled with the index of the displacement.
        If given integer values, the values from the Delta
        file will be used, for float a linear interpolation
        is used instead. For file i element (i, 0) gives the
        geometric displacment and element (i, 1) gives the
        vibrational displacement. If the array was falttened,
        geometric displacement is element (2*i) and
        vibrational displacement is element (2*i + 1).
        shape=(n_files, 2)

    Returns
    ----------
    intensity_matrix: numpy.ndarray
        Array that contains the intensities of the different beams with 
        each energy; intensity_matrix[e_index,beam_index]
        shape=(n_energies,n_beams)
    """
    # may be necessary if array was flattened
    delta_steps = delta_steps.reshape(n_files, 2)

    sum_z_disp = 0
    conc = 1    # conc will be taken as input parameter as well

    n_geo_max = geo_delta.shape[1] # I think?

    # below could also be optimized with array operations
    # if uppermost atom in crystal gets displaced up,
    # inner potential also starts higher up (higher z)
    for atom_index, at_surface in enumerate(is_surface_atom):
        if not at_surface:
            continue

        # out of surface component of geometric displacment
        # not less efficient because numpy knows which segment to pick
        z_values = geo_delta[atom_index, :n_geo[atom_index], 0]
        z_grid = np.arange(n_geo[atom_index]) 
        z_fraction = delta_steps[atom_index, 0]

        # use numpy interpolation - this is probably the way to go
        geo_interp = np.interp(z_fraction, z_grid, z_values)

        sum_z_disp += conc * geo_interp

    CXDisp = min(sum_z_disp, 1000)

    intensity_matrix = np.zeros((len(e_kin_array), n_beams))

    trar = np.empty(shape=(2,2), dtype="float")
    trar[0,:] = trar1
    trar[1,:] = trar2

    e_kin = e_kin_array
    n_en = len(e_kin)
    v_real = v_inner_array.real
    v_imag = v_inner_array.imag
    
    j = complex(0,1)

    # incident wave vector
    in_k =np.sqrt(np.maximum(0, 2*(e_kin-v_real)))
    in_k_par = in_k*np.sin(theta) # parallel component
    bk_2 = in_k_par*np.cos(phi) # shape=(n_energy)
    bk_3 = in_k_par*np.sin(phi) # shape=(n_energy)
    bk_z = np.empty_like(e_kin, dtype="complex64")
    bk_z = 2*e_kin - bk_2**2 - bk_3**2 - 2*j*v_imag
    bk_z = np.sqrt(bk_z)

    # outgoing wave vector components
    bk_components = np.stack((bk_2,bk_3))                               # shape=(n_en, 2)
    bk_components = np.outer(bk_components, np.ones(shape=(n_beams,))
                             ).reshape((n_en,2,n_beams))                # shape=(n_en, 2, n_beams)
    out_wave_vec = np.dot(beam_indices, trar)                           # shape=(n_beams, 2)
    out_wave_vec = np.outer(np.ones_like(e_kin), out_wave_vec
                            ).reshape((n_en,2,n_beams))                 # shape=(n_en, n_beams)
    out_components = bk_components + out_wave_vec

    # out k vector
    out_k = (2*np.outer(e_kin, np.ones(shape=(n_beams,))) # 2*E
            + bk_components[:,0,:]**2                     # + h**2
            + bk_components[:,1,:]**2                     # + k**2
            ).astype(dtype="complex64")
    # z component
    out_k_z = np.empty_like(out_k, dtype="complex64")                   # shape=(n_en, n_beams)
    out_k_z = np.sqrt(out_k - 2*j*np.outer(v_imag, np.ones(shape=(n_beams,))))
    # perpendicular component
    out_k_perp = out_k - 2*np.outer(v_real, np.ones(shape=(n_beams,)))  # shape=(n_en, n_beams)

    # interpolation bounds
    _left = np.floor(delta_steps).astype("int32")
    _right = np.ceil(delta_steps).astype("int32")
    geo_id_left, vib_id_left = _left.T
    geo_id_right, vib_id_right = _right.T

    # prefactors (refaction) from amplitudes to intensities
    a = np.sqrt(out_k_perp)
    c = in_k*np.cos(theta)
    prefactor = abs(np.exp(
        -1j*CXDisp*(
        np.outer(bk_z, np.ones(shape=(n_beams,))) + out_k_z
        )))**2 * a/np.outer(c, np.ones(shape=(n_beams,))).real          # shape=(n_en, n_beams)

    # AMI, TODO:
    # What is left in the loops below could be vectorized.
    # However, we can leave them as is for now, since they will be optimized by numba.

    delta_amplitudes = np.empty(shape=(n_files, ), dtype="complex64")
    # Loop over energies
    for e_index in prange(len(e_kin_array)):
        # Loop over beams
        for beam_index in range(n_beams):
            if out_k_perp[e_index, beam_index] <= 0:
                intensity = 0
            else:
                # get reference amplitudes
                amplitude = amplitudes_ref[e_index, beam_index]
                # loop over delta amplitudes
                delta_amplitudes[...] = 0
                for file in prange(n_files):
                    # interpolation of geometric and vibrational displacements
                    del_l_l = amplitudes_del[file, e_index,
                                            vib_id_left[file], geo_id_left[file], beam_index]
                    del_l_r = amplitudes_del[file, e_index,
                                            vib_id_left[file], geo_id_right[file], beam_index]
                    del_r_l = amplitudes_del[file, e_index,
                                            vib_id_right[file], geo_id_left[file], beam_index]
                    del_r_r = amplitudes_del[file, e_index,
                                            vib_id_right[file], geo_id_right[file], beam_index]
                    # values to interpolate to
                    vib_fraction = delta_steps[file, 1]
                    geo_fraction = delta_steps[file, 0]

                    delta_amplitudes[file] = bilinear_interpolation_np(
                        xy = (vib_fraction, geo_fraction),
                        x1x2=(vib_id_left[file], vib_id_right[file]),
                        y1y2=(geo_id_left[file], geo_id_right[file]),
                        f11f12f21f22=(del_l_l, del_l_r, del_r_l, del_r_r)
                    )

                amplitude += np.sum(delta_amplitudes)
                # AMI, TODO:
                # We should remove this saturation check here, since it is very slow (branch in loop).
                # If needed, rather do it outside as array operation - though I don't see the need for it anyhow...
                amp_abs = abs(amplitude).real
                if amp_abs > 10e10:
                    intensity = 10e20
                else:
                    intensity = amp_abs**2 * prefactor[e_index, beam_index]
                intensity_matrix[e_index, beam_index] = intensity.real

    return intensity_matrix


#@njit(fastmath = True)
def bilinear_interpolation_np(xy, x1x2, y1y2, f11f12f21f22):
    """Bilinear interpolation based on numpy functions.

    Parameter
    ---------
    xy: (float, float)
        floats for the fraction in x and y
        (fractions in the two dimensions)
    x1x2: (int, int)
        x1 and x2, floor and ceiling of x
        floor and ceiling of the first dimension
    y1y2: (int, int)
        y1 and y1, floor and ceiling of y
        floor and ceiling of the second dimension
    f11f12f21f22: (float, float, float, float)
        Values of the four grid points (x1,x2,y1,y2)
        the values are indexed as f_xi_yi

    Returns
    -------
    f_x_y: float
        Interpolated value between the four grid points.
    """
    x, y = xy
    x1, x2 = x1x2
    y1, y2 = y1y2
    f_x1_y1, f_x1_y2, f_x2_y1, f_x2_y2 = f11f12f21f22

    # linear interpolation in 1st coordinate
    f_x_y1 = np.interp(x, (x1, x2), (f_x1_y1, f_x2_y1))
    f_x_y2 = np.interp(x, (x1, x2), (f_x1_y2, f_x2_y2))

    # linear interpolation in 2nd coordinate
    f_x_y = np.interp(y, (y1, y2), (f_x_y1, f_x_y2))
    return f_x_y


class AtomDeltas:
    
    def __init__(self, n_energies) -> None:
        self.n_energies = n_energies
        self.delta_files = []

    
    def add_file(self, file):
        new_file = DeltaFile(self, file, n_energies=self.n_energies)
        self.delta_files.append(new_file)

class DeltaFile:
    # one Delta file
    def __init__(self, delta_atom, path, n_energies) -> None:
        self.delta_atom = delta_atom
        self.file_path = path
        self.filename = None # get filename from path
        
        self.read(self.file_path, n_energies)
        
        if self.n_vib == 1:
            self.single_vib = True
            self.min_vib = self.max_vib = None
        else:
            self.single_vib = False
            self.min_vib = np.min(self.vib_delta)
            self.max_vib = np.max(self.vib_delta)
        if self.n_geo == 1:
            self.single_geo = True
            self.min_geo = self.max_geo = self.geo_step = None
        else:
            self.min_geo = np.min(self.geo_delta)
            self.max_geo = np.max(self.geo_delta)
            # individual delta files are regularly spaced...
            self.geo_step = self.geo_delta[1] - self.geo_delta[0]


    def read(self, file, n_energies):
        """Reads the file associated with DeltaFile object.

        Args:
            n_energies (int): number of energies to read
        """
        (
            (self.phi, self.theta),
            self.trar,
            (self.n_beams, self.n_geo, self.n_vib),
            self.beam_indices,
            self.geo_delta,
            self.vib_delta,
            self.e_kin_array,
            self.v_inner_array,
            self.amplitudes_ref,
            self.amplitudes_del,
        ) = read_delta_file(file, n_energies=n_energies, read_header_only=False)
    
