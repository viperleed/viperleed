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
    """This function reads in one file of data and stores the data in arrays, which can be used in later functions (ex.: GetInt)
    
    Parameters
    ----------
    filename : string
    The filename describes which file you want to read in (path from the function location to the file)
    
    n_E : int
    Number of different energy levels. This should be a known factor for your files
    
    
    Returns
    -------
    (phi, theta): tuple of float
    Angles of how the beam hits the sample
    
    (trar1, trar2) : tuple of ndarray
    reciprocal lattice vectors
    
    n_beams : int
    Number of beams that are reflected by the sample 
    
    nc_steps: ndarray
    Number of permutations between direction deltas and vibration deltas
    
    e_kin_array : ndarray
    Array that contains the kinetic energies of the elastically scatted
    electrons inside the crystal. (Not the incidence energy!)
    
    v_imag_array : ndarray
    Imaginary part of the inner potential of the surface
    
    VV_array : ndarray
    Real part of the inner potential of the surface
    
    beam_indices : ndarray
    Array with the order of beams
    
    Cundisp : ndarray
    TBD always 0
    
    CDisp : ndarry
    Geometric displacement of given delta
    
    amplitudes_ref : ndarray
    Array that contains all values of the reference amplitudes
    
    amplitudes_del : ndarray
    Array that contains all values of the delta amplitudes
    """

    # Lists and Arrays needed; E,VPI,VV are numpy arrays that get returned, the others are lists that help saving the other data
    HeaderBlock1 = []
    HeaderBlock2 = []
    e_kin_array = np.full(n_energies, fill_value=np.nan)
    v_inner_array = np.full(n_energies, fill_value=np.nan, dtype=np.complex128)

    # we need two fortran format readers
    ff_reader_6E13_7 = ff.FortranRecordReader("6E13.7")
    ff_reader_10F10_5 = ff.FortranRecordReader("10F10.5")
    ff_reader_10F7_4 = ff.FortranRecordReader("10F7.4")

    # Reading in the data of a file
    try:
        with open(filename, mode="r") as file:
            content = file.readlines()
    except Exception:
        raise RuntimeError(f"Unable to read Delta file: {filename}")
    # make into an iterator
    file_lines = iter(content)

    # 1.Block of Header - only 1 line - theta, phi, trar1, trar2 variables
    line = next(file_lines)
    HeaderBlock1 = ff_reader_6E13_7.read(line)
    
    # surface unit cell vectors - what used to be trar1 is now trar[:,0].
    # Similarly trar2 -> trar[:,1]
    trar = np.full(shape=(2, 2), dtype=np.float64, fill_value=np.nan)
    
    theta, phi, trar[0, 0], trar[0, 1], trar[1, 0], trar[1, 1] = HeaderBlock1


    # 2.Block of Header - also only 1 line - n_beams, n_atoms, n_geo_vib_grid variables
    line = next(file_lines)
    z2_elemente = line.split()
    for part in z2_elemente:
        HeaderBlock2.append(int(part))
    
    if len(HeaderBlock2) == 2:
        n_beams, n_geo_vib_grid = HeaderBlock2
    elif len(HeaderBlock2) == 3:
        n_beams, n_atoms, n_geo_vib_grid = HeaderBlock2
        if (n_atoms != 1):
            raise ValueError("A value NATOMS in Delta file was read as != 1. This is not"
                                "supported. Rerun Refcalc and Delta calulation with a newer"
                                "TensErLEED version.")
    else:
        raise RuntimeError("Error reading Delta file header block."
                           "Check Delta file format.")
    

    # 3.Block of Header - (h,k) indices of the beams
    beam_indices = read_block(reader=ff_reader_10F10_5, lines=file_lines, shape=(n_beams, 2))

    # TODO: if we decide to throw CUNDISP out of TensErLEED entirely, this block needs to become optional
    # 4.Block of Header - position of undisplaced atom (Coordinates UNDISPlaced)
    # Unused quantity - only check if it is zero (as it should be)
    pos_undisplaced = read_block(reader=ff_reader_10F7_4, lines=file_lines, shape=(3,))
    if (np.linalg.norm(pos_undisplaced) > 1e-6):
        raise ValueError("A non-zero value of CUNDISP (undisplaced atom positions) "
                         "was read from Delta file. This quantity is currently unused "
                         "and should always be zero. Rerun Refcalc and Delta calulation "
                         "with a newer TensErLEED version.")

    # 5.Block of Header - geometric displacements (Coordinates DISPlaced)
    # For now, this contains, along the first axis, n_vib repetitions of the same
    # displacements. We will figure out n_vib firther below, then reshape this
    geo_delta = read_block(reader=ff_reader_10F7_4, lines=file_lines, shape=(n_geo_vib_grid, 3))

    # 6.Block of Header - list of (vib, 0,0,0,0,..., vib, 0,0,0,0,...)
    vib_delta = read_block(reader=ff_reader_10F7_4, lines=file_lines, shape=(n_geo_vib_grid,))
    n_vib = sum(abs(v)>1e-4 for v in vib_delta)
    n_geo = n_geo_vib_grid // n_vib
    assert n_geo_vib_grid % n_vib == 0
    geo_delta = geo_delta.reshape(n_vib, n_geo, 3)[0, :, :].reshape(n_geo,3)
    # throw out the zeros from array vib_delta
    vib_delta = vib_delta[::n_geo]
    
    if read_header_only:
        amplitudes_ref = None
        amplitudes_del = None
    else:

        # Initialize arrays for reference and delta amplitudes
        amplitudes_ref = np.full(shape=(n_energies, n_beams), fill_value=np.nan, dtype=np.complex128)
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
        (n_beams, n_geo, n_vib),
        beam_indices,
        geo_delta,
        vib_delta,
        e_kin_array,
        v_inner_array,
        amplitudes_ref,
        amplitudes_del,
    )

def read_block(reader, lines, shape, dtype=np.float64):
    llist = []
    len_lim = np.prod(shape)
    for line in lines:
        llist.extend((v for v in reader.read(line) if v is not None))
        if len(llist) >= len_lim:
            break
    return np.array(llist, dtype=dtype).reshape(shape)


@njit(fastmath=True, parallel=True, nogil=True)
def calc_delta_intensities(
    phi,
    theta,
    trar1,
    trar2,
    n_beams,
    beam_indices,
    geo_delta, # formerly CDISP in TensErLEED
    e_kin_array,
    v_inner_array,
    amplitudes_ref,
    amplitudes_del,
    n_files,
    n_geo,
    is_surface_atom,
    delta_steps,
):
    """This function reads in the values of the function Transform and uses them to get the ATSAS_matrix
    
    Parameters
    ----------
    n_files: int
    Number of files
    
    phi, theta:
    Angles of how the beam hits the sample
    
    trar1, trar2:
    reciprocal unit vectors
    
    n_beams:
    todo
    
    beam_indices:
    Array with the order of beams
    
    ph_CDisp:
    Geometric displacements of the atom
    
    E_kin_array:
    Array that contains all the energies of the file
    
    VPI_array:
    Imaginary part of the inner potential of the surface
    
    VV_array:
    Real part of the inner potential of the surface
    
    amplitudes_ref:
    Array that contains all values of the reference amplitudes
    
    amplitudes_del:
    Array that contains all values of the delta amplitudes
    
    nc_surf: np.array of bool
    Bool array with flags that decide if atom is considered to be at the surface.
    
    delta_steps: np.array of float with shape (n_files, 2) of (n_files*2)
    Which displacements to use. If given integer values, the values from the Delta
    file will be used, for float a linear interpolation is used instead.
    For file i element (i, 0) gives the geometric displacment and element
    (i, 1) gives the vibrational displacement.
    If the array was falttened, geometric displacement is element (2*i)
    and vibrational displacement is element (2*i + 1).
    
    number_z_steps:
    Total number of different delta_z values
    
    
    Returns
    ----------
    ATSAS_matrix:
    Array that contains the intensities of the different beams with each energy
    """
    
    # may be necessary if array was flattened
    delta_steps = delta_steps.reshape(n_files, 2)

    # Conc will probably be taken as Input parameter aswell

    XDisp = 0
    Conc = 1
    number_z_steps = amplitudes_del.shape[2] -1 # not needed any more?
    
    n_geo_max = geo_delta.shape[1] # I think?


    # below could also be optimized with array operations
    for i in range(n_files):
        
        
        if is_surface_atom[i]:
            # out of surface component of geometric displacment
            interp_y_values = geo_delta[i, :n_geo[i], 0]
            interp_x_values = np.arange(n_geo[i]) # there is at least one geo displacement (=0.0)
            z_fraction = delta_steps[i, 0]
            
            # use numpy interpolation - this is probably the way to go
            geo_interp = np.interp(z_fraction, interp_x_values, interp_y_values)

            XDisp += Conc * geo_interp
        
        CXDisp = np.min(np.array((XDisp, 1000)))

    # if no surface_atoms were changed, do not change onset height of damping
    if not np.any(is_surface_atom): # isn't this redundant?
        CXDisp = 0


    intensity_matrix = np.zeros((len(e_kin_array), n_beams))
    
    # many optimizations possible here...
    
    
    # Loop over energies
    for e_index in prange(len(e_kin_array)):
        # Definieren von Variablen, die in der jeweiligen Energie gleichbleiben
        E = e_kin_array[e_index]
        VV = v_inner_array[e_index].real
        VPI = v_inner_array[e_index].imag
        # all of the below can be transformed into a neater form with matrix operations
        AK = sqrt(max(2 * E - 2 * VV, 0))
        C = AK * cos(theta)
        BK2 = AK * sin(theta) * cos(phi)
        BK3 = AK * sin(theta) * sin(phi)
        BKZ = cmath.sqrt(complex(2 * E - BK2 ** 2 - BK3 ** 2, -2 * VPI))
        
        # Loop über die Beams
        for beam_index in range(n_beams):
            # Variablen per Beam
            h = beam_indices[beam_index, 0]
            k = beam_indices[beam_index, 1]
            # could be done in matrix form - not sure if that gives better performance
            AK2 = BK2 + h * trar1[0] + k * trar2[0]
            AK3 = BK3 + h * trar1[1] + k * trar2[1]
            AK = 2 * E - AK2 ** 2 - AK3 ** 2
            AKZ = complex(AK, -2 * VPI) #TODO: this is missing a sqrt()!! (line 2816 in lib.search in TensErLEED)
            A_perpendicular = AK - 2 * VV

            amplitude = amplitudes_ref[e_index, beam_index]
            for i in range(n_files):
                # interpolation of geometric and vibrational displacements
                
                geo_fraction = delta_steps[i, 0]
                geo_id_left = np.int32(np.floor(geo_fraction))
                geo_id_right = np.int32(np.ceil(geo_fraction))
                
                vib_fraction = delta_steps[i, 1]
                vib_id_left = np.int32(np.floor(vib_fraction))
                vib_id_right = np.int32(np.ceil(vib_fraction))

                del_l_l = amplitudes_del[i, e_index,
                                         vib_id_left, geo_id_left, beam_index]
                del_l_r = amplitudes_del[i, e_index,
                                         vib_id_left, geo_id_right, beam_index]
                del_r_l = amplitudes_del[i, e_index,
                                         vib_id_right, geo_id_left, beam_index]
                del_r_r = amplitudes_del[i, e_index,
                                         vib_id_right, geo_id_right, beam_index]

                amp_intpol_value = bilinear_interpolation_np(
                    xy = (vib_fraction, geo_fraction),
                    x1x2=(vib_id_left, vib_id_right),
                    y1y2=(geo_id_left, geo_id_right),
                    f11f12f21f22=(del_l_l, del_l_r, del_r_l, del_r_r)
                )

                amplitude += amp_intpol_value

            if A_perpendicular > 0:
                A = sqrt(A_perpendicular)
                PRE = (BKZ + AKZ) * CXDisp
                PRE = cmath.exp(complex(0, -1) * PRE)
                amp_abs = abs(amplitude)
                if amp_abs > 10e10:
                    intensity = 10e20 #saturation condition - do we need this? - can we instead check somewhere else if intensity is too high?
                else:
                    intensity = amp_abs * amp_abs * abs(PRE) * abs(PRE) * A / C
            else:
                intensity = 0

            intensity_matrix[e_index, beam_index] = intensity

    return intensity_matrix


@njit(fastmath = True)
def bilinear_interpolation_np(xy, x1x2, y1y2, f11f12f21f22):
    """Bilinear interpolation based on numpy functions.

    Args:
        xy (_type_): _description_
        x1x2 (_type_): _description_
        y1y2 (_type_): _description_
        f11f12f21f22 (_type_): _description_

    Returns:
        _type_: _description_
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
    

