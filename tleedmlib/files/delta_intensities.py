# -*- coding: utf-8 -*-
"""

@author: Tobias Hable, Alexander M. Imre

Reads in delta files.
"""
from re import I
import sys
import numpy as np
from numpy import sin, cos, sqrt
import fortranformat as ff
import matplotlib.pyplot as plt
import scipy
import os
from tqdm import tqdm
import cmath
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
            AKZ = complex(AK, -2 * VPI)
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

    beam_indices,
    ph_CDisp,
    E_kin_array,
    VPI_array,
    VV_array,
    amplitudes_ref,
    amplitudes_del,
    n_files,
    is_surface_atom,
    delta_steps,
    number_z_steps,
):
    """This function reads in the values of the function Transform and uses them to get the ATSAS_matrix
    
    Parameters
    ----------
    n_files: int
    Number of files
    
    phi, theta:
    Angles of how the beam hits the sample
    
    trar1, trar2:
    Vectors of the normal and the reciprocal unit cell of the sample
    
    Beam_variables:
    The variables int0, n_atoms, nc_steps for each file stored in an array
    
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
    
    delta_steps:
    List of numbers that decide which geometric displacement this atom has
    
    number_z_steps:
    Total number of different delta_z values
    
    
    Returns
    ----------
    ATSAS_matrix:
    Array that contains the intensities of the different beams with each energy
    """

    # Conc will probably taken as Input parameter aswell
    CXDisp = 1000
    XDisp = 0
    Conc = 1
    for i in range(n_files):
        if is_surface_atom[i]:
            n_beams, n_atoms, nc_steps = Beam_variables[i, :]
            for j in range(n_atoms):
                delta_fraction = delta_steps[i]  # probably besserer Name als fill finden
                delta_step_left = int(delta_fraction // 1)
                # necessary so that the x2 doesnt go out of the index range
                if delta_step_left == number_z_steps - 1:
                    delta_step_left -= 1
                delta_step_right = delta_step_left + 1
                x = delta_fraction - delta_step_left
                del_left = ph_CDisp[i, delta_step_left, j, 0]
                del_right = ph_CDisp[i, delta_step_right, j, 0]
                CDisp = x * (del_right - del_left) + del_left
                XDisp = Conc * CDisp
                if XDisp < CXDisp:
                    CXDisp = XDisp

    ATSAS_matrix = np.zeros((len(E_kin_array), n_beams))
    
    # many optimizations possible here...
    
    
    
    # Loop over energies
    for e_index in prange(len(E_kin_array)):
        # Definieren von Variablen, die in der jeweiligen Energie gleichbleiben
        E = E_kin_array[e_index]
        VV = VV_array[e_index]
        VPI = VPI_array[e_index]
        AK = sqrt(max(2 * E - 2 * VV, 0))
        C = AK * cos(theta)
        BK2 = AK * sin(theta) * cos(phi)
        BK3 = AK * sin(theta) * sin(phi)
        BKZ = sqrt(complex(2 * E - BK2 ** 2 - BK3 ** 2, -2 * VPI))

        # Loop über die Beams
        for beam_index in range(n_beams):
            # Variablen per Beam
            h = beam_indices[beam_index, 0]
            k = beam_indices[beam_index, 1]
            # could be done in matrix form - not sure if that gives better performance
            AK2 = BK2 + h * trar1[0] + k * trar2[0]
            AK3 = BK3 + h * trar1[1] + k * trar2[1]
            AK = 2 * E - AK2 ** 2 - AK3 ** 2
            AKZ = complex(AK, -2 * VPI)
            A_perpendicular = AK - 2 * VV

            # Herausfinden welche NCStep deltas man nimmt mit delta_step matrix
            DelAct = amplitudes_ref[e_index, beam_index]
            for i in range(n_files):
                # noch genau schauen wie genau gewollt
                # Interpolation of float delta_step values
                delta_fraction = delta_steps[i]
                delta_step_left = int(delta_fraction)
                delta_step_right = min(delta_step_left+1, number_z_steps)
                del_left = amplitudes_del[i, e_index, delta_step_left, beam_index]
                del_right = amplitudes_del[i, e_index, delta_step_right, beam_index]
                # 1D interpolation formula
                intpol_value = (delta_fraction - delta_step_left) * (del_right - del_left) + del_left
                # ###
                DelAct += intpol_value

            if A_perpendicular > 0:
                A = sqrt(A_perpendicular)
                PRE = (BKZ + AKZ) * CXDisp
                PRE = np.exp(complex(0, -1) * PRE)
                amp_abs = abs(DelAct)
                if amp_abs > 10e10:
                    ATSAS = 10e20
                else:
                    ATSAS = amp_abs ** 2 * abs(PRE) ** 2 * A / C
            else:
                ATSAS = 0

            ATSAS_matrix[e_index, beam_index] = ATSAS

    return ATSAS_matrix


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


@njit(fastmath=True, parallel=True, nogil=True)
def calc_delta_intensities_2D(
    phi,
    theta,
    trar1,
    trar2,
    Beam_variables,
    beam_indices,
    ph_CDisp,
    E_array,
    VPI_array,
    VV_array,
    amplitudes_ref,
    amplitudes_del,
    n_files,
    nc_surf,
    delta_steps,
    number_z_steps,
    number_vib_steps,
):

    """This function reads in the values of the function Transform and uses them to get the ATSAS_matrix
    
    Parameters
    ----------
    phi, theta:
    Angles of how the beam hits the sample
    
    trar1, trar2:
    Vectors of the normal and the reciprocal unit cell of the sample
    
    Beam_variables:
    The variables int0, n_atoms, nc_steps for each file stored in an array
    
    beam_indices:
    Array with the order of beams
    
    ph_CDisp:
    Geometric displacements of the atom
    
    E_array:
    Array that contains all the energies of the file
    
    VPI_array:
    Imaginary part of the inner potential of the surface
    
    VV_array:
    Real part of the inner potential of the surface
    
    amplitudes_ref:
    Array that contains all values of the reference amplitudes
    
    amplitudes_del:
    Array that contains all values of the delta amplitudes
    
    n_files: int
    Number of files
    
    nc_surf: np.array of bool
    Bool array with flags that decide if atom is considered to be at the surface.
    
    delta_step:
    List of numbers that decide which geometric displacement this atom has
    
    number_z_steps:
    Total number of different delta_z values
    
    number_vib_steps:
    Total number of different delta_vib values
    
    Returns
    ----------
    ATSAS_matrix:
    Array that contains the intensities of the different beams with each energy
    """

    # Conc will probably taken in as an Input parameter
    CXDisp = 1000
    XDisp = 0
    Conc = 1
    for i in range(n_files):
        if nc_surf[i]:
            int0, n_atoms, nc_steps = Beam_variables[i, :]
            for j in range(n_atoms):
                z_step = delta_steps[
                    i, 0
                ]
                x1 = int(z_step)
                x2 = min(x1+1, number_z_steps)
                x = z_step - x1
                y1 = ph_CDisp[i, x1, j, 0]
                y2 = ph_CDisp[i, x2, j, 0]
                CDisp = x * (y2 - y1) + y1
                XDisp = Conc * CDisp
                if XDisp < CXDisp:
                    CXDisp = XDisp

    ATSAS_matrix = np.zeros((len(E_array), int0))

    # Loop over the energies
    for e_index in prange(len(E_array)):
        # defining variables that stay the same for their respective energy
        E = E_array[e_index]
        VV = VV_array[e_index]
        VPI = VPI_array[e_index]

        AK = sqrt(max(2 * E - 2 * VV, 0))
        C = AK * cos(theta)
        BK2 = AK * sin(theta) * cos(phi)
        BK3 = AK * sin(theta) * sin(phi)
        BKZ = sqrt(complex(2 * E - BK2 ** 2 - BK3 ** 2, -2 * VPI))

        # Loop over the beams
        for b_index in range(int0):
            # variables that differ with each beam
            h = beam_indices[b_index, 0]
            k = beam_indices[b_index, 1]
            AK2 = BK2 + h * trar1[0] + k * trar2[0]
            AK3 = BK3 + h * trar1[1] + k * trar2[1]
            AK = 2 * E - AK2 ** 2 - AK3 ** 2
            AKZ = complex(AK, -2 * VPI)
            APERP = AK - 2 * VV

            # finding out which NCStep you take with the delta_step matrix
            DelAct = amplitudes_ref[e_index, b_index]
            for i in range(n_files):
                # Interpolation of the float delta_step values

                z = delta_steps[i, 0]
                v = delta_steps[i, 1]
                z1 = int(z)
                z2 = min(z1+1, number_z_steps)
                v1 = int(v) 
                v2 = min(v1+1, number_vib_steps)

                da_z1_v1 = amplitudes_del[i, e_index, v1 * number_vib_steps + z1, b_index]
                da_z1_v2 = amplitudes_del[i, e_index, v2 * number_vib_steps + z1, b_index]
                da_z2_v1 = amplitudes_del[i, e_index, v1 * number_vib_steps + z2, b_index]
                da_z2_v2 = amplitudes_del[i, e_index, v2 * number_vib_steps + z2, b_index]

                intpol_value = bilinear_interpolation_unit_square(
                    (z1, z2, z),
                    (v1, v2, v), 
                    (da_z1_v1, da_z1_v2, da_z2_v1, da_z2_v2)
                )
                DelAct += intpol_value

            if APERP > 0:
                A = sqrt(APERP)
                PRE = (BKZ + AKZ) * CXDisp
                PRE = np.exp(complex(0, -1) * PRE)
                amp_abs = abs(DelAct)
                if amp_abs > 10e10:
                    ATSAS = 10e20
                else:
                    RPRE = abs(PRE)
                    ATSAS = amp_abs ** 2 * RPRE ** 2 * A / C
            else:
                ATSAS = 0

            ATSAS_matrix[e_index, b_index] = ATSAS

    return ATSAS_matrix


def hartree_to_eV(E_hartree):
    return E_hartree * 27.211396


def Transform(n_E, directory):
    """This function transforms the read in data to a form, where it can be read in by the GetInt function

    Parameters
    ----------
    n_E:
    Number of different energies of one file
    
    directory:
    Relative path to the file that gets read in
    
    
    Returns 
    ----------
    phi, theta:
    Angles of how the beam hits the sample
    
    trar1, trar2:
    Vectors of the normal and the reciprocal unit cell of the sample
    
    Beam_variables:
    The variables int0, n_atoms, nc_steps for each file stored in an array
    
    beam_indices:
    Array with the order of beams
    
    CDisp:
    Geometric displacements of the atom
    
    E_array:
    Array that contains all the energies of the file
    
    VPI_array:
    Imaginary part of the inner potential of the surface
    
    VV_array:
    Real part of the inner potential of the surface
    
    amplitudes_ref:
    Array that contains all values of the reference amplitudes
    
    amplitudes_del:
    Array that contains all values of the delta amplitudes
    
    filename_list:
    List of the filenames that contain the data
    """

    filename_list = []
    data_list_all = {}
    for filenames in os.walk(directory):
        filename_list.append(filenames)
    filename_list = filename_list[0][2]

    for name in tqdm(filename_list):
        filename = directory + name
        data_list_all[name] = read_delta_file(filename, n_E)

    # constant variables
    phi, theta = data_list_all[filename_list[0]][0]
    trar1, trar2 = data_list_all[filename_list[0]][1]
    beam_indices = data_list_all[filename_list[0]][3]
    E_array, VPI_array, VV_array = data_list_all[filename_list[0]][6]
    amplitudes_ref = data_list_all[filename_list[0]][7]
    n_files = len(filename_list)

    # int0, n_atoms, nc_steps in an array(nc_steps can change, n_atoms maybe too)
    Beam_variables = np.full(shape=[n_files, 3], fill_value=np.NaN, dtype=int)
    for i, name in enumerate(filename_list):
        Beam_variables[i, :] = data_list_all[name][2]

    int0 = int(np.max(Beam_variables[:, 0]))
    n_atoms_max = int(np.max(Beam_variables[:, 1]))
    nc_steps_max = int(np.max(Beam_variables[:, 2]))

    # saving the changing data in arrays
    CDisp = np.full(shape=(n_files, nc_steps_max, n_atoms_max, 3), fill_value=np.NaN)
    amplitudes_del = np.full(
        shape=(n_files, n_E, nc_steps_max, int0), fill_value=np.NaN, dtype=complex
    )
    for i, name in enumerate(filename_list):
        int0, n_atoms, nc_steps = Beam_variables[i]
        int0 = int(int0)
        n_atoms = int(n_atoms)
        nc_steps = int(nc_steps)
        CDisp[i, 0:nc_steps, 0:n_atoms, :] = data_list_all[name][5]
        amplitudes_del[i, :, 0:nc_steps, :] = data_list_all[name][8]

    return (
        phi,
        theta,
        trar1,
        trar2,
        Beam_variables,
        beam_indices,
        CDisp,
        E_array,
        VPI_array,
        VV_array,
        amplitudes_ref,
        amplitudes_del,
        filename_list,
    )


def PlotMaker(
    phi,
    theta,
    trar1,
    trar2,
    Beam_variables,
    beam_indices,
    CDisp,
    E_array,
    VPI_array,
    VV_array,
    amplitudes_ref,
    amplitudes_del,
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
    
    n_files,
    NCSurf,
    delta_step_einzeln,
    phi_g,
    theta_g,
    trar1_g,
    trar2_g,
    Beam_variables_g,
    beam_indices_g,
    CDisp_g,
    E_array_g,
    VPI_array_g,
    VV_array_g,
    amplitudes_ref_g,
    amplitudes_del_g,
    n_files_g,
    delta_step_g,
    int0,
    number_z_steps,
    number_vib_steps,
):

    # This function determines the difference between a file, where 2 contributors to the intensities are considered, and 2 files, where only one of these
    # contributors each is considered. That way we can determine, if it is viable to use linear combination of the 3 dimensions + the vib files to describe
    # each delta of the atom.
    #
    # INPUT:
    #
    # phi, theta, phi_g, theta_g:
    # Angles of how the beam hits the sample
    #
    # trar1, trar2, trar1_g, trar2_g:
    # Vectors of the normal and the reciprocal unit cell of the sample
    #
    # Beam_variables, Beam_variables_g:
    # The variables int0, n_atoms, nc_steps for each file stored in an array
    #
    # beam_indices, beam_indices_g:
    # Array with the order of beams
    #
    # CDisp, CDisp_g:
    # Geometric displacements of the atom
    #
    # E_array, E_array_g:
    # Array that contains all the energies of the file
    #
    # VPI_array, VPI_array_g:
    # Imaginary part of the inner potential of the surface
    #
    # VV_array, VV_array_g:
    # Real part of the inner potential of the surface
    #
    # amplitudes_ref, amplitudes_ref_g:
    # Array that contains all values of the reference amplitudes
    #
    # amplitudes_del, amplitudes_del_g:
    # Array that contains all values of the delta amplitudes
    #
    # filename_list, filename_list_g:
    # List of the filenames that contain the data
    #
    # NCSurf:
    # List of 0 and 1 to decide which file takes part in creating the CXDisp
    #
    # delta_step, delta_step_g:
    # List of numbers that decide which geometric displacement this atom has
    #
    # int0:
    # number of beams
    #
    # number_z_steps:
    # Total number of different delta_z values
    #
    # number_vib_steps:
    # Total number of different delta_vib values
    #
    # OUTPUT:
    # ATSAS_matrix:
    # Array that contains the intensities of the different beams with each energy

    matrix_differenz = np.full(shape=(len(E_array), int0), fill_value=np.NaN)
    array = np.full(shape=(101, 101), fill_value=np.NaN)
    for i in tqdm(
        range(101)
    ):  # gerade weird dringend drüber nachdenken, 0.1 schritte ist dann diese range sinnvoll?
        for j in range(101):
            delta_step_einzeln = np.array(
                [j / 10, 5.0, 5.0, 5.0, 5.0, i / 10, 5.0, 5.0, 5.0, 5.0]
            )

            delta_step_g = np.full(shape=(5, 2), fill_value=5.0)
            delta_step_g[0, 0] = j / 10
            delta_step_g[0, 1] = i / 10

            matrix_z_vib_einzeln = calc_delta_intensities(
                phi,
                theta,
                trar1,
                trar2,
                Beam_variables,
                beam_indices,
                CDisp,
                E_array,
                VPI_array,
                VV_array,
                amplitudes_ref,
                amplitudes_del,
                n_files,
                NCSurf,
                delta_step_einzeln,
                number_z_steps,
            )

            matrix_z_vib_g = calc_delta_intensities_2D(
                phi_g,
                theta_g,
                trar1_g,
                trar2_g,
                Beam_variables_g,
                beam_indices_g,
                CDisp_g,
                E_array_g,
                VPI_array_g,
                VV_array_g,
                amplitudes_ref_g,
                amplitudes_del_g,
                n_files_g,
                NCSurf,
                delta_step_g,
                number_z_steps,
                number_vib_steps,
            )

            for e_index in range(len(E_array)):
                for b_index in range(int0):
                    matrix_differenz[e_index, b_index] = (
                        matrix_z_vib_g[e_index, b_index]
                        - matrix_z_vib_einzeln[e_index, b_index]
                    )

            Integral = 0
            for b in range(int0):
                Integral = Integral + abs(np.trapz(matrix_differenz[:, b]))

            array[j, i] = Integral

    return array


