# -*- coding: utf-8 -*-

"""
Created on Jul 06 2022

@author: Alexander M. Imre

New search
"""

import os
import sys
import logging
import shutil

import numpy as np

import copy

import scipy

import viperleed.tleedmlib.files.iosearch as io
import viperleed.tleedmlib.files.iorfactor as iorfactor
import viperleed.tleedmlib as tl
from viperleed.tleedmlib.files.displacements import readDISPLACEMENTS_block
from viperleed.tleedmlib.files.delta_intensities import Transform, calc_delta_intensities
from viperleed.tleedmlib.wrapped.rfactor import r_factor_new as rf
from viperleed.tleedmlib.wrapped.error_codes import error_codes, check_ierr
from viperleed.tleedmlib.sections.rfactor import get_general_rfactor_params, perform_full_intpol_deriv, get_rfactor_energy_grid

import numba
from numba import njit
from pathlib import Path


logger = logging.getLogger("tleedm.search")


def new_search(sl, rp):

    rp.searchResultConfig = None
    if rp.domainParams:
        initToDo = [(dp.rp, dp.sl, dp.workdir) for dp in rp.domainParams]
    else:
        initToDo = [(rp, sl, ".")]

    for (rpt, slt, path) in initToDo:
        # read DISPLACEMENTS block
        if not rpt.disp_block_read:
            readDISPLACEMENTS_block(rpt, slt,
                                    rpt.disp_blocks[rpt.search_index])
            rpt.disp_block_read = True
        # get Deltas
        if 2 not in rpt.runHistory:
            if "Tensors" in rpt.manifest:
                logger.error("New tensors were calculated, but no new delta "
                             "files were generated. Cannot execute search.")
                raise RuntimeError("Delta calculations was not run for "
                                   "current tensors.")



    print("In new search")

    # TODO: read in the Delta files
    # just: rp.generateSearchPars(sl) ?


    # get relation between independent dependent search pars
    rp.generateSearchPars(sl) # TODO: needs inflated DEL_* files in work
    (real_params, real_params_bounds, geo_idp_params, geo_dep_params, geo_all_params, params_relation) = create_search_params_array(
        rp)

    # get flags for surface atoms
    flag_surface_atom = get_surface_atom_flags(rp, sl, geo_all_params)
    
    # number of energy steps
    E_start, E_stop, E_step = rp.THEO_ENERGIES
    n_E = int((E_stop-E_start)/E_step)
    working_energy_grid = None
    
    exp_grid, exp_id_start, exp_n_E_beams, exp_intensities_in = iorfactor.beamlist_to_array(
        rp.expbeams
    )

    delta_dir = Path(path)/"temp"/"Deltas_unpacked"
    
    # TODO dont do it twice...
    tl.leedbase.getDeltas(rpt.TENSOR_INDEX, basedir=path,
                        targetdir=delta_dir, required=True)
    tl.leedbase.getDeltas(rpt.TENSOR_INDEX, basedir=path,
                          targetdir=path, required=True)
    
    
    # read in using Tobi's function
    delta_data = Transform(n_E, str(delta_dir)+"/")
        
    for_error = False  # flag that sets if the R-factor is used for error calculation
    
    # general parameters for R-factors
    deg, v0i, n_beams, which_r, n_derivs, real_iv_shift = get_general_rfactor_params(rp, for_error)

    # energy_grid of Deltas
    theo_grid = np.arange(start=rp.THEO_ENERGIES[0], stop=rp.THEO_ENERGIES[1], step=rp.THEO_ENERGIES[2])
    
    # TODO: make into parameter or derive from other parameters
    intpol_step = 0.5
    
    # interpolated grid
    working_grid = get_rfactor_energy_grid(rp, real_iv_shift, intpol_step)
    
    # Deal with issue of beam correspondeces between experiment and theory
    assert len(rp.ivbeams) == len(rp.expbeams) # this does not necessarily hold. What then?
    corr_list = iorfactor.getBeamCorrespondence(sl, rp)
    
    # unused
    """
    inverted_corr_explist = []
    for ii in range(len(rp.ivbeams)):
        if ii in corr_list:
            inverted_corr_explist.append(rp.ivbeams.index(ii))
        else:
            inverted_corr_explist.append(None)
    """

    
    # TODO check if this is correct up to this point!!!

    # interpolate & get derivative for experimental beams
    (
        exp_e_start_beams_out,
        exp_n_e_beams_out,
        exp_intpol_int,
        exp_intpol_deriv_y
    ) = perform_full_intpol_deriv(deg, n_derivs, working_grid, n_beams, exp_grid, exp_id_start, exp_n_E_beams, exp_intensities_in)

    # get experimental y function
    
    rf.pendry_y_beamset(
        intensities=exp_intpol_int,
        deriv_y_func=exp_intpol_deriv_y,
        beams_id_start=exp_e_start_beams_out,
        beams_n_e=exp_n_e_beams_out,
        v0i=v0i
        )
    
    # calculate the reference intensities & Y-function
    n_E_intpol = working_grid.shape[0]
    ref_start_id, ref_n_E, ref_intens, ref_deriv_y, ierr = rf.alloc_beams_arrays(
        n_beams, n_E_intpol
    )
    check_ierr(ierr, logger)
    
    # TODO: fix urgent
    n_files = 5
    
    # create buffer arrays for trial structures
    work_start_id, work_n_E, work_intens, work_deriv_y, ierr = rf.alloc_beams_arrays(
        n_beams, n_E_intpol
    )
    check_ierr(ierr, logger)
    
    ## Now we are ready to start
    
    theo_calc_params = {
        "deg" : deg,
        "n_beams":n_beams,
        "n_derivs" : n_derivs,
        "n_E_intpol": n_E_intpol,
        "v0i" : v0i,
        "which_r" : which_r,
        "intpol_step" : intpol_step,
        "intpol_grid" : working_grid,
        "input_grid_theo" : theo_grid,
        "flag_surface_atom":flag_surface_atom,
        "n_files" : n_files, # fix
        "n_z_steps": 30, # fix
        "exp_start_id": exp_e_start_beams_out,
        "exp_n_E": exp_n_e_beams_out,
        "exp_y" : exp_intpol_deriv_y,
        "correspondence": np.array(corr_list),
        
    }
    
    # Optimization loop ->
    
    # TODO: fix
    offsets = np.zeros_like(geo_all_params)

    # compare intensities now
    rfac = calc_rfac_for_config(
                            real_params, 
                            params_relation = params_relation,
                            offsets=offsets,
                            theo_calc_params=theo_calc_params,
                            delta_data=delta_data, 
                            work_n_E=work_n_E,
                            work_start_id=work_start_id,
                            work_intens= work_intens,
                            work_deriv_y= work_deriv_y)
    logger.info(f"Before search start: R = {rfac}")
    
    args = (params_relation,
            offsets,
            theo_calc_params,
            delta_data,
            work_n_E,
            work_start_id,
            work_intens,
            work_deriv_y)


    # Log seed used for search for debugging purposes
    seed = 0
    logger.debug(f"Random seed used in search: {seed}")
    result = scipy.optimize.differential_evolution(
        func=calc_rfac_for_config,
        args = args,
        x0 = real_params,
        bounds=real_params_bounds,
        seed = seed,
        init= 'random',
        atol= 1e-9,
        tol= 1e-8,
        maxiter=50, # fix
        disp=True, 
        mutation=1.3,
        popsize=10,
    )
    
    logger.info(result)
    
    pass

def create_search_params_array(rp):
    geo_idp_params = []
    geo_all_params = []
    geo_dep_params = []

    for parameter in rp.searchpars:
        # for now we only work with geometrical deltas
        if parameter.mode in ['occ', 'dom']:
            continue
        geo_all_params.append(parameter)
        if parameter.linkedTo is None:
            # independent parameter
            geo_idp_params.append(parameter)
        else:
            geo_dep_params.append(parameter)

    #sanity checks
    assert(len(geo_idp_params) == rp.indyPars)
    n_indep_params = len(geo_idp_params)
    n_depen_params = len(geo_dep_params)
    n_total_params = n_indep_params + n_depen_params

    assert(n_total_params == len(rp.search_atlist))  # debug

    # make array of parameters
    real_params = np.zeros([n_indep_params])
    # bounds as required by scipy.optimize
    real_params_bounds = []
    for idp_id, idp_par in enumerate(geo_idp_params):
        # starts at 1  because center is a Fortran index
        real_params[idp_id] = idp_par.center
        # starts at 1 because center is a Fortran index
        real_params_bounds.append((0, idp_par.steps-1))

    params_relation = np.zeros(
        [len(geo_idp_params), n_total_params], dtype="int32")

    #create the correspondence matrix between dependent and independent atoms
    for par_id, param in enumerate(geo_all_params):
        for idp_id, idp_param in enumerate(geo_idp_params):
            if param.linkedTo == idp_param or param == idp_param:
                params_relation[idp_id, par_id] = 1

    # sum must be same as number of parameters in total
    assert(np.sum(params_relation) == len(geo_all_params))
    
    # make 

    return real_params, real_params_bounds, geo_idp_params, geo_dep_params, geo_all_params, params_relation


# equivalent to NC_surf from Fortran code
def get_surface_atom_flags(rp, sl, geo_all_params):
    # flag if atom is considered to be at the surface or not
    is_surface_atom = np.full([len(geo_all_params)], False, dtype="bool")
    surface_atoms = sl.getSurfaceAtoms(rp)
    for id, par in enumerate(geo_all_params):
        if par.atom in surface_atoms:
            is_surface_atom[id] = True
    return is_surface_atom

# needs jit
#@njit(fastmath=True, parallel=True, nogil=True)
def expand_params(real_params, params_relation, offsets):
    """extends the set of real pars into the set of all pars"""
    all_pars = np.dot(real_params, params_relation) + offsets
    return all_pars

# needs jit
#@njit(fastmath=True, parallel=True, nogil=True)
def reconstruct_delta_intensities_from_reduced_parameters(indep_params, params_relation, offsets,
                                                          theo_calc_params, delta_data,
                                                          work_start_id, work_n_E, work_intens, work_deriv_y):
    
    # expand indep params to all params
    all_params = expand_params(indep_params, params_relation, offsets)

    delta_intes = calc_delta_intensities(
        delta_data["phi"],
        delta_data["theta"],
        delta_data["surface_vectors"][0],
        delta_data["surface_vectors"][1],
        delta_data["Beam_variables"],
        delta_data["beam_indices"],
        delta_data["CDisp"],
        delta_data["E_array"],
        delta_data["VPI_array"],
        delta_data["VV_array"],
        delta_data["amplitudes_ref"],
        delta_data["amplitudes_del"],
        theo_calc_params["n_files"],
        theo_calc_params["flag_surface_atom"],
        all_params,
        theo_calc_params["n_z_steps"]
        )
    
    # This is always gonna be full range anyways - can we just optimize this away ?
    delta_id_start = np.array([1]*theo_calc_params["n_beams"], dtype="int32")
    delta_n_E = np.array([len(delta_data["E_array"]),]*theo_calc_params["n_beams"], dtype="int32")
    
    ierr = rf.prepare_beams(
        e_grid_in=theo_calc_params["input_grid_theo"],
        intensities=delta_intes,
        e_start_beams=delta_id_start,  # Fortran indices
        n_e_beams=delta_n_E,
        deg=theo_calc_params["deg"],
        n_derivs=theo_calc_params["n_derivs"],
        e_grid_out=theo_calc_params["intpol_grid"],
        e_start_beams_out=work_start_id,
        n_e_beams_out=work_n_E,
        intpol_intensity=work_intens,
        intpol_derivative=work_deriv_y,
    )
    
    
    return # what do we need to return? just Y or what? nothing cause buffer arrays?


# needs jit

def calc_rfac_for_config(real_params,
                         params_relation,
                         offsets,
                         theo_calc_params, delta_data,
                         work_start_id, work_n_E, work_intens, work_deriv_y):
    # get the interpolated intensities
    reconstruct_delta_intensities_from_reduced_parameters(
        real_params,
        params_relation=params_relation,
        offsets=offsets,
        theo_calc_params=theo_calc_params,
        delta_data=delta_data,
        work_start_id=work_start_id,
        work_n_E=work_n_E,
        work_intens=work_intens,
        work_deriv_y=work_deriv_y
        )
    
    # calc Y-function
    rf.pendry_y_beamset(intensities=work_intens,
        deriv_y_func=work_deriv_y,
        beams_id_start=work_start_id,
        beams_n_e=work_n_E,
        v0i=theo_calc_params["v0i"]
        )
    
    # reshuffle into correct order
    work_deriv_y = work_deriv_y[:, theo_calc_params["correspondence"]]
    
    # calc R-factor vs experimental
    
    rfac, *_ = rf.r_pendry_beamset_y(
        theo_calc_params["intpol_step"],
        theo_calc_params["exp_y"],
        work_deriv_y,
        theo_calc_params["exp_start_id"],
        work_start_id,
        theo_calc_params["exp_n_E"],
        work_n_E,
        4
        )
    
    return rfac


def debug_quick_plot(work_deriv_y, theo_calc_params):
    import matplotlib
    import matplotlib.pyplot as plt
    
    plt.plot(work_deriv_y)
    plt.plot(theo_calc_params["exp_y"])
    plt.savefig("test.png", dpi=300)
    

def test(rp, other_x, otherbeams, exp = None, id =0):
    import matplotlib
    import matplotlib.pyplot as plt
    theo_grid, theo_id_start, theo_n_E_beams, theo_intensities_in = iorfactor.beamlist_to_array(
        rp.theobeams["refcalc"]
    )

    plt.cla()

    if exp is not None:
        plt.plot(theo_grid, theo_intensities_in[:,id], label = "theobeams orig")
        plt.plot(other_x, otherbeams[:,id], label = "reconstructed")
        plt.plot(other_x, exp[:,id]/np.nanmax(exp[:,id])*np.nanmax(theo_intensities_in[:,id]), label = "exp normalized")

    else:
        plt.plot(theo_grid, theo_intensities_in[:,id], label = "theobeams orig")
        plt.plot(other_x, otherbeams[:,id], label = "reconstructed")
    plt.legend()
    plt.savefig("test.png", dpi = 300)
    
def resconstruction_difference(rp,  theo_calc_params, work_deriv_y, work_start_id, work_n_E, v0r_shift=0):
    import matplotlib
    import matplotlib.pyplot as plt

    theo_grid, theo_id_start, theo_n_E_beams, theo_intensities_in = iorfactor.beamlist_to_array(
        rp.theobeams["refcalc"]
    )
    
    # interpolate & get derivative for experimental beams
    (
        theo_e_start_beams_out,
        theo_n_e_beams_out,
        theo_intpol_int,
        theo_intpol_deriv_y
    ) = perform_full_intpol_deriv(theo_calc_params["deg"], theo_calc_params["n_derivs"], theo_calc_params["intpol_grid"], theo_calc_params["n_beams"], theo_grid, theo_id_start, theo_n_E_beams, theo_intensities_in)

    # get experimental y function
    
    rf.pendry_y_beamset(
        intensities=theo_intpol_int,
        deriv_y_func=theo_intpol_deriv_y,
        beams_id_start=theo_e_start_beams_out,
        beams_n_e=theo_n_e_beams_out,
        v0i=theo_calc_params["v0i"]
        )
    
    # shuffle
    work_deriv_y = work_deriv_y[:, theo_calc_params["correspondence"]]
    theo_intpol_deriv_y = theo_intpol_deriv_y[:,
                                              theo_calc_params["correspondence"]]
        
    rfac, *_ = rf.r_pendry_beamset_y(
        theo_calc_params["intpol_step"],
        theo_intpol_deriv_y,
        work_deriv_y,
        theo_e_start_beams_out,
        work_start_id,
        theo_n_e_beams_out,
        work_n_E,
        v0r_shift
        )
    print(f"R fac difference: {rfac}")
    
    rfac_work_exp, *_ = rf.r_pendry_beamset_y(
        theo_calc_params["intpol_step"],
        theo_calc_params["exp_y"],
        work_deriv_y,
        theo_calc_params["exp_start_id"],
        work_start_id,
        theo_calc_params["exp_n_E"],
        work_n_E,
        v0r_shift
    )
    print(f"R fac work exp: {rfac_work_exp}")
    
    rfac_theo_exp, *_ = rf.r_pendry_beamset_y(
        theo_calc_params["intpol_step"],
        theo_calc_params["exp_y"],
        theo_intpol_deriv_y,
        theo_calc_params["exp_start_id"],
        theo_e_start_beams_out,
        theo_calc_params["exp_n_E"],
        theo_n_e_beams_out,
        v0r_shift
    )
    print(f"R fac theo exp: {rfac_theo_exp}")

