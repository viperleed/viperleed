
import numpy as np



# read in delta amplitudes into numpy array


# get intensities on the fly for parameter combination




def read_delta(rp, sl, delta_file):
    # reading of delta files by Tobias

    return
    #

# implement MPI eventually

# what is the grid?
# deal with grid, etc.

# Interpolate the deltas
# either numpy interp or scipy.interpolate. If n-dimensional, only linear interpolation is possible


#

# Optimization loop

#@jit
def rfactor_form_search_pars(real_pars, params_realation, offsets):
    """Callable R-factor function for the structure search which is passed to
     optimization algorithms found in scipy.minimize. Returns the R-factor of a given set
     of experimental and theoretical beams (with deltas)."""



    #expand_params
    #call rfactor
    R = 1

    return R




# what are the independent parameters; how do I get them; how do I get the others from there?

#symref is geometric relation, but what gives the index realation ?


def create_search_params_array(rp):
    #
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

    n_pars = len(geo_idp_params) + len(geo_dep_params)

    assert(n_pars == len(rp.search_atlist))

    # make array of parameters
    real_params = np.zeros([len(geo_idp_params)])
    # bounds as required by scipy.optimize
    real_params_bounds = []
    for idp_id, idp_par in enumerate(geo_idp_params):
        real_params[idp_id] = idp_par.center - 1 # -1 here because center is a Fortran index
        real_params_bounds.append((0, idp_par.steps -1)) # -1 here because center is a Fortran index

    params_relation = np.zeros([len(geo_idp_params), n_pars], dtype="int32")

    #create the correspondence matrix between dependent and independent atoms
    for par_id, param in enumerate(geo_all_params):
        for idp_id, idp_param in enumerate(geo_idp_params):
            if param.linkedTo == idp_param or param == idp_param:
                params_relation[idp_id, par_id] = 1

    # sum must be same as number of parameters in total
    assert(np.sum(params_relation) == len(geo_all_params))

    return real_params, real_params_bounds, params_relation, geo_idp_params, geo_all_params

def get_nc_surf(rp, sl, geo_all_params):
    nc_surf = np.full([len(geo_all_params)], False, dtype="bool") # flag if atom is considered to be at the surface or not
    surface_atoms = sl.getSurfaceAtoms(rp)
    for id, par in enumerate(geo_all_params):
        if par.atom in surface_atoms:
            nc_surf[id] = True

    return nc_surf

#@jit
def expand_params(real_params, params_relation, offsets):
    """extends the set of real pars into the set of all pars"""
    all_pars = np.dot(real_params, params_relation) + offsets

    return all_pars


