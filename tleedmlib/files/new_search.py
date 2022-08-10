
import numpy as np

from viperleed.tleedmlib.files.delta_intensities import calc_delta_intensities


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


class Search():
    
    def __init__(self, run_parameters, slab):
        self.__create_search_params_array(run_parameters)
        
        # Create bool flags for surface atoms
        self.__get_surface_atom_flags(run_parameters, slab)
    
    def __create_search_params_array(self, rp):
        self.geo_idp_params = []
        self.geo_all_params = []
        self.geo_dep_params = []

        for parameter in rp.searchpars:
            # for now we only work with geometrical deltas
            if parameter.mode in ['occ', 'dom']:
                continue
            self.geo_all_params.append(parameter)
            if parameter.linkedTo is None:
                # independent parameter
                self.geo_idp_params.append(parameter)
            else:
                self.geo_dep_params.append(parameter)

        #sanity checks
        assert(len(self.geo_idp_params) == rp.indyPars)
        self.n_indep_params = len(self.geo_idp_params)
        self.n_depen_params = len(self.geo_dep_params)
        self.n_total_params = self.n_indep_params + self.n_depen_params

        assert(self.n_total_params == len(rp.search_atlist)) # debug

        # make array of parameters
        self.real_params = np.zeros([self.n_indep_params])
        # bounds as required by scipy.optimize
        self.real_params_bounds = []
        for idp_id, idp_par in enumerate(self.geo_idp_params):
            self.real_params[idp_id] = idp_par.center # starts at 1  because center is a Fortran index
            self.real_params_bounds.append((0, idp_par.steps)) # starts at 1 because center is a Fortran index

        self.params_relation = np.zeros([len(self.geo_idp_params), self.n_total_params], dtype="int32")

        #create the correspondence matrix between dependent and independent atoms
        for par_id, param in enumerate(self.geo_all_params):
            for idp_id, idp_param in enumerate(self.geo_idp_params):
                if param.linkedTo == idp_param or param == idp_param:
                    self.params_relation[idp_id, par_id] = 1

        # sum must be same as number of parameters in total
        assert(np.sum(self.params_relation) == len(self.geo_all_params))

    
    def __get_surface_atom_flags(self, rp, sl):
        self.is_surface_atom = np.full([len(self.geo_all_params)], False, dtype="bool") # flag if atom is considered to be at the surface or not
        surface_atoms = sl.getSurfaceAtoms(rp)
        for id, par in enumerate(self.geo_all_params):
            if par.atom in surface_atoms:
                self.is_surface_atom[id] = True

    def expand_params(self, real_params):
        """extends the set of real pars into the set of all pars"""
        all_pars = np.dot(real_params, self.params_relation) + self.offsets
        return all_pars
    
    def set_various_params(phi, theta, trar1, trar2):
        pass
    
    def reconstruct_delta_intensities_from_reduced_parameters(self, indep_params):
        all_params = self.expand_params(indep_params)
        
        intpol_intensity = calc_delta_intensities()
        
    def calc_rfac_for_config(self, indep_params):
        intpol_intensity = reconstruct_delta_intensities_from_reduced_parameters(indep_params)
        
        
        
        
        