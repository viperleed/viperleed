# -*- coding: utf-8 -*-
"""
Created on June 4 2021

@author: Florian Kraushofer

Testing tleedm execution directly from a database structure.
Requires some manual modification before passing to tleedm because of symmetric
slab and incorrect c-axis.
"""

import ase.db
import sys
import os
import numpy as np
import shutil

# path to directory containing viperleed source - define explicitly
vpr_path = '/home/path/to/source/'       # without final /viperleed

if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))

from viperleed.tleedm import run_tleedm
from viperleed.tleedmlib import Slab, Rparams
from viperleed.tleedmlib.files.poscar import writePOSCAR

def main():
    # download structure as ase.Atoms from ind_label
    ind_label = '4x1_SrTiO3_pristine_Founder'
    db = ase.db.connect(
        'postgresql://hylozoism:Go8J0fnlRuij@lugus.theochem.tuwien.ac.at:5432/clinamen'
    )
    atoms = (next(db.select(ind_label=ind_label))
             .toatoms())

    # create tleedm.Slab from ase.Atoms
    slab = Slab(ase_atoms=atoms)
    rp = Rparams()  # only needed formally, no function here

    # some simple manipulation to make the slab conform to standards:
    # # switch b and c
    # --> needed only for evolutionary runs  and 4x1_SrTiO3_pristine
    #     4x1_SrTiO3_perturbed. Newer slabs are already correct.
    if (ind_label.startswith("4x1_SrTiO3_perturbed")
            or ind_label.startswith("4x1_SrTiO3_pristine")):
        slab.ucell = np.diag(np.diag(slab.ucell)[[0, 2, 1]])
        for at in slab.atlist:
            at.pos[[1, 2]] = at.pos[[2, 1]]
        slab.getCartesianCoordinates(updateOrigin=True)

    # cut below c = 0.4
    slab.atlist = [at for at in slab.atlist if at.pos[2] > 0.4]
    slab.updateAtomNumbers()
    slab.updateElementCount()
    
    # # MR: here one can scale the unit cell to fit the experimental
    # # lattice constants of the bulk.  This is rather easy if only
    # # isotropic scaling is necessary (otherwise one needs to get an
    # # appropriate matrix transformation):
    # a_bulk_experiment = 3.905           # SrTiO3 bulk
    # a_bulk_dft = 16.0542637142856996/4  # b/4 for the 4x1
    # slab.ucell *= a_bulk_experiment / a_bulk_dft
    # slab.getCartesianCoordinates(updateOrigin=True)

    # # test: output POSCAR, just to check
    # writePOSCAR(slab, filename='POSCAR_turned_cut')

    # figure out surface sites
    site_def = {}
    surface_atoms = slab.getSurfaceAtoms(rp)
    for el in slab.elements:
        atn = [at.oriN for at in surface_atoms if at.el == el]
        if atn:
            site_def[el] = {'surf': atn}

    # create work directory if necessary
    # alternative: could use a tempfile.TemporaryDirectory()
    work_path = os.path.abspath(os.path.join('.', 'work_' + ind_label))
    os.makedirs(work_path, exist_ok=True)

    # copy input files to work directory
    # NOTE: PARAMETERS should NOT contain SITE_DEF flags.
    # VIBROCC should contain *_surf sites for all elements.
    for file in ['PARAMETERS', 'VIBROCC', 'IVBEAMS', 'PHASESHIFTS']:
        try:
            shutil.copy2(file, os.path.join(work_path, file))
        except FileNotFoundError:
            pass

    # go to work directory, execute there
    home = os.path.abspath('.')
    os.chdir(work_path)
    run_tleedm(slab=slab, preset_params={'SITE_DEF': site_def},
               source=os.path.join(vpr_path, 'viperleed'))

    # copy back THEOBEAMS.csv
    # alternative: read in as data, store back to database
    if not os.path.isfile('THEOBEAMS.csv'):
        print('ERROR: No THEOBEAMS.csv output found for ' + ind_label)
        return
    theo_name = 'THEOBEAMS_' + ind_label + '.csv'
    try:
        shutil.copy2('THEOBEAMS.csv', os.path.join(home, theo_name))
    except Exception as e:
        print('Error copying ' + theo_name + ' to home directory: ' + str(e))
        return

    # go back, delete work directory (or keep it for diagnosis...)
    os.chdir(home)
    try:
        shutil.rmtree(work_path)
    except Exception as e:
        print('Error deleting work directory: ' + str(e))

    return


if __name__ == '__main__':
    main()
