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

# only required if running from within 'viperleed' directory:
cd = os.path.realpath(os.path.dirname(__file__))
vpr_path = os.path.realpath(os.path.join(cd, '..'))
for import_path in (cd, vpr_path):
    if import_path not in sys.path:
        sys.path.append(import_path)

from viperleed.tleedm import run_tleedm
from viperleed.tleedmlib import Slab, Rparams
from viperleed.tleedmlib.files.poscar import writeCONTCAR

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
    # switch b and c
    slab.ucell = np.diag(np.diag(slab.ucell)[[0, 2, 1]])
    for at in slab.atlist:
        at.pos[[1, 2]] = at.pos[[2, 1]]
    slab.getCartesianCoordinates(updateOrigin=True)

    # cut below c = 0.4
    slab.atlist = [at for at in slab.atlist if at.pos[2] > 0.4]
    slab.updateAtomNumbers()
    slab.updateElementCount()

    # test: output POSCAR, just to check
    writeCONTCAR(slab, filename='POSCAR_turned_cut')

    # figure out surface sites
    site_def = {}
    surface_atoms = slab.getSurfaceAtoms(rp)
    for el in slab.elements:
        atn = [at.oriN for at in surface_atoms if at.el == el]
        if atn:
            site_def[el] = {'surf': atn}

    # here one could cd to a temp directory; don't forget to copy input
    #  files there: PARAMETERS, PHASESHIFTS, IVBEAMS, (VIBROCC)

    # run tleedm
    run_tleedm(slab=slab, preset_params={'SITE_DEF': site_def})

    # copy back THEOBEAMS.csv ?
    # ...

    return


if __name__ == '__main__':
    main()
