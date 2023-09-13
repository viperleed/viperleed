#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Attach Bulk

Takes a slab POSCAR and adds a bulk POSCAR on the bottom, rescaling the unit
cell. Very primitive script, should be updated to include more recent
functionality.
"""
import time
import logging
import copy

import numpy as np
from timeit import default_timer as timer

from viperleed.calc.files import poscar

__authors__ = ["Florian Kraushofer (@fkraushofer)",
               "Alexander M. Imre (@amimre)"]
__created__ = "2019-06-14"


###############################################
#                  MAIN                       #
###############################################

def add_cli_parser_arguments(parser):
    pass


def cleanup(consoleHandler):
    logging.getLogger().removeHandler(consoleHandler)
    logging.shutdown()


def _read_poscar_from_user_input(name):
    """Return a slab read from a file specified by the user."""
    while True:
        filename = input(f"Enter {name} POSCAR name (Ctrl+C to abort): ")
        if not filename:
            print("Input failed. Please try again.")
            continue
        try:
            slab = poscar.read(filename)
        except FileNotFoundError:
            logging.error(f"{filename} not found.")
        except Exception:
            logging.error("Exception while reading POSCAR", exc_info=True)
            return None
        else:
            logging.info(f'{name.capitalize()} POSCAR was read successfully.')
            return slab


logFormatter = logging.Formatter('%(levelname)s - %(message)s')
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)


def main(args=None):
    # start logger, write to file:
    logname = 'Combine-POSCAR.log'
    logging.basicConfig(level=logging.DEBUG, filename=logname, filemode='w',
                        format='%(levelname)s - %(message)s')
    # add console output to logger:
    logging.getLogger().addHandler(consoleHandler)
    logging.info(f"Starting new log: {logname}\nTime of execution (UTC): "
                 + time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())+"\n")
    starttime = timer()

    # open input files:
    try:
        slab = _read_poscar_from_user_input('slab')
    except KeyboardInterrupt:
        return -1
    if not slab:
        return 1  # Error reading

    try:
        bulk = _read_poscar_from_user_input('bulk')
    except KeyboardInterrupt:
        return -1
    if not slab:
        return 1  # Error reading

    # if a and b vectors differ, cancel operation
    if (np.linalg.norm(slab.ucell[:, 0] - bulk.ucell[:, 0]) > 1e-4
            or np.linalg.norm(slab.ucell[:, 1] - bulk.ucell[:, 1]) > 1e-4):
        logging.error("Slab and bulk unit cell vectors are not equal in a and "
                      "b (error > 1e-4). Stopping operation...")
        cleanup(consoleHandler)
        return(1)
    else:
        logging.debug("Slab and bulk unit cell vectors are equal in a and b "
                      "(error < 1e-4).")

    # if element lists in slab and bulk differ, cancel operation
    if not slab.elements == bulk.elements:
        logging.debug("Slab and bulk elements are not equal. Adding missing "
                      "elements.")
        for el in bulk.elements:
            if el not in slab.elements:
                slab.elements.append(el)
                slab.n_per_elem[el] = 0
    else:
        logging.debug("Slab and bulk elements are identical.")

    batoms = len(bulk.atlist)
    cfact = slab.ucell[2][2]/bulk.ucell[2][2]
    slab.ucell[:, 2] *= (cfact + 1) / cfact        # resize the slab unit cell
    # recalculate c for the slab atoms (undistort & shift)
    for atom in slab.atlist:
        atom.pos[2] = (atom.pos[2]*(cfact/(cfact+1)))+(1/(cfact+1))
        atom.oriN += batoms
    # copy atoms from bulk and add them to the slab
    for atom in bulk.atlist:
        newat = copy.copy(atom)
        # recalculate bulk atom c in the new slab unit cell
        newat.pos[2] /= cfact+1
        slab.atlist.append(newat)
        newat.slab = slab
    for el in slab.elements:
        if el in bulk.elements:
            slab.n_per_elem[el] += bulk.n_per_elem[el]

    slab.sort_by_element()

    try:
        poscar.write(slab, reorder=False)
    except Exception:
        logging.error("Exception while writing CONTCAR:", exc_info=True)
#    print(cfact)

    endtime = timer()
    et = endtime-starttime
    if et >= 3600:
        elapsedTimeStr = (str(int(et/3600))
                          + ":" + (str(int(et/60) % 60)).zfill(2) + " hours.")
    elif et >= 60:
        elapsedTimeStr = (str(int(et/60))
                          + ":" + (str(int(et) % 60)).zfill(2) + " minutes.")
    else:
        elapsedTimeStr = str(et) + " seconds."
    logging.info("Finishing "+logname+" at "
                 + time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
                 + ". \nElapsed time: "+elapsedTimeStr+"\n")
    cleanup(consoleHandler)


if __name__ == "__main__":
    main()
