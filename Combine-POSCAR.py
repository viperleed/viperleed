# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 13:58:38 2019

@author: Florian Kraushofer

Takes a slab POSCAR and adds a bulk POSCAR on the bottom, rescaling the unit cell.
"""

import time
import logging
import copy
# import tleedmlib as tl
import numpy as np
from timeit import default_timer as timer

from tleedmlib.files.poscar import readPOSCAR, writeCONTCAR

###############################################
#                  MAIN                       #
###############################################

def cleanup(consoleHandler):
    logging.getLogger().removeHandler(consoleHandler)
    logging.shutdown()

def main():
    #start logger, write to file:
    logname = 'Combine-POSCAR.log'
    logging.basicConfig(level=logging.DEBUG, filename=logname, filemode='w', format='%(levelname)s - %(message)s')
    #add console output to logger:
    logFormatter = logging.Formatter('%(levelname)s - %(message)s')
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logging.getLogger().addHandler(consoleHandler)
    logging.info("Starting new log: "+logname+"\nTime of execution (UTC): "+time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())+"\n")
    starttime = timer()

    # open input files:
    filename = ""
    while filename == "":
        filename = input("Enter slab POSCAR name: ")
        if not filename:
            print("Input failed. Please try again.")
        else:
            try:
                slab = readPOSCAR(filename)
                logging.info('Slab POSCAR was read successfully.')
            except FileNotFoundError:
                logging.error(filename+" not found.")
                filename = ""
            except:
                logging.error("Exception while reading POSCAR", exc_info=True)
    
    filename = ""
    while filename == "":
        filename = input("Enter bulk POSCAR name: ")
        if not filename:
            print("Input failed. Please try again.")
        else:
            try:
                bulk = readPOSCAR(filename)
                logging.info('Bulk POSCAR was read successfully.')
            except FileNotFoundError:
                logging.error(filename+" not found.")
                filename = ""
            except:
                logging.error("Exception while reading POSCAR", exc_info=True)
    
    # if a and b vectors differ, cancel operation
#    logging.debug("a_bulk = "+str(bulk.ucell[:,0])+", a_slab = "+str(slab.ucell[:,0]))
#    logging.debug("b_bulk = "+str(bulk.ucell[:,1])+", b_slab = "+str(slab.ucell[:,1]))
    if np.linalg.norm(slab.ucell[:,0]-bulk.ucell[:,0]) > 1e-4 or np.linalg.norm(slab.ucell[:,1]-bulk.ucell[:,1]) > 1e-4:
        logging.error("Slab and bulk unit cell vectors are not equal in a and b (error > 1e-4). Stopping operation...")
        cleanup(consoleHandler)
        return(1)
    else:
        logging.debug("Slab and bulk unit cell vectors are equal in a and b (error < 1e-4).")
        
    # if element lists in slab and bulk differ, cancel operation
    if not slab.elements == bulk.elements:
        logging.debug("Slab and bulk elements are not equal. Adding missing elements.")
        for el in bulk.elements:
            if not el in slab.elements:
                slab.elements.append(el)
                slab.nperelem.append(0)
    else:
        logging.debug("Slab and bulk elements are identical.")
     
    batoms = len(bulk.atlist)
    cfact = slab.ucell[2][2]/bulk.ucell[2][2]
    slab.ucell[:,2] *= (cfact+1)/cfact        #resize the slab unit cell
    for atom in slab.atlist:                #recalculate c for the slab atoms (undistort & shift)
        atom.pos[2] = (atom.pos[2]*(cfact/(cfact+1)))+(1/(cfact+1))
        atom.oriN += batoms
    for atom in bulk.atlist:                #copy atoms from bulk Slab object and add them to the slab Slab (sorry)
        newat = copy.copy(atom)
        newat.pos[2] /= cfact+1             #recalculate bulk atom c in the new slab unit cell
        slab.atlist.append(newat)
        newat.slab = slab
    for i in range(0,len(slab.nperelem)):
        if slab.elements[i] in bulk.elements:
            slab.nperelem[i] += bulk.nperelem[i]
    
    slab.sortByEl()
    
    try:
        writeCONTCAR(slab, reorder=False)
    except:
        logging.error("Exception occured:", exc_info=True)
#    print(cfact)
    
    endtime = timer()
    et = endtime-starttime
    if et >= 3600:
        elapsedTimeStr = str(int(et/3600)) + ":"+(str(int(et/60)%60)).zfill(2) +" hours."
    elif et >= 60:
        elapsedTimeStr = str(int(et/60)) + ":"+(str(int(et)%60)).zfill(2) +" minutes."
    else:
        elapsedTimeStr = str(et) +" seconds."
    logging.info("Finishing "+logname+" at "+time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())+". \nElapsed time: "+elapsedTimeStr+"\n")
    #clean up
    cleanup(consoleHandler)

if __name__ == "__main__":
    main()
