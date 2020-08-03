# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

For testing different functionalities
"""

#import sys
import time
from timeit import default_timer as timer
import logging
import TensErLeedModules as tl
#import copy
#import fortranformat as ff
#import numpy as np


###############################################
#                  MAIN                       #
###############################################

def main():
    #start logger, write to file:
    logname = 'Initialize.log'
    logging.basicConfig(level=logging.DEBUG, filename=logname, filemode='w', format='%(levelname)s - %(message)s')
    #add console output to logger:
    #add console output to logger:
    logFormatter = logging.Formatter('%(levelname)s - %(message)s')
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logging.getLogger().addHandler(consoleHandler)
    logging.info("Starting new log: "+logname+"\nTime of execution (UTC): "+time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())+"\n")
    starttime = timer()
    
    #read all required files
    try:
        sl = tl.readPOSCAR()
    except:
        logging.error("Exception while reading POSCAR", exc_info=True)

    try:
        rp = tl.readPARAMETERS(slab=sl)
    except:
        logging.error("Exception while reading PARAMETERS file", exc_info=True)
    
    sl.fullUpdate(rp)   #gets PARAMETERS data into slab
    
    sl.findSymmetry(rp)
    sl.enforceSymmetry(rp)
    
    sl.sortOriginal()
    
    try:
        tl.writeCONTCAR(sl, filename='AUXPOSCAR', comments='all')
    except:
        logging.error("Exception occured while writing AUXPOSCAR:", exc_info=True)
    
    sl.revertUnitCell()
    
    try:
        tl.writeCONTCAR(sl, filename='AUXPOSCAR_oricell', comments='nodir')
    except:
        logging.error("Exception occured while writing AUXPOSCAR_oricell:", exc_info=True)
   
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
    logging.getLogger().removeHandler(consoleHandler)
    logging.shutdown()

if __name__ == "__main__":
    main()