# -*- coding: utf-8 -*-
"""
Created on Oct 07 2019

@author: Florian Kraushofer

For testing different functionalities
"""

#import sys
import time
from timeit import default_timer as timer
import logging
import TensErLeedModules as tl
import os
#import copy
#import fortranformat as ff
#import numpy as np


###############################################
#                  MAIN                       #
###############################################

def main():
    #start logger, write to file:
    logname = 'TestPhaseshiftgen.log'
    logging.basicConfig(level=logging.DEBUG, filename=logname, filemode='w', format='%(levelname)s - %(message)s')
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
    
    try:
        tl.readVIBOCCIN(rp,sl)  #gets VIBOCCIN data into slab
    except:
        logging.error("Exception while reading VIBOCCIN file", exc_info=True)
    
    
    try:
#        dirname = os.path.dirname(os.path.abspath(__file__))
        dirname = '.'
        rundgrenpath = os.path.join(dirname, 'EEASiSSS.x')
        serneliuspath = os.path.join(dirname, 'seSernelius')
        tl.runPhaseshiftGen(sl,rp,psgensource = rundgrenpath, excosource = serneliuspath)
        logging.info("Generated phaseshift file")
    except:
        logging.error("Exception occured while calling phaseshiftgen:", exc_info=True)
    
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