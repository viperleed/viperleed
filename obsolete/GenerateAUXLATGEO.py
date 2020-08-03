# -*- coding: utf-8 -*-
"""
Created on Jul 07 2019

@author: Florian Kraushofer

Needs input files: 
    - POSCAR: atom positions per element
    - PARAMETERS: user-defined parameters
Writes output file:
    - AUXLATGEO: bulk and surface unit cell, energy range, inner potential onset
    - GenerateAUXLATGEO.log: Logfile of this script
"""

#import sys
import time
from timeit import default_timer as timer
import logging
import TensErLeedModules as tl
import fortranformat as ff
import numpy as np
import os


###############################################
#                  MAIN                       #
###############################################

def main():
    #start logger, write to file:
    logname = 'GenerateAUXLATGEO.log'
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
    
    # start generating output
    output = ''
    dirpath = os.getcwd()
    if os.name == 'nt':     #windows - replace the stupid backslashes...
        dirpath = dirpath.replace('\\','/')
    folders = dirpath.split('/')
    if len(folders) > 1:
        foldername = folders[-2]
    else:
        foldername = 'unknown'
    
    output += foldername+' '+time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())+'\n'       #this originally contained V0i, maybe add again once : V0i=5.10eV db=2.2656A
    f72x3 = ff.FortranRecordWriter('3F7.2')
    ol = f72x3.write(rp.TEOENRANGE)
    ol = ol.ljust(24)
    output += ol + 'EI,EF,DE\n'
    f74x2 = ff.FortranRecordWriter('2F7.4')
    ucsurf = sl.ucell[:2,:2]
    ucbulk = np.dot(np.linalg.inv(rp.SUPERLATTICE),ucsurf)
    ol = f74x2.write(ucbulk[0])
    ol = ol.ljust(24)
    output += ol + 'ARA1\n'
    ol = f74x2.write(ucbulk[1])
    ol = ol.ljust(24)
    output += ol + 'ARA2\n'
    output += ' 0.0    0.0             SS1\n'
    output += ' 0.0    0.0             SS2\n'
    output += ' 0.0    0.0             SS3\n'
    output += ' 0.0    0.0             SS4\n'
    ol = f74x2.write(ucsurf[0])
    ol = ol.ljust(24)
    output += ol + 'ARB1\n'
    ol = f74x2.write(ucsurf[1])
    ol = ol.ljust(24)
    output += ol + 'ARB1\n'
    output += ' 0.0    0.0             SO1\n'
    output += ' 0.0    0.0             SO2\n'
    output += ' 0.0    0.0             SO3\n'
    ol = f74x2.write([0.0, rp.INPOT_Z])
    ol = ol.ljust(24)
    output += ol + 'FR ASE\n'
    
    try:
        wf = open('AUXLATGEO', 'w')
        wf.write(output)
        wf.close()
    except:
        logging.error("Failed to write AUXLATGEO:", exc_info=True)
    logging.info("Wrote to AUXLATGEO successfully.")
   
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