# -*- coding: utf-8 -*-
"""
Created on Jul 07 2019

@author: Florian Kraushofer

Needs input files: 
    - POSCAR: atom positions per element
    - PARAMETERS: user-defined parameters
    - _PHASESHIFTS: phaseshifts per L and element for each energy
    - _BEAMLIST: original list of all beams
    - IVBEAMS: list of beams requested as output
Writes output file:
    - AUXNONSTRUCT: non-structural parameters for FIN
    - GenerateAUXNONSTRUCT.log: Logfile of this script
"""

#import sys
import time
from timeit import default_timer as timer
import logging
import TensErLeedModules as tl
import fortranformat as ff
#import numpy as np


###############################################
#                  MAIN                       #
###############################################

def main():
    #start logger, write to file:
    logname = 'GenerateAUXNONSTRUCT.log'
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
    
    try:
        beamnums,_,_=tl.writeAUXBEAMS(write=False)
    except:
        logging.error("Exception while reading _BEAMLIST/IVBEAMS", exc_info=True)
    
    # start generating output
    output = ''
    f74 = ff.FortranRecordWriter('F7.4')
    ol = f74.write([rp.INTERLAYER_MAXATTENUATION])
    output += ol+'           >>>>> ! <<<<<              TST\n'
    i4x15 = ff.FortranRecordWriter('15I4')
    ol = i4x15.write(beamnums)
    output += ol+'\n'
    f72f61 = ff.FortranRecordWriter('(F7.2, F6.1)')  #formatting here is a bit crazy... if ever the fortran side is changed, have a look at this
    ol = f72f61.write([rp.THETA, rp.PHI])
    ol = ol.ljust(45)
    output += ol+'THETA FI\n'
    ol = f74.write([rp.BULKDOUBLING_EPS])
    ol = ol.ljust(45)
    output += ol+'EPS\n'
    i3 = ff.FortranRecordWriter('I3')
    ol = i3.write([rp.BULKDOUBLING_ITER])
    ol = ol.ljust(45)
    output += ol+'LITER\n'
    ol = i3.write([rp.LMAX])
    ol = ol.ljust(45)
    output += ol+'LMAX\n'
    
    try:
        wf = open('AUXNONSTRUCT', 'w')
        wf.write(output)
        wf.close()
    except:
        logging.error("Failed to write AUXNONSTRUCT:", exc_info=True)
    logging.info("Wrote to AUXNONSTRUCT successfully.")
   
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