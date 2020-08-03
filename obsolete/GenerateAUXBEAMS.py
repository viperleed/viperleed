# -*- coding: utf-8 -*-
"""
Created on Jun 21 2019

@author: Florian Kraushofer

Needs input files: 
    - _BEAMLIST: original list of all beams
    - IVBEAMS: list of beams requested as output
Writes output file:
    - AUXBEAMS: Beams specified in IVBEAMS, formatted like in BEAMLIST, ready as input for FIN
    - GenerateAUXBEAMS.log: Logfile of this script
    
"""

#import sys
import time
from timeit import default_timer as timer
import logging
import TensErLeedModules as tl
#import fortranformat as ff


###############################################
#                  MAIN                       #
###############################################

def main():
    #start logger, write to file:
    logname = 'GenerateAUXBEAMS.log'
    logging.basicConfig(level=logging.DEBUG, filename=logname, filemode='w', format='%(levelname)s - %(message)s')
    #add console output to logger:
    logFormatter = logging.Formatter('%(levelname)s - %(message)s')
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logging.getLogger().addHandler(consoleHandler)
    logging.info("Starting new log: "+logname+"\nTime of execution (UTC): "+time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())+"\n")
    starttime = timer()
    
    try:
        tl.writeAUXBEAMS()
    except:
        logging.error("Exception while writing AUXBEAMS", exc_info=True)
    
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