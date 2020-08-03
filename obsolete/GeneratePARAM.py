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
    - PARAM: parameter input file for refcalc
    - GeneratePARAM.log: Logfile of this script
"""

#import sys
import time
from timeit import default_timer as timer
import logging
import TensErLeedModules as tl
#import fortranformat as ff
import numpy as np


###############################################
#                  MAIN                       #
###############################################

def main():
    #start logger, write to file:
    logname = 'GeneratePARAM.log'
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
        beamlist,beamblocks,beamN=tl.writeAUXBEAMS(write=False)
    except:
        logging.error("Exception while reading _BEAMLIST/IVBEAMS", exc_info=True)
        
    # define Clebsh-Gordon coefficient tables:
    mnlmo = [1, 70, 264, 759, 1820, 3836, 7344, 13053, 21868, 34914, 53560, 79443, 114492, 160952, 221408] 
    mnlm = [1, 76, 284, 809, 1925, 4032, 7680, 13593, 22693, 36124, 55276, 81809, 117677, 165152, 226848]
    
    # start generating output
    output = 'C  Dimension statements for Tensor LEED reference calculation, \nC  version v1.2\n\n'
    output += 'C  1. lattice symmetry\n\n'
    m = rp.SUPERLATTICE.copy()
    if m[0,1] != 0:
        f = tl.lcm(abs(int(m[0,1])),abs(int(m[1,1])))
        m[0] *= f/m[0,1]
        m[1] *= f/m[1,1]
        m[0] -= m[1]*np.sign(m[0,1])*np.sign(m[1,1])
    nl1 = int(round(m[0,0]))
    nl2 = int(round(abs(np.linalg.det(rp.SUPERLATTICE))/nl1))
    ideg = 1
    if sl.planegroup in ['p2','pmm','pmg','pgg','cmm','rcmm']:
        ideg = 2
    elif sl.planegroup in ['p3','p3m1','p31m']:
        ideg = 3
    elif sl.planegroup in ['p4','p4m','p4g']:
        ideg = 4
    elif sl.planegroup in ['p6','p6m']:
        ideg = 6
    output += '      PARAMETER (MIDEG='+str(ideg)+',MNL1='+str(nl1)+',MNL2='+str(nl2)+')\n'
    output += '      PARAMETER (MNL = MNL1*MNL2)\n'
    output += '\nC  2. General calculational quantities\n\n'
    output += '      PARAMETER (MKNBS = '+str(beamblocks)+')\n'
    output += '      PARAMETER (MKNT =  '+str(beamN)+')\n'
    output += '      PARAMETER (MNPUN = '+str(len(beamlist))+', MNT0 = '+str(len(beamlist))+')\n'
    output += '      PARAMETER (MNPSI = '+str(len(rp.phaseshifts))+', MNEL = '+str(len(rp.phaseshifts[0][1]))+')\n'
    output += '      PARAMETER (MLMAX = '+str(rp.LMAX)+')\n'
    output += '      PARAMETER (MNLMO = '+str(mnlmo[rp.LMAX-1])+', MNLM = '+str(mnlm[rp.LMAX-1])+')\n'
    output += '\nC  3. Parameters for (3D) geometry within (2D) unit mesh\n\n'
    output += '      PARAMETER (MNSITE  = '+str(len(sl.sitelist))+')\n'
    output += '      PARAMETER (MNLTYPE = '+str(len(sl.layers))+')\n'
    mnbrav = 0
    mnsub = 0
    mnstack = 0
    for layer in sl.layers:
        if len(layer.atlist) == 1:
            mnbrav += 1
        if len(layer.atlist) > mnsub:
            mnsub = len(layer.atlist)
        if not layer.isBulk:
            mnstack += 1
    output += '      PARAMETER (MNBRAV  = '+str(mnbrav)+')\n'
    output += '      PARAMETER (MNSUB   = '+str(mnsub)+')\n'
    output += '      PARAMETER (MNSTACK = '+str(mnstack)+')\n'
    output += '\nC  4. some derived quantities that must be treated explicitly (dummy dimensions for\n'
    output += 'C     special cases necessary\n\n'
    output += '      PARAMETER (MLMAX1=MLMAX+1)\n'
    output += '      PARAMETER (MLMMAX = MLMAX1*MLMAX1)\n\n'
    output += '      PARAMETER (MNBRAV2 = '+('MNBRAV' if mnbrav > 0 else '1')+')\n\n'
    output += '      PARAMETER (MNCOMP= '+('MNLTYPE-MNBRAV' if mnsub > 1 else '1 ')+')\n'
    output += '      PARAMETER (MLMT  = '+('MNSUB*MLMMAX' if mnsub > 1 else '1 ')+')\n'
    output += '      PARAMETER (MNSUB2= '+('MNSUB * (MNSUB-1)/2' if mnsub > 1 else '1 ')+')\n'
    output += '      PARAMETER (MLMG  = '+('MNSUB2*MLMMAX*2' if mnsub > 1 else '1 ')+')\n'
    output += '      PARAMETER (MLMN  = '+('MNSUB * MLMMAX' if mnsub > 1 else '1 ')+')\n'
    output += '      PARAMETER (MLM2N = '+('2*MLMN' if mnsub > 1 else '1 ')+')\n'
    output += '      PARAMETER (MLMNI = '+('MNSUB*MLMMAX' if mnsub > 1 else '1 ')+')\n'
        
    try:
        wf = open('PARAM', 'w')
        wf.write(output)
        wf.close()
    except:
        logging.error("Failed to write PARAM:", exc_info=True)
    logging.info("Wrote to PARAM successfully.")
   
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