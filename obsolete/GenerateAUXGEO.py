# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

Needs input files: 
    - POSCAR: atom positions per element
    - PARAMETERS: user-defined parameters
    - VIBOCCIN: user-defined vibrational amplitudes and occupations per site
Writes output file:
    - AUXGEO: list of sites, layers, atoms per layer, and layer stacking, ready as input for FIN
    - GenerateAUXGEO.log: Logfile of this script
"""

#import sys
import time
from timeit import default_timer as timer
import logging
import TensErLeedModules as tl
import fortranformat as ff
import numpy as np


###############################################
#                  MAIN                       #
###############################################

def main():
    #start logger, write to file:
    logname = 'GenerateAUXGEO.log'
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
        tl.readVIBOCCIN(rp,sl)  #gets VIBOCCIN data into slab
    except:
        logging.error("Exception while reading VIBOCCIN file", exc_info=True)
    
    # start generating output
    output = ''
    output += '-------------------------------------------------------------------\n'
    output += '--- define chem. and vib. properties for different atomic sites ---\n'
    output += '-------------------------------------------------------------------\n'
    i3 = ff.FortranRecordWriter('I3')
    ol = i3.write([len(sl.sitelist)])
    ol = ol.ljust(26)
    output += ol + 'NSITE: number of different site types\n'
    f74x2 = ff.FortranRecordWriter('2F7.4')
    for i, site in enumerate(sl.sitelist):
        output += '-   site type '+str(i+1)+' ---\n'
        
        for el in sl.elements:
            # this reproduces the order of blocks contained in _PHASESHIFTS:
            if el in rp.ELSPLIT:
                chemelList = rp.ELSPLIT[el]
            else:
                chemelList = [el]
            siteList = [s for s in sl.sitelist if s.el == el]
            for cel in chemelList:
                for s in siteList:
                    if s.isEquivalent(site):
                        occ, vib = site.occ[cel], site.vibamp[cel]
                        comment = ('Occ & VibAmp for '+cel+' in '+site.label
                                   +' site')
                    else:
                        occ, vib = 0.0, 0.0
                        comment = ''
                    try:
                        ol = f74x2.write([occ, vib])
                    except:
                        logging.error("Exception while trying to write \
                            occupation / vibrational amplitude for site "
                            +site.label, exc_info=True)
                    ol = ol.ljust(26)
                    output += ol + comment + '\n'
        output += (f74x2.write([0.0,0.0]).ljust(26)
                   +'Occ & VibAmp for vacancy\n')
        
#        for el in sl.chemelem:     # TODO: this should go over the phaseshift elements, not the chemelem; in principle, the site names will contain the chemelem names, so the dict can be adressed in the same way
#            try:
#                ol = f74x2.write([site.occ[el], site.vibamp[el]])
#            except:
#                logging.error("Exception while trying to write occupation / vibrational amplitude for site "+site.label, exc_info=True)
#            ol = ol.ljust(26)
#            output += ol+'Occ & VibAmp for '+el+' in '+site.label+' site\n'
    output += '-------------------------------------------------------------------\n'
    output += '--- define different layer types                                ---\n'
    output += '-------------------------------------------------------------------\n'
    ol = i3.write([len(sl.layers)])
    ol = ol.ljust(26)
    output += ol + 'NLTYPE: number of different layer types\n'
    f74x3 = ff.FortranRecordWriter('3F7.4')
    blayers = [l for l in sl.layers if l.isBulk]
    nblayers = [l for l in sl.layers if not l.isBulk]
    for i, layer in enumerate(sl.layers):
        output += '-   layer type '+str(i+1)+' ---\n'
        if layer.isBulk:
            output += '  2                       LAY = 2: layer type no. '+str(i+1)+' has bulk lateral periodicity\n'
        else:
            output += '  1                       LAY = 1: layer type no. '+str(i+1)+' has overlayer lateral periodicity\n'
        if layer.isBulk:
            bulkUnique = []
            for atom in layer.atlist:
                found = False
                for x in range(-1,2):   # sample all surrounding cells to find atoms in bulk cell:
                    for y in range(-1,2):
                        spos = np.dot(rp.SUPERLATTICE,(atom.pos[:2]+np.array([x,y])))
                        if np.amin(spos) >= 0.0 and np.amax(spos) < 1.0:
                            found = True
                if found:
                    bulkUnique.append(atom)
            ol = i3.write([len(bulkUnique)])
            #sanity check: ratio of unit cell areas (given simply by SUPERLATTICE) should match ratio of written vs skipped atoms:
            arearatio = 1/np.linalg.det(rp.SUPERLATTICE)
            atomratio = len(bulkUnique)/len(layer.atlist)
            if abs(arearatio-atomratio) > 1e-3:
                logging.warning('Ratio of bulk atoms inside/outside the bulk unit cell differs from bulk/slab unit cell size ratio. This means that the actual periodicity of the POSCAR does not match the periodicity given in the SUPERLATTICE parameter. Check SUPERLATTICE parameter and bulk symmetry!')
        else:
            ol = i3.write([len(layer.atlist)])
        ol = ol.ljust(26)
        output += ol+'number of Bravais sublayers in layer '+str(i+1)+'\n'
        if layer.isBulk:
            writelist = bulkUnique
        else:
            writelist = layer.atlist
        writelist.sort(key=lambda atom: -atom.pos[2])
        for atom in writelist:
            if not rp.GEO_VERTSTACK:
                writepos = atom.posInLayer
            else:
                writepos = atom.cartpos - np.array([nblayers[-1].cartori[0], nblayers[-1].cartori[1], atom.layer.cartori[2]])
            ol = i3.write([sl.sitelist.index(atom.site)+1]) + f74x3.write([writepos[2], writepos[0], writepos[1]])
            ol = ol.ljust(26)
            output += ol+atom.el+str(atom.oriN)+'\n'
    output += '-------------------------------------------------------------------\n'
    output += '--- define bulk stacking sequence                               ---\n'
    output += '-------------------------------------------------------------------\n'
    output += '  0                       TSLAB = 0: compute bulk using layer doubling\n' #change if this is supposed to be a parameter...
    if rp.ASAZ == tl.DEFAULT:
        # assume that interlayer vector from bottom non-bulk to top bulk layer is the same as between bulk units
        bvectors_ASA = blayers[0].cartori - sl.layers[blayers[0].num-1].cartori
        bvectors_ASA[2] = blayers[0].cartori[2] - sl.layers[blayers[0].num-1].cartbotz
    else:
        zdiff = rp.ASAZ+blayers[0].cartbotz-blayers[0].cartori[2]
        bvectors_ASA = [-sl.ucell[0][2]*zdiff/sl.ucell[2][2], -sl.ucell[1][2]*zdiff/sl.ucell[2][2], rp.ASAZ]
    if rp.BLAY == 2:
        bvectors_ASBULK = blayers[1].cartori - blayers[0].cartori
        bvectors_ASBULK[2] = blayers[1].cartori[2] - blayers[0].cartbotz
        bl2num = blayers[1].num
    else:
        bvectors_ASBULK = bvectors_ASA
        bl2num = blayers[0].num
    ol = f74x3.write([bvectors_ASA[2], bvectors_ASA[0], bvectors_ASA[1]])
    ol = ol.ljust(26)
    output += ol + 'ASA interlayer vector between different bulk units\n'
    ol = i3.write([blayers[0].num+1])
    ol = ol.ljust(26)
    output += ol + 'top layer of bulk unit: layer type '+str(blayers[0].num+1)+'\n'
    ol = i3.write([bl2num+1])
    ol = ol.ljust(26)
    output += ol + 'bottom layer of bulk unit: layer type '+str(bl2num+1)+'\n'
    ol = f74x3.write([bvectors_ASBULK[2], bvectors_ASBULK[0], bvectors_ASBULK[1]])
    ol = ol.ljust(26)
    output += ol + 'ASBULK between the two bulk unit layers (may differ from ASA)\n'
    output += '-------------------------------------------------------------------\n'
    output += '--- define layer stacking sequence and Tensor LEED output       ---\n'
    output += '-------------------------------------------------------------------\n'
    nonbulk = len(sl.layers)-rp.BLAY
    if len(rp.TOUTPUT) < nonbulk:    #check TOUTPUT parameter
        if rp.TOUTPUT:
            logging.warning('Parameters TOUTPUT is defined, but contains less values than there are non-bulk layers. Missing values will be set to 1.')
        for i in range(0,nonbulk-len(rp.TOUTPUT)):
            rp.TOUTPUT.append(1)
    if len(rp.TOUTPUT) > nonbulk:
        logging.warning('Parameters TOUTPUT is defined, but contains more values than there are non-bulk layers. Excess values will be ignored.')
    ol = i3.write([len(sl.layers)-rp.BLAY])
    ol = ol.ljust(26)
    output += ol + 'NSTACK: number of layers stacked onto bulk\n'
    for layer in list(reversed(nblayers)):
        n = layer.num+1
        if layer == nblayers[-1] or not rp.GEO_VERTSTACK:
            v = sl.layers[n].cartori-layer.cartori
        else:
            v = np.array([0.0,0.0,0.0])
        v[2] = sl.layers[n].cartori[2]-layer.cartbotz
        ol = i3.write([n]) + f74x3.write([v[2],v[0],v[1]])
        ol = ol.ljust(26)
        output += ol + 'layer '+str(n)+': layer type '+str(n)+', interlayer vector below\n'     #at the moment, every layer is also a layer type
        ol = i3.write([rp.TOUTPUT[layer.num]])
        ol = ol.ljust(26)
        output += ol+'0/1: Tensor output is required for this layer (TOUTPUT)\n'
        i = 1
        for atom in layer.atlist:
            ol = 'T_'+atom.el+str(atom.oriN)
            ol = ol.ljust(26)
            output += ol + 'Tensor file name, current layer, sublayer '+str(i)+'\n'
            i += 1
    output += '-------------------------------------------------------------------\n'
    output += '--- end geometrical input                                       ---\n'
    output += '-------------------------------------------------------------------\n' 
    
    try:
        wf = open('AUXGEO', 'w')
        wf.write(output)
        wf.close()
    except:
        logging.error("Failed to write AUXGEO:", exc_info=True)
    logging.info("Wrote to AUXGEO successfully.")
   
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