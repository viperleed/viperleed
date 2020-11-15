# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 09:44:38 2020

@author: Florian Kraushofer

Functions for reading and writing beams files: BEAMLIST, IVBEAMS, EXPBEAMS, 
AUXEXPBEAMS
"""

import logging
import numpy as np
import re
import fortranformat as ff

import tleedmlib as tl
from guilib.base import project_to_first_domain

logger = logging.getLogger("tleedm.files.logging")

def readBEAMLIST(filename="_BEAMLIST"):
    """Reads the _BEAMLIST file and returns the contents as a list of the
    lines as strings."""
    beamlist = []
    try:
        with open(filename, "r") as rf:
            beamlist = rf.readlines()
        logger.debug("_BEAMLIST file was read successfully")
    except:
        logger.error("Error opening _BEAMLIST file.")
        raise
    return beamlist

def readIVBEAMS(filename='IVBEAMS'):
    """Reads an IVBEAMS file and returns a list of beams (using Beam class)"""
    # open input file
    try:
        with open(filename, 'r') as rf:
            ivbeamlines = rf.readlines()
    except:
        raise
    #read IVBEAMS:
    linenum = 1		# iterates the current line being read
    hklist = []
    for line in ivbeamlines:
        # ignore brackets and vbars, except as spacers
        line = line.replace("(", " ")
        line = line.replace(")", " ")
        line = line.replace("|", " ")
        llist = line.split()
        if len(llist) == 1 and linenum != 1:
            logger.warning('A line with only one element was found in '
                            'IVBEAMS and will be skipped: '+line)
        elif len(llist) >= 2:
            f = [None,None]
            for i in range(0,2):
                try:
                    f[i] = float(llist[i])
                except ValueError:
                    if '/' in llist[i]:
                        try:
                            f[i] = (float(llist[i].split('/')[0])
                                    / float(llist[i].split('/')[1]))
                        except ValueError:
                            if linenum != 1:
                                logger.error('Error reading IVBEAMS line: '
                                              +line)
                                raise
                    else:
                        if linenum != 1:
                            logger.error('Error reading IVBEAMS line: '+line)
                            raise
            if linenum != 1:
                if not (f[0],f[1]) in hklist and not None in f:
                    hklist.append((f[0],f[1]))
            else:
                # check whether there is data in first line by mistake
                if not None in f:  # data was read
                    logger.warning('It looks like the first line in the '
                        'IVBEAMS file may contain data. Note that the first '
                        'line should be a header line, so the data will not '
                        'be read.')
        linenum += 1
    beams = []
    for hk in hklist:
        beams.append(tl.Beam(hk))
    logger.debug("IVBEAMS file was read successfully")
    return beams

def sortIVBEAMS(sl, rp):
    """Sorts the beams in IVBEAMS such that they appear in the same order as
    in the _BEAMLIST. Returns the sorted list."""
    # read BEAMLIST
    if rp.beamlist == []:
        logger.warning("sortIVBEAMS routine: no beamlist passed, "
            "attempting to read _BEAMLIST directly.")
        try:
            with open('_BEAMLIST', 'r') as rf:
                rp.beamlist = rf.readlines()
        except FileNotFoundError:
            logger.error("_BEAMLIST not found.")
            raise
    err = 1e-3           #since beams are saved as floats, give error tolerance
    symeq = tl.leedbase.getSymEqBeams(sl, rp)
    # first, get beamlist as floats
    blfs = []
    for line in rp.beamlist:
        llist = line.split()
        if len(llist) > 1:
            fl = []
            for i in range(0,2):
                try:
                    fl.append(float(llist[i]))
                except:
                    pass
            if len(fl) == 2:
                blfs.append(fl)
    # now figure out if there are beams in IVBEAMS without a partner
    for (ind, ib) in enumerate(rp.ivbeams):
        found = any([ib.isEqual_hk(lb, eps=err) for lb in blfs])
        if found:
            continue
        eqbl = []   # hk of equivalent beams
        for (hk, i) in symeq:
            if ib.isEqual_hk(hk, eps=err):
                eqbl.extend([hk2 for (hk2, j) in symeq if i == j])
        for eqb in eqbl:
            for lb in blfs:
                if abs(eqb[0]-lb[0]) < err and abs(eqb[1]-lb[1]) < err:
                    # rename the beam
                    if any([ib2.isEqual_hk(eqb, eps=err) for ib2 
                                                          in rp.ivbeams]):
                        logger.debug("Beam " + ib.label + " is not in "
                        "BEAMLIST, but an equivalent beam already is. "
                        "Beam will be dropped.")
                    else:
                        b = tl.Beam(eqb)
                        logger.debug("Beam " + ib.label + " is not in "
                            "BEAMLIST, renaming to equivalent beam "
                             + b.label + ".")
                        rp.ivbeams[ind] = b
                    found = True
                    break
            if found:
                break
        if not found:
            logger.warning('IVBEAMS contains beam ' + ib.label + ', which '
                'was not found in the _BEAMLIST file. Beam will be dropped.')
    # now sort
    ivsorted = []
    for lb in blfs:
        ivsorted.extend([ib for ib in rp.ivbeams if ib.isEqual_hk(lb, 
                                                                  eps=err)])
    return ivsorted

def readOUTBEAMS(filename="EXPBEAMS.csv", sep=";", enrange=[]):
    """Reads beams from an EXPBEAMS.csv or THEOBEAMS.csv file. Returns a list
    of Beam objects. The 'sep' parameter defines the separator. If an energy
    range 'enrange' is passed, beams that contain no data within that range
    will be filtered out before returning. If a Slab and Rparams object are 
    passed, will also check whether any of the experimental beams should be 
    equivalent, and if they are, warn, discard one and raise the halting 
    level."""
    beams = []
    try:
        with open(filename, 'r') as rf:
            lines = [l[:-1] for l in rf.readlines()]
    except FileNotFoundError:
        logger.error("Error reading "+filename)
        raise
    firstline = True
    rgx = re.compile(r'[\*\(\s]*(?P<h>[-0-9/]+)\s*\|\s*(?P<k>[-0-9/]+)')
    for line in lines:
        if firstline and sep not in line:   # try some other separators
            for sep2 in [s for s in [";", ","] if s != sep]:
                if sep2 in line and len(line.split(sep2)) > 2:
                    logger.info("Found separator '"+sep2+"' in "+filename
                                 +", expected '"+sep+"'. Attempting to read "
                                 "with '"+sep2+"'.")
                    sep = sep2
                    break
        llist = line.split(sep)
        if llist[-1].strip() == "":
            llist = llist[:-1]
        if firstline:
            firstline = False
            for label in llist[1:]:
                m = rgx.match(label)
                if m == None:
                    logger.error("readOUTBEAMS: Could not parse h/k in "
                                  "label: "+label)
                    return []
                sh = m.group("h")   #string h
                sk = m.group("k")   #string k
                try:
                    if "/" in sh:
                        h = int(sh.split("/")[0]) / int(sh.split("/")[1])
                    else:
                        h = int(sh)
                    if "/" in sk:
                        k = int(sk.split("/")[0]) / int(sk.split("/")[1])
                    else:
                        k = int(sk)
                except:
                    logger.error("readOUTBEAMS: Could not parse h/k in "
                                    "label: "+label)
                    return []
                beams.append(tl.Beam((h,k)))
        elif len(line) > 1:
            try:
                en = float(llist[0])
            except:
                logger.error("readOUTBEAMS: Could not parse "+llist[0]+"as "
                              "an energy")
                return []
            for i in range(0,len(llist)):
                try:
                    f = float(llist[i+1])
                    if not np.isnan(f):
                        beams[i].intens[en] = f
                except:
                    f = None
    if len(enrange) == 2:
        remlist = []
        for b in beams:
            if ((enrange[0] > 0 and max(b.intens) < enrange[0]) or
                (enrange[1] > 0 and min(b.intens) > enrange[1])):
                # beam has no data in given interval; remove
                remlist.append(b)
                logger.warning("Experimental beam "+b.label+" contains no "
                                "data in the given energy range. The beam "
                                "will be ignored.")
        for b in remlist:
            beams.remove(b)
    totalrange = 0.
    minmax = enrange[:]
    if len(minmax) < 2:
        minmax = [0, 1e5]
    if minmax[0] < 0:
        minmax[0] = 0
    if minmax[1] <= 0:
        minmax[1] = 1e5
    for b in beams:
        totalrange += (min(max(b.intens), minmax[1]) 
                       - max(min(b.intens), minmax[0]))
    logger.info("Loaded "+filename+" file containing {} beams (total energy "
                "range: {:.2g} eV).".format(len(beams), totalrange))
    return beams

def checkEXPBEAMS(sl, rp):
    remlist = []
    symeq = tl.leedbase.getSymEqBeams(sl, rp)
    for (bi, b) in enumerate(rp.expbeams):
        if b in remlist:
            continue
        eqbl = []   # hk of equivalent beams
        for (hk, i) in symeq:
            if b.isEqual_hk(hk):
                eqbl.extend([hk2 for (hk2, j) in symeq if i == j])
                break
        for b2 in rp.expbeams[bi+1:]:
            for hk in eqbl:
                if b2.isEqual_hk(hk):
                    remlist.append(b2)
                    w = ("Experimental beam "+b2.label+" is "
                        "equivalent to experimental beam "+b.label+". ")
                    if rp.HALTING > 2:
                        w += "Beam "+b2.label+" will be discarded."
                    else:
                        w += ("Remove one of them, or average them as "
                              "appropriate.")
                    logger.warning(w)
                    rp.setHaltingLevel(2)
                    break
    for b in remlist:
        rp.expbeams.remove(b)

def readAUXEXPBEAMS(filename="AUXEXPBEAMS", interactive=False):
    """Reads beams from an AUXEXPBEAMS file, which already has the formatting
    required by TensErLEED. Returns a list of Beam objects. Only works if the
    comment lines are of format *( h k )*, including the stars, with h and k
    separated by whitespace. If 'interactive' is True, then the user will be 
    asked to clarify labels that cannot be read."""
    expbeams = []
    try:
        with open(filename, 'r') as rf:
            lines = [l[:-1] for l in rf.readlines()]
    except FileNotFoundError:
        logger.error("Error reading AUXEXPBEAMS.")
        raise
    read = False
    rf62x12 = ff.FortranRecordReader('12F6.2')
    rgx = re.compile(r'[\*\(\s]*(?P<h>[-0-9/]+)\s+(?P<k>[-0-9/]+)')
    for line in lines:
        if "*" in line:
            read = True
            topline = True  #next line contains number of beams and scaling
            failedToRead = False
            m = rgx.match(line)
            if m != None:
                sh = m.group("h")   #string h
                sk = m.group("k")   #string k
                try:
                    h = tl.base.parseMathSqrt(sh)
                    k = tl.base.parseMathSqrt(sk)
                except:
                    failedToRead = True
                else:
                    newbeam = tl.Beam((h, k))
                    expbeams.append(newbeam)
            else:
                failedToRead = True
            if failedToRead:
                if not interactive:
                    logger.error("readAUXEXPBEAMS: Could not parse h/k in "
                                  "line:\n"+line)
                    return []
                else:
                    print("readAUXEXPBEAMS: Could not parse h/k in "
                                  "line:\n"+line)
                    while True:
                        try:
                            hks = input("Please manually input the beam "
                                        "indices (space-separated): ")
                        except:
                            return []
                        if hks.lower() in ["quit", "exit"]:
                            return []
                        if hks and len(hks.split()) > 1:
                            try:
                                h = tl.base.parseMathSqrt(hks.split()[0])
                                k = tl.base.parseMathSqrt(hks.split()[1])
                            except:
                                print("Could not parse h/k")
                            else:
                                newbeam = tl.Beam((h, k))
                                expbeams.append(newbeam)
                                break
        elif read and topline:
            topline = False
            llist = line.split()
            try:
                scaling = float(llist[1])
            except:
                logger.error("readAUXEXPBEAMS: Could not parse number of "
                              "beams or scaling factor in line:\n"+line)
                return []
        elif read:
            newvals = rf62x12.read(line)
            i = 0
            while i+1 < len(newvals):
                if newvals[i] == None:
                    break
                newbeam.intens[newvals[i]] = newvals[i+1] / scaling
                i += 2
    return expbeams

def writeIVBEAMS(sl, rp, filename="IVBEAMS"):
    """Writes an IVBEAMS file based on rp.exbeams. Returns 
    those beams in IVBEAMS form."""
    d = tl.leedbase.getLEEDdict(sl, rp)
    if d is None:
        logger.error("Failed to write IVBEAMS")
        return []
    makebeams = project_to_first_domain([b.hkfrac for b in rp.expbeams], d)
    output = "This IVBEAMS file was automatically generated from EXPBEAMS\n"
    ivbeams = [tl.Beam(hk) for hk in makebeams]
    for b in ivbeams:
        output += "{: 10.6f} {: 10.6f}\n".format(b.hk[0], b.hk[1])
    try:
        with open(filename, "w") as wf:
            wf.write(output)
        logger.debug("Wrote IVBEAMS file successfully.")
    except:
        logger.error("Exception while writing IVBEAMS file: ", 
                      exc_info=True)
        raise
    return ivbeams

def writeOUTBEAMS(beams, filename="THEOBEAMS.csv", sep="; "):
    """Takes a list of Beam objects and writes them to a comma-separated 
    file."""
    nan = "NaN" #what to put when no value
    output = "E".rjust(7)+sep
    w = max(11, beams[0].lwidth*2 + 3)
    energies = []
    for b in beams:
        output += b.label.rjust(w)+sep
        energies.extend([k for k in b.intens if k not in energies])
    output = output[:-len(sep)]
    output += "\n"
    energies.sort()
    for en in energies:
#        output += str(round(en,2)).rjust(7)+sep
        output += '{:7.2f}'.format(en) + sep
        for b in beams:
            if en in b.intens:
                output += '{:0.5E}'.format(b.intens[en]).rjust(w) + sep
            else:
                output += nan.rjust(w) + sep
        output = output[:-len(sep)]
        output += "\n"
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except:
        logger.error("Failed to write "+filename)
        raise
    logger.debug("Wrote to "+filename+" successfully")
    return 0

def writeAUXBEAMS(ivbeams=[], beamlist=[], beamsfile='IVBEAMS', 
                  readfile='_BEAMLIST', writefile='AUXBEAMS', write=True):
    """"Reads from a _BEAMLIST file (full list of beams for calculation), 
    finds the beams listed in IVBEAMS (if not passed as a list 'beams', will 
    attempt to call readIVBEAMS directly) and copies the corresponding lines 
    to AUXBEAMS. Returns a list of the corresponding beam numbers in _BEAMLIST.
    """
    if not readfile == "_BEAMLIST":
        blstr = "Beam list (filename "+readfile+")"
    else:
        blstr = "_BEAMLIST"
        
    if ivbeams == []:         #if 'ivbeams' is empty, try to fill it
        ivbeams = readIVBEAMS(beamsfile)
        
    output = '   1               IFORM\n'      # !!! WHAT IS THIS VALUE? WHERE TO GET IT FROM?
    
    # read BEAMLIST
    if beamlist == []:
        logger.warning("writeAUXBEAMS routine: no beamlist passed, "
            "attempting to read _BEAMLIST directly.")
        try:
            beamlist = readBEAMLIST(readfile)
        except:
            logger.error("Error getting beamlist not found.")
            raise
    
    err = 1e-4           #since beams are saved as floats, give some error 
                         #  tolerance when comparing e.g. 1/3
    #read BEAMLIST - very little error handling here, I'm assuming 
    #  the beamlists are safe; could be added later
    foundbeams = []
    numlist = []
    blocks = 0
    totalbeams = 0
    for line in beamlist:
        llist = line.split()
        if len(llist) > 1:
            totalbeams += 1
            for b in ivbeams:
                if b.isEqual_hk((float(llist[0]), float(llist[1])), eps=err):
                    output += line
                    foundbeams.append(b)
                    numlist.append(int(line.split('.')[-1]))
        else:
            blocks += 1
    for beam in [b for b in ivbeams if b not in foundbeams]:
        logger.warning('IVBEAMS contains beam '+ beam.label +', which was '
                        'not found in '+blstr)
    #write output
    if write:
        if not writefile == 'AUXBEAMS':
            wfstr = 'AUXBEAMS file (filename '+writefile+')'
        else:
            wfstr = 'AUXBEAMS'
        try:
            with open(writefile, 'w') as wf:
                wf.write(output)
        except:
            logger.error("Failed to write "+wfstr+" file")
            raise
        logger.debug("Wrote to "+wfstr+" successfully.")
    return numlist, blocks, totalbeams

def writeAUXEXPBEAMS(beams, filename="AUXEXPBEAMS", header="Unknown system", 
                     write=True, numbers=False):
    """Takes a list of Beam objects and writes them in the format required by 
    TensErLEED for experimental beam data. Returns the whole output as a 
    string. 'numbers' defines whether a sequence of beam numbers in line 2 
    is expected."""
    output = header+"\n"
    if numbers:
        i3x25 = ff.FortranRecordWriter("25I3")
        output += i3x25.write([n+1 for n in range(0, len(beams))]) + "\n"
    output += " (12F6.2)\n"
    f62x12 = ff.FortranRecordWriter('12F6.2')
    i4 = ff.FortranRecordWriter('I4')
    for beam in beams:
        # renormalize
        minintens = min(beam.intens.values())
        offset = max(0, -minintens)  # if beams contain negative values, offset 
        scaling = 999.99 / (max(beam.intens.values()) + offset)
        for k in beam.intens:
            beam.intens[k] = (beam.intens[k] + offset) * scaling
        # write
        output += "*"+beam.label.replace("|"," ") + "*\n"
        ol = i4.write([len(beam.intens)]).ljust(8)
        output += ol + '{:10.4E}\n'.format(scaling)
        outlist = [val for tup in beam.intens.items() for val in tup]
                             # zip & flatten energies and values into one list
        output += f62x12.write(outlist)
        if output[-1] != "\n":
            output += "\n"
    #write output
    if write:
        try:
            with open(filename, 'w') as wf:
                wf.write(output)
        except:
            logger.warning("Failed to write "+filename)
        logger.debug("Wrote to "+filename+" successfully")
    return output