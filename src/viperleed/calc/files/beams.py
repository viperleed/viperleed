"""Functions for reading and writing beams files.

This includes: BEAMLIST, IVBEAMS, EXPBEAMS and AUXEXPBEAMS.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-19'
__license__ = 'GPLv3+'

import copy
from io import StringIO
import logging
import os
import re

import fortranformat as ff
import numpy as np

from viperleed.calc.classes.beam import Beam
from viperleed.calc.lib import leedbase
from viperleed.calc.lib.base import parseMathSqrt
from viperleed.guilib import project_to_first_domain


logger = logging.getLogger(__name__)


def averageBeams(beams, weights=None):
    """Takes a list of percentages and a list of lists of Beam objects.
    Returns a new list of Beam objects with weighted averaged intensities."""
    if beams is None:
        raise ValueError("averageBeams: No beams passed.")
    if weights is None or len(weights) == 0:
        weights = [1/len(beams)] * len(beams)
    avbeams = copy.deepcopy(beams[0])
    for (i, b) in enumerate(avbeams):
        if not all([beams[j][i].isEqual(b) for j in range(1, len(beams))]):
            _err = "averageBeams: Beam lists are mismatched."
            logger.error(_err)
            raise ValueError(_err)
            return []
        if not all([set(beams[j][i].intens.keys()) == set(b.intens.keys())
                    for j in range(1, len(beams))]):
            _err = "averageBeams: Beams have different energy ranges."
            logger.error(_err)
            raise ValueError(_err)
            return []
        for en in b.intens:
            b.intens[en] *= weights[0]
            for j in range(1, len(beams)):
                b.intens[en] += beams[j][i].intens[en] * weights[j]
    return avbeams


def readBEAMLIST(filename="BEAMLIST"):
    """Reads the BEAMLIST file and returns the contents as a list of the
    lines as strings."""
    beamlist = []
    try:
        with open(filename, "r") as rf:
            beamlist = rf.readlines()
        logger.debug("BEAMLIST file was read successfully")
    except Exception:
        logger.error("Error opening BEAMLIST file.")
        raise
    return beamlist


def readIVBEAMS(filename='IVBEAMS'):
    """Reads an IVBEAMS file and returns a list of beams (using Beam class)"""
    # open input file
    try:
        with open(filename, 'r') as rf:
            ivbeamlines = rf.readlines()
    except Exception:
        raise
    linenum = 1		# iterates the current line being read
    hklist = []
    # TODO: skip commented out lines (same as PARAMETERS, etc.)
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
            f = [None, None]
            for i in range(0, 2):
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
                                             + line)
                                raise
                    else:
                        if linenum != 1:
                            logger.error('Error reading IVBEAMS line: '+line)
                            raise
            if linenum != 1:
                if not (f[0], f[1]) in hklist and None not in f:
                    hklist.append((f[0], f[1]))
            else:
                # check whether there is data in first line by mistake
                if None not in f:  # data was read
                    logger.warning(
                        'It looks like the first line in the '
                        'IVBEAMS file may contain data. Note that the first '
                        'line should be a header line, so the data will not '
                        'be read.')
        linenum += 1
    beams = []
    for hk in hklist:
        beams.append(Beam(hk))
    logger.debug("IVBEAMS file was read successfully")
    return beams


def sortIVBEAMS(sl, rp):
    """Sorts the beams in IVBEAMS such that they appear in the same order as
    in the BEAMLIST. Returns the sorted list."""
    # read BEAMLIST
    if rp.beamlist == []:
        logger.warning("sortIVBEAMS routine: no beamlist passed, "
                       "attempting to read BEAMLIST directly.")
        try:
            with open('BEAMLIST', 'r') as rf:
                rp.beamlist = rf.readlines()
        except FileNotFoundError:
            logger.error("BEAMLIST not found.")
            raise
    err = 1e-3          # since beams are saved as floats, give error tolerance
    symeq = leedbase.getSymEqBeams(sl, rp)
    # first, get beamlist as floats
    blfs = []
    for line in rp.beamlist:
        llist = line.split()
        if len(llist) > 1:
            fl = []
            for i in range(0, 2):
                try:
                    fl.append(float(llist[i]))
                except (ValueError, IndexError):
                    pass
            if len(fl) == 2:
                blfs.append(fl)
    # now figure out if there are beams in IVBEAMS without a partner
    not_found = []
    found_equivalent = []
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
                    if any([ib2.isEqual_hk(eqb, eps=err)
                            for ib2 in rp.ivbeams]):
                        found_equivalent.append(ib.label)
                    else:
                        b = Beam(eqb)
                        logger.debug(
                            "Beam " + ib.label + " is not in "
                            "BEAMLIST, renaming to equivalent beam "
                            + b.label + ".")
                        rp.ivbeams[ind] = b
                    found = True
                    break
            if found:
                break
        if not found:
            not_found.append(ib.label)
    if found_equivalent:
        logger.debug(
            "The following beams from IVBEAMS will be dropped because they "
            "are not in BEAMLIST, but equivalent beams already are: "
            + ", ".join(found_equivalent))
    if not_found:
        logger.warning(
            "The following beams from IVBEAMS are not in BEAMLIST and will be "
            "dropped: " + ", ".join(not_found))
    # now sort
    ivsorted = []
    for lb in blfs:
        ivsorted.extend([ib for ib in rp.ivbeams if ib.isEqual_hk(lb,
                                                                  eps=err)])
    return ivsorted


def readOUTBEAMS(filename="EXPBEAMS.csv", sep=",", enrange=None):
    #TODO: this need some optimizing; use csv and dict reader class
    """Reads beams from an EXPBEAMS.csv or THEOBEAMS.csv file.
    
    The beams are returned as a list of Beam objects.
    
    Parameters
    ----------
    filename : str, StringIO, optional
        If string, interpret as the path to the file to read. If
        a StringIO object, assume it contains already the contents
        of the file to be interpreted. Default is "EXPBEAMS.csv".
    sep: str, optional
        Separator character for splitting beams. Default is ",".
    enrange : Sequence, or None, optional
        Energy range to be used to filter beams. Beams that contain no
        data within the range will be filtered out before returning.
        If given and not None, it should be a Sequence with two elements
        (min and max).
    
    Returns
    -------
    beams : list of Beam
        Beam objects with info read from the input given.
    """
    beams = []
    # Deals with input in the form of a StingIO (can be used to feed file input as a string)
    if isinstance(filename, StringIO):
        lines = (li[:-1] for li in filename)
        filename = "StringIO"
    else:
        try:
            with open(filename, 'r') as rf:
                lines = [li[:-1] for li in rf.readlines()]
        except FileNotFoundError:
            _filename = str(filename)
            if _filename.endswith('.csv') and os.path.isfile(_filename[:-4]):
                with open(_filename[:-4], 'r') as rf:
                    lines = [li[:-1] for li in rf.readlines()]
            else:
                logger.error("Error reading "+_filename)
                raise

    firstline = True
    rgx = re.compile(r'[\*\(\s]*(?P<h>[-0-9/]+)\s*\|\s*(?P<k>[-0-9/]+)')
    for line in lines:
        if firstline and sep not in line:   # try some other separators
            for sep2 in [s for s in [",", ";"] if s != sep]:
                if sep2 in line and len(line.split(sep2)) > 2:
                    logger.info("Found separator '"+sep2+"' in "+filename
                                + ", expected '"+sep+"'. Attempting to read "
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
                if m is None:
                    logger.error("readOUTBEAMS: Could not parse h/k in "
                                 "label: "+label)
                    return []
                sh = m.group("h")   # string h
                sk = m.group("k")   # string k
                try:
                    if "/" in sh: # TODO if fractionals found - pass fractions already to beam
                        h = int(sh.split("/")[0]) / int(sh.split("/")[1])
                    else:
                        h = float(sh)
                    if "/" in sk:
                        k = int(sk.split("/")[0]) / int(sk.split("/")[1])
                    else:
                        k = float(sk)
                except (ValueError, IndexError, ZeroDivisionError):
                    logger.error("readOUTBEAMS: Could not parse h/k in "
                                 "label: "+label)
                    return []
                beams.append(Beam((h, k)))
        elif len(line) > 1:
            try:
                en = float(llist[0])
            except ValueError:
                logger.error("readOUTBEAMS: Could not parse "+llist[0]+"as "
                             "an energy")
                return []
            for i in range(0, len(llist)):
                try:
                    f = float(llist[i+1])
                    if not np.isnan(f):
                        beams[i].intens[en] = f
                except (ValueError, IndexError):
                    pass
                
    # cleanup of beams
    threshold_intes = 1e-8  # threshold value of TensErLEED (rfacsb.f in subroutine grid)
    lowest_real_en = 0
    for beam in beams:

        unfiltered_intens = copy.deepcopy(list(beam.intens.items())) 
        # find first intensity > threshold
        for en, intensity in unfiltered_intens:
            if intensity > threshold_intes:
                lowest_real_en = en
                break

        # remove parts of beam that are set to NaN or below threshhold
        for en, intensity in unfiltered_intens:
            if (np.isnan(intensity) or (en < lowest_real_en)):
                beam.intens.pop(en)

    if enrange is not None and len(enrange) == 2:
        remlist = []
        for b in beams:
            if (len(b.intens) == 0 or
                    (enrange[0] > 0 and max(b.intens) < enrange[0]) or
                    (enrange[1] > 0 and min(b.intens) > enrange[1])):
                # beam has no data in given interval; remove
                remlist.append(b)
        if remlist:
            logger.warning("The following beams contain no data in the given "
                           "energy range and will be ignored: "
                           + ", ".join([b.label for b in remlist]))
        for b in remlist:
            beams.remove(b)

    # Add an offset such that intensities are always positive, and warn.
    for beam in beams:
        min_intensity = min((0, *beam.intens.values()))
        if min_intensity < 0:
            logger.warning(f"Negative intensity encountered in beam {beam.label} while reading {filename}."
                           f" An offset was added so that the minimum intensity of this beam is 0.")
            for en in beam.intens.keys():
                beam.intens[en] += abs(min_intensity)

    if enrange is None:
        return beams
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


def checkEXPBEAMS(sl, rp, domains=False):
    remlist = []
    symeq = leedbase.getSymEqBeams(sl, rp)
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
            lines = [li[:-1] for li in rf.readlines()]
    except FileNotFoundError:
        logger.error("Error reading AUXEXPBEAMS.")
        raise
    read = False
    rf62x12 = ff.FortranRecordReader('12F6.2')
    rgx = re.compile(r'[\*\(\s]*(?P<h>[-0-9/]+)\s+(?P<k>[-0-9/]+)')
    for line in lines:
        if "*" in line:
            read = True
            topline = True  # next line contains number of beams and scaling
            failedToRead = False
            m = rgx.match(line)
            if m is not None:
                sh = m.group("h")   # string h
                sk = m.group("k")   # string k
                try:
                    h = parseMathSqrt(sh)
                    k = parseMathSqrt(sk)
                except Exception:
                    failedToRead = True
                else:
                    newbeam = Beam((h, k))
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
                        except Exception:
                            return []
                        if hks.lower() in ["quit", "exit"]:
                            return []
                        if hks and len(hks.split()) > 1:
                            try:
                                h = parseMathSqrt(hks.split()[0])
                                k = parseMathSqrt(hks.split()[1])
                            except Exception:
                                print("Could not parse h/k")
                            else:
                                newbeam = Beam((h, k))
                                expbeams.append(newbeam)
                                break
        elif read and topline:
            topline = False
            llist = line.split()
            try:
                scaling = float(llist[1])
            except Exception:
                logger.error("readAUXEXPBEAMS: Could not parse number of "
                             "beams or scaling factor in line:\n"+line)
                return []
        elif read:
            newvals = rf62x12.read(line)
            i = 0
            while i+1 < len(newvals):
                if newvals[i] is None:
                    break
                newbeam.intens[newvals[i]] = newvals[i+1] / scaling
                i += 2
    return expbeams


def writeIVBEAMS(sl, rp, filename="IVBEAMS", domains=False):
    """Writes an IVBEAMS file based on rp.exbeams. Returns
    those beams in IVBEAMS form."""
    if not domains:
        d = [leedbase.getLEEDdict(sl, rp)]
    else:
        d = [leedbase.getLEEDdict(dp.sl, dp.rp) for dp in rp.domainParams]
    if any([v is None for v in d]):
        logger.error("Failed to write IVBEAMS")
        return []
    output = ("  h         k          Beams to calculate, automatically "
              "generated from EXPBEAMS\n")
    if rp.AVERAGE_BEAMS is not False:
        makebeams = project_to_first_domain([b.hkfrac for b in rp.expbeams],
                                            *d)
        ivbeams = [Beam(hk) for hk in makebeams]
    else:
        ivbeams = [Beam(b.hkfrac) for b in rp.expbeams]
    for b in ivbeams:
        output += "{: 10.6f} {: 10.6f}\n".format(b.hk[0], b.hk[1])
    try:
        with open(filename, "w") as wf:
            wf.write(output)
        logger.debug("Wrote IVBEAMS file successfully.")
    except Exception:
        logger.error("Exception while writing IVBEAMS file: ",
                     exc_info=True)
        raise
    return ivbeams


def writeOUTBEAMS(beams, filename="THEOBEAMS.csv", sep=", ",
                  which="intensity"):
    """Takes a list of Beam objects and writes them to a comma-separated
    file. Parameter 'which' defines what to write; set to 'amp_real' or
    'amp_imag' to write real and imaginary parts of complex amplitudes."""
    nan = "NaN"  # what to put when no value
    output = "E".rjust(7)+sep
    if which == "intensity":
        minwidth = 11
    else:
        minwidth = 17
    w = max(minwidth, len(beams[0].label))
    energies = []
    for b in beams:
        output += b.label.rjust(w)+sep
        if which == "intensity":
            energies.extend([k for k in b.intens if k not in energies])
        else:
            energies.extend([k for k in b.complex_amplitude
                             if k not in energies])
    output = output[:-len(sep)]
    output += "\n"
    energies.sort()
    for en in energies:
        output += '{:7.2f}'.format(en) + sep
        for b in beams:
            if which == "intensity":
                write_dict = b.intens
            else:
                write_dict = b.complex_amplitude
            if en in write_dict:
                if which == "intensity":
                    output += '{:0.5E}'.format(b.intens[en]).rjust(w) + sep
                else:
                    if which == "amp_real":
                        v = b.complex_amplitude[en].real
                    else:
                        v = b.complex_amplitude[en].imag
                    output += '{:0.10E}'.format(v).rjust(w) + sep
            else:
                output += nan.rjust(w) + sep
        output = output[:-len(sep)]
        output += "\n"
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except OSError:
        logger.error("Failed to write "+filename)
        raise
    logger.debug("Wrote to "+filename+" successfully")
    return


def writeAUXBEAMS(ivbeams=None, beamlist=None, beamsfile='IVBEAMS',
                  readfile='BEAMLIST', writefile='AUXBEAMS', write=True):
    """"Reads from a BEAMLIST file (full list of beams for calculation),
    finds the beams listed in IVBEAMS (if not passed as a list 'beams', will
    attempt to call readIVBEAMS directly) and copies the corresponding lines
    to AUXBEAMS. Returns a list of the corresponding beam numbers in BEAMLIST.
    """
    if not readfile == "BEAMLIST":
        blstr = "Beam list (filename "+readfile+")"
    else:
        blstr = "BEAMLIST"

    if ivbeams is None:         # if 'ivbeams' is empty, try to fill it
        ivbeams = readIVBEAMS(beamsfile)

    output = '   1               IFORM\n'  # always use formatted input+output

    # read BEAMLIST
    if beamlist is None:
        logger.warning("writeAUXBEAMS routine: no beamlist passed, "
                       "attempting to read BEAMLIST directly.")
        try:
            beamlist = readBEAMLIST(readfile)
        except Exception:
            logger.error("Error getting beamlist not found.")
            raise
    # since beams are saved as floats, give some error tolerance
    err = 1e-4
    # read BEAMLIST - very little error handling here, I'm assuming
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
    not_found = [b.label for b in ivbeams if b not in foundbeams]
    if not_found:
        logger.warning('IVBEAMS contains beams that are not in ' + blstr + ": "
                       + ", ".join(not_found))
    if write:
        if not writefile == 'AUXBEAMS':
            wfstr = 'AUXBEAMS file (filename '+writefile+')'
        else:
            wfstr = 'AUXBEAMS'
        try:
            with open(writefile, 'w') as wf:
                wf.write(output)
        except Exception:
            logger.error("Failed to write "+wfstr+" file")
            raise
        logger.debug("Wrote to "+wfstr+" successfully.")
    return numlist, blocks, totalbeams


def writeAUXEXPBEAMS(beams, filename="AUXEXPBEAMS", header="Unknown system",
                     write=True, numbers=False, version=0, output_format='12E14.5'):
    """Takes a list of Beam objects and writes them in the format required by
    TensErLEED for experimental beam data. Returns the whole output as a
    string. 'numbers' defines whether a sequence of beam numbers in line 2
    is expected. 'version' is to pass the TensErLEED-version, only needed for
    new formatting styles."""
    output = header+"\n"
    if numbers:
        if version < 1.7:
            bformatter = ff.FortranRecordWriter("25I3")
        else:
            bformatter = ff.FortranRecordWriter("25I4")
        output += bformatter.write([n+1 for n in range(0, len(beams))]) + "\n"
        if len(beams) % 25 == 0:
            output += "\n"
    # Standard format used to be 12F6.2, but this may 
    # have too little dynamic range for theory-theory comparisons.
    # Instead we can use 12E14.5 which is the same as used everywhere else.
    if output_format not in ("12F6.2", "12E14.5"):
        logger.error(f"Incompatible format requested for AUXEXPEBAMS: {output_format}\n"
                     f"Possible are '12F6.2' and '12E14.5'.")
        raise RuntimeError
    
    output += f" ({output_format})\n"
    decimal_writer = ff.FortranRecordWriter(output_format)
    integer_writer = ff.FortranRecordWriter('I4')
    for beam in beams:
        # renormalize
        minintens = min(beam.intens.values())
        # Alex March 2022: This is quite important actually!!! Will mess up otherwise if negative values in EXPBEAMS
        offset = max(0, -minintens)  # if beams contain negative values, offset
        scaling = 999.99 / (max(beam.intens.values()) + offset)
        # working dict to hold normalized beams
        working_beam = {}
        for k in beam.intens:
            working_beam[k] = (beam.intens[k] + offset) * scaling
        # write
        output += "*"+beam.label.replace("|", " ") + "*\n"
        ol = integer_writer.write([len(working_beam)]).ljust(8)
        output += ol + '{:10.4E}\n'.format(scaling)
        # zip & flatten energies and values into one list
        outlist = [val for tup in working_beam.items() for val in tup]
        output += decimal_writer.write(outlist)
        if output[-1] != "\n":
            output += "\n"
    if write:
        try:
            with open(filename, 'w') as wf:
                wf.write(output)
        except Exception:
            logger.warning("Failed to write "+filename)
        logger.debug("Wrote to "+filename+" successfully")
    return output


def writeFdOut(beams, beamlist=None, filename="refcalc-fd.out",
               header="Unknown system"):
    """Writes an fd.out file for the given beams. Also returns the file "
    "contents as string."""
    out = header+"\n"
    # read BEAMLIST
    if beamlist is None:
        logger.warning("writeFdOut: no beamlist passed, attempting to read "
                       "BEAMLIST directly.")
        try:
            beamlist = readBEAMLIST("BEAMLIST")
        except Exception:
            logger.error("Error getting beamlist: file not found.")
            raise
    out += "{:>3}\n".format(len(beams))
    energies = set()
    for b in beams:
        for (i, line) in enumerate(beamlist[1:]):
            llist = line.split()
            if (len(llist) > 1
                    and b.isEqual_hk((float(llist[0]), float(llist[1])))):
                out += "{:>4}".format(i+1) + line[:23] + "\n"
                break
        else:
            logger.error("writeFdOut: passed beams contain beam {}, which "
                         "was not found in BEAMLIST.".format(b.label))
            raise ValueError("writeFdOut: beam {} is not in BEAMLIST."
                             .format(b.label))
        energies.update(b.intens.keys())
    for en in sorted(list(energies)):
        o = "{:7.2f} 0.0001".format(en)
        line_len = 1
        for b in beams:
            if en in b.intens:
                o += "{:14.5E}".format(b.intens[en])
            else:
                o += "{:14.5E}".format(0)
            line_len += 1
            if line_len == 5:
                o += "\n"
                line_len = 0
        out += o+"\n"
    try:
        with open(filename, "w") as wf:
            wf.write(out)
    except Exception:
        logger.warning("writeFdOut: Failed to write output to file {}."
                       .format(filename))
        raise
    return out
