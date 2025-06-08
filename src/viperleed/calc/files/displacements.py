"""Functions for reading and interpreting the DISPLACEMENTS file."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-08-19'
__license__ = 'GPLv3+'

import logging
import re

import numpy as np

from viperleed.calc.lib.base import splitSublists
from viperleed.calc.lib.string_utils import read_int_range
from viperleed.calc.symmetry import enforceSymmetry
from viperleed.calc.symmetry import setSymmetry


logger = logging.getLogger(__name__)

class DisplacementsError(Exception):
    """Base exception for DISPLACEMENTS-related errors."""


class NoDisplacementsError(DisplacementsError):
    """No valid displacements were found in DISPLACEMENTS."""


def readDISPLACEMENTS(rp, filename="DISPLACEMENTS"):                            # TODO: This should probably set rp.fileLoaded['DISPLACEMENTS'] = True
    """
    Reads the DISPLACEMENTS file to rp.disp_blocks without interpreting it.

    Parameters
    ----------
    rp : Rparams
        The run parameters object.
    filename : str, optional
        The file name to read from. The default is "DISPLACEMENTS".

    Raises
    ------
    RuntimeError
        Raised if unexpected syntax is encountered and cannot be resolved.

    Returns
    -------
    None.

    """
    rp.disp_blocks = []
    copyblock = None
    loopStarts = []
    loopLine = False
    try:
        with open(filename, 'r') as rf:
            lines = rf.readlines()
    except FileNotFoundError:
        logger.error("DISPLACEMENTS not found.")
        raise
    except Exception:
        logger.error("Error reading DISPLACEMENTS file.")
        raise
    for (i, line) in enumerate(lines):
        for c in ["!", "#", "%"]:    # start of comment
            line = line.split(c)[0].strip()
        if re.search(r'<\s*/?\s*loop\s*>', line.lower()):
            if "=" in line:
                logger.warning("DISPLACEMENTS file: One line cannot contain "
                               "both an '=' sign and a <loop> flag.")
                rp.setHaltingLevel(2)
                continue
            loopLine = True
            ltags = re.findall(r'<\s*/?\s*loop\s*>', line.lower())
            n_blocks_real = len(rp.disp_blocks)
            if copyblock:
                n_blocks_real += 1
            for t in ltags:
                if re.match(r'<\s*loop\s*>', t):
                    loopStarts.append(n_blocks_real)
                else:
                    if len(loopStarts) == 0:
                        logger.warning("DISPLACEMENTS file: Unmatched "
                                       "</loop> flag found, skipping.")
                        rp.setHaltingLevel(2)
                        continue
                    if loopStarts[-1] == n_blocks_real:
                        logger.warning("DISPLACEMENTS file: Found empty "
                                       "loop.")
                        loopStarts.pop()
                        rp.setHaltingLevel(1)
                        continue
                    rp.disp_loops.append((loopStarts.pop(), n_blocks_real - 1))
            continue
        if '=' not in line:
            continue
        if line.startswith('=='):
            if re.match(r'==\s*s', line.lower()):  # search block
                loopLine = False
                try:
                    name = (line[re.match(r'==\s*s(earch)?\s+',
                                          line.lower()).span()[1]:].strip())
                except Exception:
                    name = ""
                names = [n for (_, n) in rp.disp_blocks]
                names.extend([n[:-6] for (_, n) in rp.disp_blocks
                              if n.endswith(" [0 1]") or n.endswith(" [1 0]")])
                if not name or name in names:  # get unique name
                    if not name:
                        i = 1
                    else:
                        i = 2
                    nn = (name + " {}".format(i)).strip()
                    while nn in names:
                        i += 1
                        nn = (name + " {}".format(i)).strip()
                    name = nn
                if not rp.disp_blocks or len(rp.disp_blocks[-1][0]) != 0:
                    if copyblock:
                        rp.disp_blocks.append(copyblock)
                        copyblock = None
                    rp.disp_blocks.append(([], name))
                else:
                    rp.disp_blocks[-1] = ([], name)
                continue
            elif line.split('==')[1].strip().lower()[0] != "d":
                logger.warning("Found DISPLACEMENTS line starting with '==',"
                               "but could not interpret it as the start of a "
                               "search or domain block.")
                rp.setHaltingLevel(2)
                continue
        elif loopLine and not len(rp.disp_blocks) == 0:
            logger.error("DISPLACEMENTS file line {}: Line with <loop> or "
                         "</loop> flags must be followed by a new search "
                         "block or end of file!".format(i+1))
            raise RuntimeError("Misplaced <loop> flag in DISPLACEMENTS file.")
        loopLine = False
        if len(rp.disp_blocks) == 0:
            rp.disp_blocks.append(([], ""))
        p = line.split("=")[0]
        if ("xy" in p or "ab" in p) and "[" not in p:  # shorthand for 2 blocks
            if not copyblock:
                if rp.disp_blocks[-1][1]:
                    names = (rp.disp_blocks[-1][1]+" [1 0]",
                             rp.disp_blocks[-1][1]+" [0 1]")
                else:
                    names = ("", "")
                rp.disp_blocks[-1] = (rp.disp_blocks[-1][0], names[0])
                copyblock = (rp.disp_blocks[-1][0][:], names[1])
            if "xy" in p:
                s = "xy"
            else:
                s = "ab"
            copyblock[0].append(line.replace(s, s+"[0 1]", 1))
            line = line.replace(s, s+"[1 0]", 1)

        elif copyblock:
            copyblock[0].append(line)
        rp.disp_blocks[-1][0].append(line)
    if copyblock:
        rp.disp_blocks.append(copyblock)
    if len(rp.disp_blocks[-1][0]) == 0 and len(rp.disp_blocks) > 1:
        rp.disp_blocks = rp.disp_blocks[:-1]
    if len(loopStarts) != 0:
        logger.warning("DISPLACEMENTS file: Unmatched <loop> flags found, "
                       "loops are still open at end of file.")
        rp.setHaltingLevel(2)
    if not rp.domainParams:
        return
    # in case of domains, now split blocks to domains
    for dp in rp.domainParams:
        dp.rpars.disp_blocks = []
    for (lines, blockname) in rp.disp_blocks:
        d = ""
        dlines = {}
        dnum = 1
        names = [dp.name for dp in rp.domainParams]
        skipblock = False
        for name in names:
            dlines[name] = []
        for line in lines:
            if '=' not in line:
                continue
            if re.match(r'==\s*d', line.lower()):  # start domain block
                skipblock = False
                try:
                    d = (line[re.match(r'==\s*d(omain)?\s+',
                                       line.lower()).span()[1]:].strip())
                except Exception:
                    if str(dnum) in names:
                        d = str(dnum)
                        dnum += 1
                    else:
                        logger.warning(
                            "DISPLACEMENTS file: found start of domain block "
                            "with missing domain name. Block will be skipped.")
                        rp.setHaltingLevel(2)
                        d = ""
                        skipblock = True
                if d and d not in names:
                    logger.warning(
                        "DISPLACEMENTS file: found start of domain "
                        "block with invalid domain name {}. Block will be "
                        "skipped.".format(d))
                    d = ""
                    skipblock = True
            elif skipblock:
                continue
            else:
                if not d:
                    logger.warning(
                        "DISPLACEMENTS file: found line outside "
                        "of a domain block, line will be skipped: "+line)
                else:
                    dlines[d].append(line)
        for dp in rp.domainParams:
            dp.rpars.disp_blocks.append((dlines[dp.name], blockname))
    return


def readDISPLACEMENTS_block(rp, sl, dispblock, only_mode=""):
    """
    Reads a block from the DISPLACEMENTS file and adds the information to
    all atoms in the slab.

    Parameters
    ----------
    rp : Rparams
        The run parameters object.
    sl : Slab
        The slab object containing atom and site information, to which the
        displacements will be written.
    dispblock : tuple (lines, name)
        The information in the DISPLACEMENTS block to be interpreted, as read
        be readDISPLACEMENTS.
    only_mode : str, optional
        If set to 'vib', 'geo' or 'occ', all other displacements
        will be skipped. This is meant to be used for error calculations, where
        only one-dimensional deltas are desired.

    Raises
    ------
    RuntimeError
        Raised for invalid or inconsistent input.
    ValueError
        Raised if displacement lists have inconsistent length, or if
        site occupation exceeds 1.

    Returns
    -------
    deltas_required : bool
        True if any displacements were assigned, False otherwise

    """
    deltas_required = False
    abst = sl.ab_cell.T
    # if the unit cell gets modified by SYM_DELTA, restore it afterwards
    uCellState = sl.ucell_mod
    (lines, name) = dispblock
    mode = 0  # 0: not reading; 1: reading geo, 2: reading vib, 3: reading occ
    regex = False   # read regular expressions as-is or not
    constraints = []  # collect constraints here first, interpret after loop
    for line in lines:
        if "!" in line:
            line = line.split("!")[0]
        line = line.strip()
        if '=' in line:
            if line[0] == '=':
                llist = line[1:].split()
                if llist[0][0].lower() == 'g':
                    mode = 1  # geometry
                    if only_mode and only_mode != "geo":
                        mode = 0
                elif llist[0][0].lower() == 'v':
                    mode = 2  # vibration
                    if only_mode and only_mode != "vib":
                        mode = 0
                elif llist[0][0].lower() == 'o':
                    mode = 3  # occupation
                    if only_mode and only_mode != "occ":
                        mode = 0
                elif llist[0][0].lower() == 'c':
                    mode = 4  # constraint
                elif llist[0][0].lower() == 'r':
                    regex = True
                    if len(llist) >= 2:
                        if llist[1].lower() == 'off':
                            regex = False
                else:
                    logger.warning("DISPLACEMENTS: Found line starting with "
                                   "'=', but didn't recognize flag.")
                    rp.setHaltingLevel(1)
                continue
            elif mode == 0:
                continue
            else:
                pside = line.split('=')[0].strip()
                if pside:
                    # get rid of spaces and check the leftmost entry.
                    plist = pside.split()
                    if plist:
                        param = plist[0]
                        if param[0] == '!':
                            continue
                else:
                    continue
        else:
            continue
        value = line.split('=')[1].strip()
        try:
            llist = value.split()
        except IndexError:
            logger.warning("DISPLACEMENTS file: " + param + " appears to "
                           "have no value")
            rp.setHaltingLevel(1)
            continue
        if not llist:
            logger.warning("DISPLACEMENTS file: " + param + " appears to "
                           "have no value")
            rp.setHaltingLevel(1)
            continue
        if param.lower() == "sym_delta":
            s = llist[0].lower()
            grouplist = [
                "p1", "p2", "pm", "pg", "cm", "rcm", "pmm", "pmg", "pgg",
                "cmm", "rcmm", "p4", "p4m", "p4g", "p3", "p3m1", "p31m",
                "p6", "p6m"]  # TODO: use gui or put as const in leedbase
            targetsym = ""
            if s[0] == "t":
                # True - go to highest symmetry
                targetsym = "found"
            elif s[0] == "f":
                # False - set p1
                targetsym = "p1"
            elif s in grouplist:
                if s not in ["cm", "pmg"]:
                    targetsym = s
                else:
                    logger.warning('DISPLACEMENTS file: SYM_DELTA: For '
                                   'group '+s+', direction needs to be '
                                   'specified. Input will be ignored.')
                    rp.setHaltingLevel(2)
            elif s[:2] in ["pm", "pg", "cm"] or s[:3] == "pmg":
                # regex to read
                rgx = re.compile(
                    r'\s*(?P<group>(pm|pg|cm|rcm|pmg))\s*'
                    r'\[\s*(?P<i1>[-012]+)\s+(?P<i2>[-012]+)\s*\]')
                m = rgx.match(value.lower())
                if not m:
                    logger.warning(
                        'DISPLACEMENTS file: SYM_DELTA: Could not parse given '
                        'value. Input will be ignored.')
                    rp.setHaltingLevel(2)
                else:
                    i1 = i2 = -2
                    group = m.group('group')
                    try:
                        i1 = int(m.group('i1'))
                        i2 = int(m.group('i2'))
                    except ValueError:
                        logger.warning('DISPLACEMENTS file: SYM_DELTA: '
                                       'Could not parse given value. '
                                       'Input will be ignored.')
                        rp.setHaltingLevel(1)
                    if (group in ["pm", "pg", "cm", "rcm", "pmg"]
                            and i1 in range(-1, 3) and i2 in range(-1, 3)):
                        targetsym = m.group(0)
                    else:
                        logger.warning('DISPLACEMENTS file: SYM_DELTA: '
                                       'Could not parse given value. '
                                       'Input will be ignored.')
                        rp.setHaltingLevel(2)
            else:
                logger.warning('DISPLACEMENTS file: SYM_DELTA: Could not '
                               'parse given value. Input will be ignored.')
                rp.setHaltingLevel(2)
            if targetsym != "":
                sl.revert_unit_cell(uCellState)
                setSymmetry(sl, rp, targetsym)
                enforceSymmetry(sl, rp, movement=False, rotcell=False)
            continue
        # if we are still here, then it's a "normal" line
        ats = []
        if mode == 1:
            rgx = re.compile(
                r'\s*(?P<nums>(([-:0-9]|((l|L)\([-:0-9\s]+\))'
                r')\s*)*)(?P<dir>(z|((ab|xy)\s*)?\[\s*[-0-9]+\s+[-0-9]+\s*\]|'
                r'(azi|r)\(?((ab|xy)\s*)?\[\s*[-0-9\.]+\s+[-0-9\.]+\s*\]\)?)?'
                r'(\s*offset)?)')
        else:
            rgx = re.compile(r'(?P<nums>(([-:0-9]|((l|L)\([-:0-9\s]+\)))\s*)*)'
                             r'(\s*offset)?')
        if mode == 4:
            if param.lower() == "geo":
                ctype = 1
            elif param.lower() == "vib":
                ctype = 2
            elif param.lower() == "occ":
                ctype = 3
            else:
                logger.warning("DISPLACEMENTS file: CONSTRAINT flag " + param +
                               " not recognized. Input will be ignored.")
                rp.setHaltingLevel(2)
                continue
            pspec = pside.split(maxsplit=1)[1]
        else:
            pspec = pside
        _break = True
        for s in pspec.split(","):
            s = s.strip()
            sname = s.split()[0]
            spec = s.split(sname, 1)[1].strip().lower()
            m = rgx.match(spec)
            if not m or (mode == 1 and not m.group("dir")):
                logger.warning('DISPLACEMENTS file: could not parse as '
                               'numbers or direction, skipping line: '+pside)
                rp.setHaltingLevel(2)
                break
            ats.append((sname, m.group("nums")))
            if mode == 1:
                dr = m.group("dir")
                if not dr:
                    logger.warning('DISPLACEMENTS file: could not parse '
                                   'direction, skipping line: '+pside)
                    rp.setHaltingLevel(2)
                    break
                # quick consistency check on direction:
                if "[" in dr:
                    if "azi" not in dr and "r" not in dr:
                        try:
                            dirints = (
                                int(dr.split("[")[1].strip().split()[0]),
                                int(dr.split("]")[0].strip().split()[-1]))
                        except ValueError:
                            logger.warning(
                                'DISPLACEMENTS file: could not parse '
                                'direction "'+dr+'", skipping line: '+pside)
                            rp.setHaltingLevel(2)
                            break
                    else:
                        try:
                            dirfloats = (
                                float(dr.split("[")[1].strip().split()[0]),
                                float(dr.split("]")[0].strip().split()[-1]))
                        except ValueError:
                            logger.warning(
                                'DISPLACEMENTS file: could not parse '
                                'direction "'+dr+'", skipping line: '+pside)
                            rp.setHaltingLevel(2)
                            break
        else:
            _break = False
        if _break:
            rp.setHaltingLevel(2)
            continue
        # get list of atoms to manipulate
        targetAtEls = []
        _break = True
        for (sname, nums) in ats:
            numlist = nums.split()
            if numlist:
                # get proper numlist as integers
                # first need to recombine expressions like "L(2 3)"
                nl = []
                partstring = ""
                bracketcount = 0
                for s in numlist:
                    partstring += s
                    bracketcount += s.count("(") - s.count(")")
                    if bracketcount == 0:
                        nl.append(partstring)
                        partstring = ""
                    else:
                        partstring += " "
                if bracketcount != 0:
                    logger.warning("DISPLACEMENTS file: Bracket mismatch in "
                                   "line " + pside)
                    if len(partstring) > 0:
                        nl.append(partstring)
                numlist = []
                for s in nl:
                    if "l(" not in s:
                        try:
                            numlist.extend(read_int_range(s))
                        except ValueError:
                            logger.warning(
                                'DISPLACEMENTS file: could not parse '
                                f'integer range, skipping line: {pside}'
                                )
                            rp.setHaltingLevel(2)
                            break
                    else:
                        m = re.match(r'l\((?P<laynum>[-:0-9\s]+)\)', s)
                        if not m:
                            logger.warning(
                                'DISPLACEMENTS file: could not parse layer '
                                f'expression, skipping line: {pside}'
                                )
                            rp.setHaltingLevel(2)
                            break
                        try:
                            il = read_int_range(m['laynum'])
                        except ValueError:
                            logger.warning(
                                'DISPLACEMENTS file: could not parse layer '
                                f'expression, skipping line: {pside}'
                                )
                            rp.setHaltingLevel(2)
                            break
                        for ln in il:
                            if ln > sl.n_layers:
                                logger.warning(
                                    'DISPLACEMENTS file: layer number out of '
                                    'bounds, skipping line: '+pside)
                                rp.setHaltingLevel(2)
                                break
                            numlist.extend([at.num for at in sl
                                            if at.layer == sl.layers[ln-1]])
                else:  # loop finished without beak
                    _break = False
                if _break:
                    break
                _break = True
            # interpret label:
            targetsites = []
            # first try with POSCAR element:
            for site in sl.sitelist:
                s = sname
                if "_" not in s:
                    if regex:
                        s += ".*"
                    else:
                        s += "*"
                if not regex:
                    s = re.escape(s)   # double-slash non-literal characters
                    # if regular expressions are not enabled, we want to
                    #  still interpret * as "any number of any characters"
                    s = s.replace('\\*', '.*')
                m = re.match(s, site.label)
                if m:
                    # if the length of the matched text == the site label,
                    #  it's a real match
                    if m.end(0) == len(site.label):
                        targetsites.append(site)
            if len(targetsites) > 0:
                targetel = ""   # use ALL elements of site
            elif mode == 3:
                # if we're reading occupation, only POSCAR elements are
                #   allowed -> error
                logger.warning(
                    'DISPLACEMENTS file: could not parse '+sname+' as POSCAR '
                    'element or site label, skipping line.')
                rp.setHaltingLevel(2)
                break
            else:
                # check ELEMENT_MIX elements:
                el = sname.split("_")[0].capitalize()
                s = el
                if not regex:
                    s = re.escape(s)   # double-slash
                    s = s.replace('\\*', '.*')
                cels = [cel for cel in sl.chemelem if re.match(s, cel)]
                if el not in sl.elements and len(cels) > 0:
                    pel = ""
                    for k in rp.ELEMENT_MIX:
                        if any(e in rp.ELEMENT_MIX[k] for e in cels):
                            pel = k
                            targetel = [e for e in cels
                                        if e in rp.ELEMENT_MIX[k]][0]
                            break
                    s = sname.replace(sname.split("_")[0], pel, 1)
                    if s in sl.elements:  # is element, all sites
                        targetsites.extend([site for site in sl.sitelist if
                                            site.el == s])
                    else:
                        for site in sl.sitelist:
                            if not regex:
                                s = re.escape(s)   # double-slash
                                s = s.replace('\\*', '.*')
                            m = re.match(s, site.label)
                            if m:
                                if m.end(0) == len(site.label):
                                    targetsites.append(site)
                else:
                    logger.warning(
                        'DISPLACEMENTS file: could not parse '+sname+' as '
                        'element or site label, skipping line.')
                    rp.setHaltingLevel(2)
                    break
                if len(targetsites) == 0:
                    logger.warning(
                        'DISPLACEMENTS file: could not parse '+sname+' as '
                        'element or site label, skipping line.')
                    rp.setHaltingLevel(2)
                    break
            for at in sl:
                if ((at.num in numlist or len(numlist) == 0)
                        and at.site in targetsites and not at.is_bulk):
                    targetAtEls.append((at, targetel))
        else:  # loop finished without break
            _break = False
        if _break:
            continue  # error message has already happened
        if len(targetAtEls) == 0:
            logger.warning('DISPLACEMENTS file: no atoms found for '
                           'line, line will be skipped: '+pside)
            rp.setHaltingLevel(1)
            continue
        if mode == 2 or (mode == 1 and "offset" not in dr):
            # geometric or vibration displacement, get range:
            try:
                fl = [float(s) for s in llist]
            except ValueError:
                if ("offset" in pside and llist[0].lower() in ["clear",
                                                               "original"]):
                    fl = [0.]
                elif "offset" in pside and llist[0].lower().startswith("cont"):
                    continue
                    # instruction to use value from previous search
                    # -> default, ignore
                else:
                    fl = []
                    logger.warning(
                        'DISPLACEMENTS file: could not parse ' + value +
                        ' as list of floats, skipping line.')
                    rp.setHaltingLevel(1)
                    continue
            if len(fl) < 3:
                if mode == 2 and len(fl) == 1:  # interpret as static offset
                    fl = [fl[0], fl[0], 1]
                else:
                    logger.warning('DISPLACEMENTS file: too few values '
                                   'found, skipping line: '+pside)
                    rp.setHaltingLevel(1)
                    continue
            #  OLD: force uneven number of steps
            # steps = abs(int(round((fl[1]-fl[0]) / fl[2])))+1
            # mid = (fl[1]+fl[0]) / 2
            # if steps % 2 == 0:      # even number of steps, extend range
            #     steps += 1
            # drange = np.arange(mid - ((steps-1)/2*fl[2]),
            #                    mid + ((steps-1)/2*fl[2])+1e-6, fl[2])
            if fl[0] < fl[1]:
                drange = np.arange(fl[0], fl[1] + 1e-6, abs(fl[2]))
            else:
                drange = np.arange(fl[0], fl[1] - 1e-6, -abs(fl[2]))
            if min([abs(v) for v in drange]) > 5e-5:
                logger.warning(
                    "DISPLACEMENTS: A range does not contain zero. This means "
                    "that keeping the input configuration unchanged will not "
                    "be possible during the search: " + value)
                rp.setHaltingLevel(1)
        if mode == 1 and "offset" in dr:
            # geo offset, get value
            try:
                offval = float(llist[0])
            except (ValueError, IndexError):
                if llist[0].lower() in ["clear", "original"]:
                    offval = 0.
                elif llist[0].lower().startswith("cont"):
                    continue
                else:
                    logger.warning('DISPLACEMENTS file: could not parse '
                                   + value + ' as float, skipping line.')
                    rp.setHaltingLevel(1)
                    continue
            drange = [offval]
        if mode == 1:
            if dr.strip() == "offset":
                # No direction
                if offval != 0.:
                    logger.warning('DISPLACEMENTS file: cannot assign '
                                   'geo offset : no direction given: '+pside)
                    rp.setHaltingLevel(1)
                    continue
                else:
                    for (at, targetel) in targetAtEls:
                        at.clearOffset(1, targetel)
                    continue
            # geometric displacement, get direction vector:
            drvec = np.array([0., 0., 0.])
            if "azi" not in dr and "r" not in dr:
                if "z" in dr:
                    drvec[2] = -1.  # invert dir: positive = away from bulk
                elif "xy" in dr:
                    v = np.array(dirints)
                    v = v / np.linalg.norm(v)
                    drvec[:2] = v
                else:
                    v = dirints[0]*abst[0] + dirints[1]*abst[1]
                    v = v / np.linalg.norm(v)
                    drvec[:2] = v
                disprange = [f*drvec for f in drange]
                for (at, targetel) in targetAtEls:
                    # check whether displacement is allowed
                    if "z" in dr:
                        allowed = True
                    elif type(at.freedir) == int:
                        if at.freedir == 0:
                            allowed = False
                        else:
                            allowed = True
                    else:
                        freev = at.freedir[0]*abst[0] + at.freedir[1]*abst[1]
                        freev = freev / np.linalg.norm(freev)
                        if abs(abs(np.dot(drvec[:2], freev))-1) < 1e-5:
                            allowed = True
                        else:
                            allowed = False
                    if allowed:
                        if "offset" not in dr:
                            at.assignDisp(mode, disprange, targetel,
                                          disp_label=dr,
                                          disp_lin_steps=drange)
                        else:
                            at.clearOffset(1, targetel)
                            at.assignDisp(4, disprange, targetel)
                        deltas_required = True
                    else:
                        logger.warning(
                            f'In-plane displacement assignment for {at} is '
                            'forbidden by symmetry and will be skipped. See '
                            '"FreeDir" in POSCAR file. To apply this '
                            'displacement, use either the SYMMETRY_FIX '
                            'parameter to lower the symmetry, or use '
                            'SYM_DELTA in the DISPLACEMENTS file to allow '
                            'symmetry breaking for this atom.'
                            )
            else:
                if "xy" in dr:
                    c = np.array(dirfloats)
                else:
                    c = dirfloats[0]*abst[0] + dirfloats[1]*abst[1]
                for (at, targetel) in targetAtEls:
                    if "r" in dr:
                        v = at.cartpos[:2] - c
                        v = v / np.linalg.norm(v)
                        allowed = True
                        if type(at.freedir) == int and at.freedir == 0:
                            allowed = False
                        elif type(at.freedir) != int:
                            freev = (at.freedir[0]*abst[0] +
                                     at.freedir[1]*abst[1])
                            freev = freev / np.linalg.norm(freev)
                            if abs(abs(np.dot(v, freev))-1) > 1e-5:
                                allowed = False
                        if allowed:
                            drvec = np.array([0., 0., 0.])
                            drvec[:2] = v
                            disprange = [f*drvec for f in drange]
                            if "offset" not in dr:
                                at.assignDisp(mode, disprange, targetel,
                                              disp_label=dr,
                                              disp_lin_steps=drange)
                            else:
                                at.assignDisp(4, disprange, targetel)
                            deltas_required = True
                    else:  # azi
                        # allowed only for completely free atoms
                        if type(at.freedir) == int and at.freedir == 1:
                            r = at.cartpos[:2] - c
                            disprange = []
                            for d in drange:
                                # disprange based on rotating vector around c
                                ar = d / np.linalg.norm(r)  # rotation angle
                                tm = np.array([[np.cos(ar)-1, -np.sin(ar)],
                                               [np.sin(ar), np.cos(ar)-1]])
                                disprange.append(np.append(np.dot(tm, r), 0.))
                            if "offset" not in dr:
                                at.assignDisp(mode, disprange, targetel,
                                              disp_label=dr,
                                              disp_lin_steps=drange)
                            else:
                                at.clearOffset(1, targetel)
                                at.assignDisp(4, disprange, targetel)
                            deltas_required = True
                        else:
                            logger.warning(
                                'In-plane azimuthal displacement assignment '
                                f'for {at} is forbidden by symmetry and will '
                                'be skipped. See "FreeDir" in POSCAR file. '
                                'To apply this displacement, use either the '
                                'SYMMETRY_FIX parameter to lower the symmetry,'
                                ' or use SYM_DELTA in the DISPLACEMENTS file '
                                'to allow symmetry breaking for this atom.'
                                )
        elif mode == 2:
            # vibration displacement, apply:
            for (at, targetel) in targetAtEls:
                if "offset" in pside:
                    at.clearOffset(mode, targetel)
                if not (len(drange) == 1 and drange[0] == 0.):
                    at.assignDisp(mode, drange, targetel,
                                  disp_label="N/A",
                                  disp_lin_steps=drange)
                deltas_required = True
        elif mode == 3:
            # occupations, get ranges:
            sublists = splitSublists(llist, ',')
            maxsteps = 0
            _break = True
            for subl in sublists:
                # get range
                try:
                    fl = [float(s) for s in subl[-3:]]
                except (ValueError, IndexError):
                    fl = []
                    if len(subl) == 2:
                        try:
                            fl = [float(subl[-1])]
                        except Exception:
                            pass
                    if len(fl) == 1 or (len(fl) == 3 and fl[2] < 5e-5):
                        fl = [fl[0], fl[0], 1]
                    if len(fl) < 3:
                        logger.warning(
                            'DISPLACEMENTS file: could not parse ' + value
                            + ' as list of floats, skipping line.')
                        rp.setHaltingLevel(1)
                        break
                if len(fl) < 3 or (len(fl) == 3 and fl[2] < 5e-5):
                    if len(fl) == 1 or (len(fl) == 3 and fl[2] < 5e-5):
                        # interpret as static offset
                        fl = [fl[0], fl[0], 1]
                    else:
                        logger.warning('DISPLACEMENTS file: too few values '
                                       'found, skipping line: '+pside)
                        rp.setHaltingLevel(1)
                        break
                steps = abs(int(round((fl[1]-fl[0]) / fl[2])))+1
                if steps > maxsteps:
                    maxsteps = steps
            else:
                _break = False
            if _break:
                continue
            for subl in sublists:
                try:
                    fl = [float(s) for s in subl[-3:]]
                except (ValueError, IndexError):
                    f = float(subl[-1])  # checked before, will work
                    fl = [f, f, 1]
                if fl[2] < 5e-5:
                    fl = [fl[0], fl[0], 1]
                steps = abs(int(round((fl[1]-fl[0]) / fl[2])))+1
                if steps > maxsteps:
                    steps = maxsteps
                    logger.warning(
                        'DISPLACEMENTS file: inconsistent step '
                        'numbers for occupancies of '+pside+'. Decreasing '
                        'spacing for '+subl[0]+' to make number of steps '
                        'equal.')
                if steps == 1 and maxsteps > 1:
                    drange = np.array([fl[0]] * maxsteps)
                else:
                    mid = (fl[1]+fl[0]) / 2
                    drange = np.arange(mid-((steps-1)/2*fl[2]),
                                       mid+((steps-1)/2*fl[2])+1e-5, fl[2])
                    if fl[1] < fl[0]:
                        drange = np.flip(drange)
                if subl[0].capitalize() in sl.chemelem:
                    targetel = subl[0]
                else:
                    logger.error('DISPLACEMENTS file: '+subl[0]+' not '         # TODO: This should probably raise an exception
                                   'found in element list. No assignment '
                                   'will be made.')
                    rp.setHaltingLevel(2)
                    continue
                for (at, _) in targetAtEls:
                    if targetel in at.disp_occ:
                        at.assignDisp(mode, drange, targetel,
                                      disp_label="N/A",
                                      disp_lin_steps=None)  # for occupation, we figure it out in r_error
                        deltas_required = True
                    else:
                        logger.warning(
                            'DISPLACEMENTS file: trying to address {} with '
                            'wrong element. Atom will be skipped.'.format(at))
                        rp.setHaltingLevel(1)
                # else:
                #     # !!! currently unused, not in wiki - delete?             # TODO: Not true! It is in the Syntax example in the documentation of OCC_DELTA, and described in detail there!
                #     elparts = "".join(subl[:-3]).split("+")
                #     elweights = []
                #     elsum = 0.0
                #     rgx = re.compile(r'\s*(?P<number>[0-9\.]+)\s*'
                #                      +r'(?P<name>[a-zA-Z]+)')
                #     for part in elparts:
                #         m = rgx.match(part)
                #         if not m:
                #             logger.warning('DISPLACEMENTS file: could '
                #                 'not parse '+part+' as number and '
                #                 'element. No assignment will be made.')
                #             rp.setHaltingLevel(1)
                #             break
                #         try:
                #             f = float(m.group("number"))
                #         except ValueError:
                #             logger.warning('DISPLACEMENTS file: '
                #              'could not parse '+m.group("number")
                #              +' as float. No assignment will be made.')
                #             rp.setHaltingLevel(1)
                #             break
                #         el = m.group("name").capitalize()
                #         if not el in sl.chemelem:
                #             logger.warning('DISPLACEMENTS file: '+el
                #                 +' not found in element list. No '
                #                 'assignment will be made.')
                #             rp.setHaltingLevel(1)
                #             break
                #         else:
                #             elweights.append((el, f))
                #             elsum += f
                #     else:  # loop finished without break
                #         for (el, f) in elweights:
                #             ndr = [v*f/elsum for v in drange]
                #             for (at, targetel) in targetAtEls:
                #                 if el in at.disp_occ:
                #                     at.assignDisp(mode, ndr, el)
                #                 else:
                #                     logger.warning('DISPLACEMENTS file: '
                #                         'trying to address atom number '
                #                         +str(at.num)+' with wrong element. '
                #                         'Atom will be skipped.')
                #                     rp.setHaltingLevel(1)
        elif mode == 4:
            constraints.append((targetAtEls, ctype, value))
    sl.revert_unit_cell(uCellState)  # if modified by SYM_DELTA, go back
    # now read constraints
    for (targetAtEls, ctype, value) in constraints:
        if value.lower().startswith("link"):
            linkto = targetAtEls[0]
            for (at, el) in targetAtEls:  # [1:]:
                at.assignConstraint(ctype, targetel=el, linkAtEl=linkto)
        elif value.lower().startswith("ind("):
            try:
                ind = int(value.lower().split("ind(")[1].split(")")[0])
            except ValueError:
                logger.warning("DISPLACEMENTS file: Could not convert "+value
                               + " to index, skipping line.")
                continue
            for (at, el) in targetAtEls:
                at.assignConstraint(ctype, targetel=el, index=ind)
        else:
            if ctype == 3:
                occEl = value.split()[0]
                value = value.split(maxsplit=1)[-1]
            try:
                f = float(value.split()[0])
            except ValueError:
                logger.warning("DISPLACEMENTS file: Could not convert "+value
                               + " to float, skipping line.")
                continue
            for (at, el) in targetAtEls:
                if ctype == 3:
                    el = occEl
                at.assignConstraint(ctype, targetel=el, value=f)
    # Reading done, consistency check...
    o = ""
    if name:
        o = " "+name
    logger.debug("DISPLACEMENTS block{} was read successfully".format(o))
    for at in sl:
        if len(at.disp_occ) > 5:
            logger.error(
                'DISPLACEMENTS file: '+str(at)+' has occupations defined for '
                'more than 5 elements, which is not supported by TensErLEED.')
            raise RuntimeError("Too many elements")
        occlists = []
        error = False
        for k in at.disp_occ:
            occlists.append(at.disp_occ[k])
        for i in range(0, len(occlists[0])):
            totalocc = 0.
            for ol in occlists:
                if len(ol) <= i:
                    error = True
                    break
                else:
                    totalocc += ol[i]
            if error:
                logger.error('DISPLACEMENTS file: Lengths of occupation '
                             'lists per element differ for '+str(at))
                raise ValueError("Occupation lists per element must have same "
                                 "length.")
            if totalocc > 1 + 1e-4:    # some tolerance
                logger.error('DISPLACEMENTS file: Occupations for ' + str(at)
                             + ' sum to more than 1 for step '+str(i+1)+'.')
                raise ValueError("Occupations must sum to 1 or less.")
            elif totalocc < 1 - 1e-4 and len(at.disp_occ) > 4:
                logger.error(
                    'DISPLACEMENTS file: Occupations for ' + str(at)
                    + ' sum to less than 1 for step '+str(i+1)+', but a '
                    'vacancy cannot be added since there are already '
                    '5 elements.')
                raise ValueError("Occupations must sum to exactly 1 if five "
                                 "elements are already defined.")
    return deltas_required
