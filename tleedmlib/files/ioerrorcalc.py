# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 17:22:20 2021

@author: Florian Kraushofer

Functions for reading and writing files relevant to the error calculation
"""

import numpy as np
import logging

logger = logging.getLogger("tleedm.files.ioerrorcalc")


def write_errors_csv(errors, filename="Errors.csv", sep=";"):
    """
    Writes errors from the error calculation into a CSV file

    Parameters
    ----------
    errors : R_Error
        Data structure from sections.errorcalc containing the information
    filename : str, optional
        Name of the csv file to write. The default is "Errors.csv".
    sep : str, optional
        The separator to use. The default is ";".

    Returns
    -------
    None.

    """
    # contents of the columns; start with only titles:
    columns = {"at": ["Atoms"],
               "mode": ["Mode"],
               "dir": ["Direction (1st atom)"],
               "disp": ["Displacement [A]"],
               "rfac": ["R"]}
    for err in errors:
        ats = ", ".join(str(at.oriN) for at in err.atoms)
        if isinstance(err.displacements[0], np.ndarray):
            dirvec = err.displacements[-1] * np.array([1, 1, -1])
            dirvec = dirvec / np.linalg.norm(dirvec)
            direction = ("[" + ", ".join([str(round(f, 4))
                                          for f in dirvec]) + "]")
        for i in range(0, len(err.rfacs)):
            columns["at"].append(ats)
            columns["mode"].append(err.mode)
            if isinstance(err.displacements[i], np.ndarray):
                columns["dir"].append(direction)
                disp = (np.linalg.norm(err.displacements[i])
                        * np.sign(np.dot(dirvec * np.array([1, 1, -1]),
                                         err.displacements[i])))
            else:
                columns["dir"].append("")
                disp = err.displacements[i]
            columns["disp"].append("{:.4f}".format(disp))
            columns["rfac"].append("{:.4f}".format(err.rfacs[i]))

    if all(s == "" for s in columns["dir"]):
        del columns["dir"]

    widths = {}
    for k in columns:
        widths[k] = max(len(s) for s in columns[k]) + 1

    output = ""
    for i in range(len(next(iter(columns.values())))):
        for k in columns:
            output += columns[k][i].rjust(widths[k]) + sep
        output = output[:-1] + "\n"

    try:
        with open(filename, "w") as wf:
            wf.write(output)
    except Exception:
        logger.warning("Failed to write "+filename)
    return
