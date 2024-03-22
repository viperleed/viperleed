"""Module beam of viperleed.calc.classes.

Defines the Beam class, a container of information concerning a
diffraction beam.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2019-12-03'
__license__ = 'GPLv3+'

import math

import numpy as np
from quicktions import Fraction


class Beam:
    """Class storing properties and data of a diffraction beam.

    Stores h,k, and can store intensity over energy.

    Attributes
    ----------
    hkfrac : tuple (Fraction, Fraction)
        The h and k coordinates of the beam as Fractions
    hk : tuple (float, float)
        The h and k coordinates of the beam as float
    intens : dict {float: float}
        contains {energy: intensity} pairs
    complex_amplitude : dict {float: complex}
        contains {energy: amplitude} pairs
    label : str
        Label for the beam
    """

    def __init__(self, hk, maxdenom=99):
        if all([isinstance(v, Fraction) for v in hk]):
            self.hkfrac = hk
            self.hk = (float(hk[0]), float(hk[1]))
        else:
            self.hk = hk
            self.hkfrac = (Fraction(hk[0]).limit_denominator(maxdenom),
                           Fraction(hk[1]).limit_denominator(maxdenom))
        self.intens = {}
        self.complex_amplitude = {}
        self.label, _ = self.getLabel()

    @property
    def energies(self):
        return np.fromiter(self.intens.keys(), float)
    
    @property
    def intensities(self):
        return np.fromiter(self.intens.values(), float)
    
    def updateIndex(self, hk, maxdenom=99):
        """Keep values but change indices"""
        if all([isinstance(v, Fraction) for v in hk]):
            self.hkfrac = hk
            self.hk = (float(hk[0]), float(hk[1]))
        else:
            self.hk = hk
            self.hkfrac = (Fraction(hk[0]).limit_denominator(maxdenom),
                           Fraction(hk[1]).limit_denominator(maxdenom))
        self.label, _ = self.getLabel()

    def getLabel(self, style="plain", lwidth=-1):
        """
        Returns a string of format '( 1/2 | -1   )', where a,b in '(a|b)'
        are justified to at least the given width. If h or k end up longer,
        the width of the other one will also be changed.

        Parameters
        ----------
        style : str, optional
            Select a different style for the label.

            plain:
                UTF8-compatible, with hyphen '-' as minus, for log
            minus:
                replace hyphen by proper minus sign, '\u2212'
            overbar:
                use latex overbars instead of minus
        lwidth: int, optional
            width of each component of the label; total width = 2*lwidth+3.
            Leave negative to determine automatically

        Returns
        -------
        out : str
            The new label.
        lwidth : int
            The width of the new label

        """
        out = ""
        i = 0
        if lwidth < 0:
            lwidth = 2
            if style == "overbar":
                lwidth = 1
        while i < 2:
            texlen = 0  # added length from tex commands
            if i == 0:
                wl = math.ceil((lwidth - 1) / 2)  # width for numerator
            v = self.hkfrac[i]
            if v.denominator == 1:
                if style != "overbar":
                    s = str(v.numerator).rjust(max(2, wl))   # use more space
                else:
                    if v.numerator < 0:
                        s = ((r"\overline{" + str(abs(v.numerator)) + "}")
                             .rjust(wl + 11))  # 11 chars for the "\overline{}"
                        texlen += 11
                    else:
                        s = str(v.numerator).rjust(wl)
            else:
                if style != "overbar":
                    s = str(v.numerator).rjust(wl)+"/"+str(v.denominator)
                else:
                    if v.numerator < 0:
                        s = ((r"\overline{" + str(abs(v.numerator)) + "/"
                              + str(v.denominator) + "}")
                             .rjust(wl + 11))  # 11 chars for the "\overline{}"
                        texlen += 11
                    else:
                        s = str(v.numerator).rjust(wl)+"/"+str(v.denominator)
            s = s.ljust(lwidth + texlen)
            if len(s) - texlen > lwidth:
                lwidth = len(s) - texlen
                i = -1     # start over
            elif i == 0:
                if style == "overbar":
                    out = "$"
                else:
                    out = ""
                out += "(" + s + "|"
            else:
                out += s + ")"
                if style == "overbar":
                    out += "$"
            i += 1
        if style == "minus":
            out = out.replace("-", "\u2212")
        return out, lwidth

    def normMax(self):
        """
        Normalizes the beam to maximum, i.e. sets the highest value to
        1.0 and rescales the others accordingly.
        """
        if len(self.intens.values()) > 0:
            m = max(self.intens.values())
            if m > 0:
                for en in self.intens:
                    self.intens[en] /= m

    def isEqual(self, beam, eps=1e-4):
        """
        Checks whether the beam is equal to another beam with a given
        tolerance. Returns True or False.
        """
        if (abs(self.hk[0] - beam.hk[0]) < eps
                and abs(self.hk[1] - beam.hk[1]) < eps):
            return True
        return False

    def isEqual_hk(self, hk, eps=1e-4):
        """
        Checks whether the beam hk is equal to a tuple hk with a given
        tolerance. Returns True or False.
        """
        hk = [float(v) for v in hk]  # in case Fractions were passed
        if (abs(self.hk[0] - hk[0]) < eps
                and abs(self.hk[1] - hk[1]) < eps):
            return True
        return False
