"""Data and functions related to the periodic table and element properties.

The functionality in this module used to be part of calc.lib.leedbase.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-06-07'
__license__ = 'GPLv3+'

PERIODIC_TABLE = (
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na',
    'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr',
    'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
    'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po',
    'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
    'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs',
    'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
    )


COVALENT_RADIUS = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96,
    "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05,
    "Cl": 1.02, "Ar": 1.06, "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60,
    "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24,
    "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20,
    "Br": 1.20, "Kr": 1.16, "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75,
    "Nb": 1.64, "Mo": 1.54, "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39,
    "Ag": 1.45, "Cd": 1.44, "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38,
    "I": 1.39, "Xe": 1.40, "Cs": 2.44, "Ba": 2.15, "La": 2.07, "Ce": 2.04,
    "Pr": 2.03, "Nd": 2.01, "Pm": 1.99, "Sm": 1.98, "Eu": 1.98, "Gd": 1.96,
    "Tb": 1.94, "Dy": 1.92, "Ho": 1.92, "Er": 1.89, "Tm": 1.90, "Yb": 1.87,
    "Lu": 1.87, "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44,
    "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "Tl": 1.45, "Pb": 1.46,
    "Bi": 1.48, "Po": 1.40, "At": 1.50, "Rn": 1.50, "Fr": 2.60, "Ra": 2.21,
    "Ac": 2.15, "Th": 2.06, "Pa": 2.00, "U": 1.96, "Np": 1.90, "Pu": 1.87,
    "Am": 1.80, "Cm": 1.69}
# from Cordero et al., 2008 (DOI: 10.1039/B801115J)


ATOMIC_MASS = {
    "H": 1.00797, "He": 4.00260, "Li": 6.941, "Be": 9.01218,
    "B": 10.81, "C": 12.011, "N": 14.0067, "O": 15.9994, "F": 18.998403,
    "Ne": 20.179, "Na": 22.98977, "Mg": 24.305, "Al": 26.98154, "Si": 28.0855,
    "P": 30.97376, "S": 32.06, "Cl": 35.453, "K": 39.0983, "Ar": 39.948,
    "Ca": 40.08, "Sc": 44.9559, "Ti": 47.90, "V": 50.9415, "Cr": 51.996,
    "Mn": 54.9380, "Fe": 55.847, "Ni": 58.70, "Co": 58.9332, "Cu": 63.546,
    "Zn": 65.38, "Ga": 69.72, "Ge": 72.59, "As": 74.9216, "Se": 78.96,
    "Br": 79.904, "Kr": 83.80, "Rb": 85.4678, "Sr": 87.62, "Y": 88.9059,
    "Zr": 91.22, "Nb": 92.9064, "Mo": 95.94, "Tc": 98, "Ru": 101.07,
    "Rh": 102.9055, "Pd": 106.4, "Ag": 107.868, "Cd": 112.41, "In": 114.82,
    "Sn": 118.69, "Sb": 121.75, "I": 126.9045, "Te": 127.60, "Xe": 131.30,
    "Cs": 132.9054, "Ba": 137.33, "La": 138.9055, "Ce": 140.12,
    "Pr": 140.9077, "Nd": 144.24, "Pm": 145, "Sm": 150.4, "Eu": 151.96,
    "Gd": 157.25, "Tb": 158.9254, "Dy": 162.50, "Ho": 164.9304, "Er": 167.26,
    "Tm": 168.9342, "Yb": 173.04, "Lu": 174.967, "Hf": 178.49, "Ta": 180.9479,
    "W": 183.85, "Re": 186.207, "Os": 190.2, "Ir": 192.22, "Pt": 195.09,
    "Au": 196.9665, "Hg": 200.59, "Tl": 204.37, "Pb": 207.2, "Bi": 208.9804,
    "Po": 209, "At": 210, "Rn": 222, "Fr": 223, "Ra": 226.0254, "Ac": 227.0278,
    "Pa": 231.0359, "Th": 232.0381, "Np": 237.0482, "U": 238.029}


def get_atomic_number(element):
    """Return atomic number for a given element symbol.

    Parameters
    ----------
    element : str
        String of the element symbol (e.g. 'Fe').

    Returns
    -------
    int
        Atomic number of the element (e.g. 26 for 'Fe').

    Raises
    ------
    ValueError
        If element is not a valid chemical element.
    """
    try:
        # Offset by one because of Python indexing
        _element = element.lower().capitalize()
        return PERIODIC_TABLE.index(_element) + 1
    except ValueError:
        raise ValueError(
            f"Unknown chemical element {element}"
            ) from None


def get_element_symbol(atomic_number):
    """Returns element symbol for given atomic number.

    Parameters
    ----------
    atomic_number : int
        Atomic number of the element (e.g. 26 for 'Fe').

    Returns
    -------
    str
        String of the element symbol (e.g. 'Fe').

    Raises
    ------
    ValueError
        If atomic_number is not a valid Z for a chemical element.
    """
    try:
        # Offset by one because of Python indexing
        return PERIODIC_TABLE[atomic_number - 1]
    except (IndexError, TypeError):
        raise ValueError(
            f"Invalid atomic number {atomic_number}."
            ) from None
