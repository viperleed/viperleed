"""Moule domain_params of viperleed.calc.classes.rparams.

Defines the DomainParameters class. Contains of information useful
when a calculation with multiple structural domains is carried out.
This module was originally (2019-06-13) part of the rparams.py module,
refactored by Michele Riva in Oct 2023.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = '2019-2024 ViPErLEED team'
__created__ = '2023-10-23'
__license__ = 'GPLv3+'

from pathlib import Path


class DomainParameters:
    """Information about one structural domain."""

    def __init__(self, workdir, name):
        """Initialize instance.

        Parameters
        ----------
        workdir : pathlike
            Path to the sub-directory of the main work directory
            where calculations for this domain are performed.
        name : str
            The name of this domain (e.g., the one defined by
            the user via the DOMAIN parameter).
        """
        self.workdir = Path(workdir).resolve()
        self.name = name
        self.sl = None
        self.rp = None
        self.refcalcRequired = False
        self.tensorDir = None

    def __str__(self):
        """Return a string representation for this domain."""
        return f'domain {self.name}'
