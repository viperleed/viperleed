# -*- coding: utf-8 -*-
"""Moule _domain_params of viperleed.calc.classes.rparams.

Defines the DomainParameters class. Contains of information useful
when a calculation with multiple structural domains is carried out.
This module was originally (2019-06-13) part of the rparams.py module,
refactored by Michele Riva in Oct 2023.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__created__ = '2023-10-23'

import logging
from pathlib import Path

_LOGGER = logging.getLogger('tleedm.rparams')


class DomainParameters:
    """Stores workdir, slab and Rparams objects for each domain"""

    def __init__(self, workdir, homedir, name):
        self.workdir = Path(workdir)  # path to sub-directory for domain calculation
        self.homedir = Path(homedir)  # path to main tleedm working directory
        self.name = name        # domain name as defined by user
        self.sl = None
        self.rp = None
        self.refcalcRequired = False
        self.tensorDir = None
