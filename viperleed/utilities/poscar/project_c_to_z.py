"""ViPErLEED utility: Project c vector to z axis.

This utility takes a POSCAR slab and projects the c vector to the z axis.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

from viperleed.utilities.poscar.base import _PoscarStreamCLI


class ProjectCToZCLI(_PoscarStreamCLI, cli_name='project_c_to_z'):
    """Remove atoms between a pair of c fractions."""

    long_name = 'project c to z axis'

    def process_slab(self, slab, args):
        """Project c vector of slab to z axis."""
        slab.project_c_to_z()
        return slab


if __name__ == '__main__':
    ProjectCToZCLI.run_as_script()
