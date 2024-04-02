"""ViPErLEED utility: Strip comments from POSCAR file."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

from viperleed.utilities.poscar.base import _PoscarStreamCLI


class StripCommentsCLI(_PoscarStreamCLI, cli_name='strip_comments'):
    """Remove ViPErLEED (or VASP) comments from a POSCAR file."""

    long_name = 'strip comments'

    def process_slab(self, slab, args):
        """Return an unchanged slab, as we only remove comments."""
        return slab


if __name__ == '__main__':
    StripCommentsCLI.run_as_script()
