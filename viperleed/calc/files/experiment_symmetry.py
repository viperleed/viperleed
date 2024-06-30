"""Functions for reading and writing experiment_symmetry.ini."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-19'
__license__ = 'GPLv3+'

import logging

from viperleed.calc.classes.slab import MissingBulkSlabError


_LOGGER = logging.getLogger(__name__)


def write(slab, rpars):                                                         # TODO: shouldn't this use getLEEDdict from leedbase? Otherwise the info in rpars may be lost.
    """Write the experiment_symmetry.ini file for the ViPErLEED GUI.

    The experiment_symmetry.ini file can be used by the ViPErLEED
    pattern-simulator GUI to display the expected LEED pattern and
    export a "pattern file" to be used in the spot tracker.

    Parameters
    ----------
    slab : Slab
        The slab whose information should be written to file.
        It must have its `.bulkslab` attribute already set to
        the correct bulk slab.
    rpars : Rparams
        The current PARAMETERS.

    Raises
    ------
    MissingBulkSlabError
        If slab has no `.bulkslab` attribute.
    OSError
        If writing to file experiment_symmetry.ini fails.
    """
    output = f'eMax = {rpars.THEO_ENERGIES.max:.2f}\n'
    mstring = '[[{}, {}], [{}, {}]]'.format(*slab.ab_cell.T.ravel())
    output += f'surfBasis = {mstring}\n'
    mstring = '[[{}, {}], [{}, {}]]'.format(
        *rpars.SUPERLATTICE.round().astype(int).ravel()
        )
    output += f'superlattice = {mstring}\n'
    pgstring = slab.planegroup
    if pgstring in {'pm', 'pg', 'cm', 'rcm', 'pmg'}:                            # TODO: shouldn't we do this whenever there is a slab.orisymplane? The current implementation excludes, e.g., cmm on hex cells.
        pgstring += str(slab.orisymplane.par)
    output += f'surfGroup = {pgstring}\n'
    if slab.bulkslab is None:
        _LOGGER.error('experiment_symmetry.ini: bulk '
                      'slab has not been initialized.')
        raise MissingBulkSlabError(
            'experiment_symmetry.write called without bulk slab.'
            )
    output += f'bulkGroup = {slab.bulkslab.foundplanegroup}\n'
    output += f'bulk3Dsym = {slab.bulkslab.get_bulk_3d_str()}'
    # write output
    filename = 'experiment_symmetry.ini'
    try:  # pylint: disable=too-many-try-statements  # Two OK for open
        with open(filename, 'w', encoding='utf-8') as file:
            file.write(output)
    except OSError:
        _LOGGER.error(f'Failed to write {filename!r}')
        raise
    _LOGGER.debug(f'Wrote to {filename!r} successfully')
