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


def write(slab, rpars, filename='experiment_symmetry.ini'):
    """Writes a experiment_symmetry.ini file that can be used by the ViPErLEED
    GUI utility to display the expected LEED pattern and show beam labelling."""
    output = f'eMax = {rpars.THEO_ENERGIES.max:.2f}\n'
    mstring = '[[{}, {}], [{}, {}]]'.format(*slab.ab_cell.T.ravel())
    output += f'surfBasis = {mstring}\n'
    mstring = ('[[{:.0f}, {:.0f}], [{:.0f}, {:.0f}]]'
               .format(rpars.SUPERLATTICE[0, 0], rpars.SUPERLATTICE[0, 1],
                       rpars.SUPERLATTICE[1, 0], rpars.SUPERLATTICE[1, 1]))
    output += f'superlattice = {mstring}\n'
    if slab.planegroup in ['pm', 'pg', 'cm', 'rcm', 'pmg']:
        pgstring = slab.planegroup+str(slab.orisymplane.par)
    else:
        pgstring = slab.planegroup
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
    try:  # pylint: disable=too-many-try-statements  # Two OK for open
        with open(filename, 'w', encoding='utf-8') as file:
            file.write(output)
    except OSError:
        _LOGGER.error(f'Failed to write {filename!r}')
        raise
    _LOGGER.debug(f'Wrote to {filename!r} successfully')
