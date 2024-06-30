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


def write(sl, rp, filename='experiment_symmetry.ini'):
    """Writes a experiment_symmetry.ini file that can be used by the ViPErLEED
    GUI utility to display the expected LEED pattern and show beam labelling."""
    output = f'eMax = {rp.THEO_ENERGIES.max:.2f}\n'
    mstring = '[[{}, {}], [{}, {}]]'.format(*sl.ab_cell.T.ravel())
    output += 'surfBasis = '+mstring+'\n'
    mstring = ('[[{:.0f}, {:.0f}], [{:.0f}, {:.0f}]]'
               .format(rp.SUPERLATTICE[0, 0], rp.SUPERLATTICE[0, 1],
                       rp.SUPERLATTICE[1, 0], rp.SUPERLATTICE[1, 1]))
    output += 'superlattice = {mstring}\n'
    if sl.planegroup in ['pm', 'pg', 'cm', 'rcm', 'pmg']:
        pgstring = sl.planegroup+str(sl.orisymplane.par)
    else:
        pgstring = sl.planegroup
    output += f'surfGroup = {pgstring}\n'
    if sl.bulkslab is None:
        _LOGGER.error('experiment_symmetry.ini: bulk '
                      'slab has not been initialized.')
        raise MissingBulkSlabError(
            'experiment_symmetry.write called without bulk slab.'
            )
    output += f'bulkGroup = {sl.bulkslab.foundplanegroup}\n'
    output += f'bulk3Dsym = {sl.bulkslab.get_bulk_3d_str()}'
    # write output
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except OSError:
        _LOGGER.error(f'Failed to write {filename!r}')
        raise
    _LOGGER.debug(f'Wrote to {filename!r} successfully')
