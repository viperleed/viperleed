"""Functions for reading and writing experiment_symmetry.ini."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-19'
__license__ = 'GPLv3+'

import logging

from viperleed.calc.classes.slab import MissingBulkSlabError
from viperleed.calc.lib.leedbase import getLEEDdict
from viperleed.gui.leedsim.classes.leedparser import LEEDParser


FILENAME = 'experiment_symmetry.ini'
_LOGGER = logging.getLogger(__name__)


def write(slab, rpars):
    """Write the experiment_symmetry.ini file for the ViPErLEED GUI.

    The experiment_symmetry.ini file can be used as input for the
    ViPErLEED pattern-simulator GUI to display the expected LEED
    pattern and export a "spot-pattern file" to be used in the
    ImageJ spot tracker.

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
        If `slab` has no `.bulkslab` yet.
    ValueError
        If `rpars.SUPERLATTICE` is not (close to) integer.
    OSError
        If writing to file fails.
    """
    if slab.bulkslab is None:
        _LOGGER.error('experiment_symmetry.ini: bulk '
                      'slab has not been initialized.')
        raise MissingBulkSlabError('experiment_symmetry.write '
                                   'called without bulk slab.')
    as_dict = getLEEDdict(slab, rpars)
    if as_dict is None:
        raise ValueError('SUPERLATTICE is not integer-valued.')

    # Convert arrays to lists for one-line display
    for key in ('SUPERLATTICE', 'surfBasis'):
        as_dict[key] = as_dict[key].tolist()

    if rpars.systemName:
        as_dict['name'] = rpars.systemName                                      # TODO: adapt for domains
    try:
        LEEDParser.write_structures(as_dict, FILENAME)
    except OSError:
        _LOGGER.error(f'Failed to write {FILENAME!r}')
        raise
    _LOGGER.debug(f'Wrote to {FILENAME!r} successfully')
