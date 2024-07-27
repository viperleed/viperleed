"""ViPErLEED utility: Attach Bulk

Takes a slab POSCAR and adds a bulk POSCAR on the bottom, rescaling the
unit cell. Very primitive script, should be updated to include more
recent functionality.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2019-06-14'
__license__ = 'GPLv3+'

import copy
import logging
import time

import numpy as np

from viperleed.calc.files import poscar
from viperleed.calc.lib.time_utils import ExecutionTimer
from viperleed.cli_base import ViPErLEEDCLI


_LOG_FORMAT = '%(levelname)s - %(message)s'
_LOG_HEAD = '''\
Starting new log: {logname}
Time of execution (UTC): {now}
'''
_LOG_TAIL = '''\
Finishing {logname} at {now}.
Elapsed time: {delta_t}
'''
_CONSOLE_HANDLER = logging.StreamHandler()
_CONSOLE_HANDLER.setFormatter(logging.Formatter(_LOG_FORMAT))


class AttachBulkCLI(ViPErLEEDCLI, cli_name='attach_bulk'):
    """Main command-line interface of this utility."""

    def __call__(self, _=None):
        """Call this utility."""
        try:
            return self._call_impl()
        except KeyboardInterrupt:
            return 1
        finally:
            _tear_down_logger()

    @staticmethod
    def _call_impl():
        """The actual implementation of the __call__ method."""
        logname = 'Combine-POSCAR.log'
        _set_up_logger(logname)
        timer = ExecutionTimer()

        # pylint: disable-next=logging-format-interpolation
        logging.info(_LOG_HEAD.format(logname=logname, now=_now()))
        try:
            slab, bulk = _read_both_poscars_from_user_input()
        except poscar.POSCARError:
            # Error reading. Already reported
            return 2

        try:  # If a and b vectors differ, cancel operation
            _check_ab_consistent(slab, bulk)
        except ValueError as exc:
            logging.error(f'{exc} Stopping operation...')
            return 2

        _attach_bulk(slab, bulk)

        try:
            poscar.write(slab)
        except OSError:
            logging.error('Exception while writing combined POSCAR:',
                          exc_info=True)
        # pylint: disable-next=logging-format-interpolation
        logging.info(_LOG_TAIL.format(logname=logname,
                                      now=_now(),
                                      delta_t=timer.how_long(as_string=True)))
        return 0


def _add_missing_elements(slab, bulk):
    """Add to slab the extra elements of bulk."""
    if slab.elements == bulk.elements:
        logging.debug('Slab and bulk elements are identical.')
        return
    logging.debug('Slab and bulk elements are not '
                  'equal. Adding missing elements.')
    for element in bulk.elements:
        slab.n_per_elem[element] = 0


def _attach_bulk(slab, bulk):
    """Attach bulk at the bottom of slab."""
    _add_missing_elements(slab, bulk)

    # Resize the slab unit cell
    cfact = slab.c_vector[2]/bulk.c_vector[2]
    slab.c_vector[:] *= (cfact + 1) / cfact
    # Recalculate c for the slab atoms (undistort & shift)
    for atom in slab:
        atom.pos[2] = (atom.pos[2]*(cfact/(cfact+1)))+(1/(cfact+1))
        atom.num += bulk.n_atoms
    # Copy atoms from bulk and add them to the slab
    for atom in bulk:
        newat = copy.copy(atom)
        # Recalculate bulk atom c in the new slab unit cell
        newat.pos[2] /= cfact+1
        slab.atlist.append(newat)
        newat.slab = slab
    slab.atlist.update_atoms_map()
    slab.update_element_count()
    slab.sort_by_element()


def _check_ab_consistent(slab, bulk):
    """Raise if a&b are not identical within 1e-4."""
    eps = 1e-4
    ab_diff = slab.ab_cell.T - bulk.ab_cell.T
    if any(np.linalg.norm(ab_diff, axis=1) > eps):
        raise ValueError('Slab and bulk unit cell vectors are not '
                         f'equal in a and b (error > {eps:g}).')
    logging.debug('Slab and bulk unit cell vectors are '
                  f'equal in a and b (error < {eps:g}).')


def _now():
    """Return the current time as a string."""
    return time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())


def _read_both_poscars_from_user_input():
    """Return slab and bulk POSCARs to be merged."""
    slab = _read_poscar_from_user_input('slab')
    bulk = _read_poscar_from_user_input('bulk')
    return slab, bulk


def _read_poscar_from_user_input(name):
    """Return a slab read from a file specified by the user."""
    while True:
        filename = input(f'Enter {name} POSCAR name (Ctrl+C to abort): ')
        if not filename:
            print('Input failed. Please try again.')
            continue
        try:
            slab = poscar.read(filename)
        except FileNotFoundError:
            logging.error(f'{filename} not found.')
        except poscar.POSCARError:
            logging.error('Exception while reading POSCAR', exc_info=True)
            raise
        else:
            logging.info(f'{name.capitalize()} POSCAR was read successfully.')
            return slab


def _set_up_logger(logname):
    """Prepare the logging module to log to logname and the console."""
    logging.basicConfig(level=logging.DEBUG, filename=logname,
                        filemode='w', format=_LOG_FORMAT)
    logging.getLogger().addHandler(_CONSOLE_HANDLER)


def _tear_down_logger():
    """Remove handlers and shutdown the logger."""
    logging.getLogger().removeHandler(_CONSOLE_HANDLER)
    logging.shutdown()


if __name__ == '__main__':
    AttachBulkCLI.run_as_script()
