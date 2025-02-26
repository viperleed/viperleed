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

import logging
from pathlib import Path
import shutil
from zipfile import BadZipFile

from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.files import iotensors
from viperleed.calc.lib import leedbase


_LOGGER = logging.getLogger(__name__)
_DOMAIN_INPUT_FILES = (
    'PARAMETERS',
    'PHASESHIFTS',
    'POSCAR',
    'VIBROCC',
    )


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
        self.tensorDir = None                                                   # TODO: this seems unused. @fkraushofer: should we remove it?

    def __str__(self):
        """Return a string representation for this domain."""
        return f'domain {self.name}'

    def collect_input_files(self, src):
        """Fetch input files from the source given by the user."""
        _LOGGER.info(f'Fetching input files for {self}')
        if not src.is_dir() and not src.is_file():
            return
        _collect = (self._collect_inputs_from_directory if src.is_dir()
                    else self._collect_inputs_from_tensor_file)
        try:
            _collect(src)
        except (OSError, BadZipFile) as exc:
            raise RuntimeError('Error getting domain input files') from exc

    def _collect_inputs_from_directory(self, src):
        """Fetch input files from a user-given source directory."""
        # Try first to pull the inputs from the most recent tensor file
        try:
            tensor_dir = self._collect_inputs_from_most_recent_tensor(src)
        except (OSError, BadZipFile):
            pass
        else:
            self.tensorDir = tensor_dir
            return

        # No usable tensors in src; fetch inputs from src directly
        self.refcalcRequired = True
        _LOGGER.info(f'No previous {DEFAULT_TENSORS} found, '
                     'reference calculation is required.')
        may_auto_generate = {'PHASESHIFTS'}
        for file in _DOMAIN_INPUT_FILES:
            try:
                shutil.copy2(src / file, self.workdir)
            except FileNotFoundError:
                if file in may_auto_generate:
                    continue
                _LOGGER.error(f'Required file {file} for {self} '
                              f'not found in origin folder {src}')
                raise
            except OSError:
                if file in may_auto_generate:
                    continue
                _LOGGER.error(f'Error copying required file {file} for '
                              f'{self} from origin folder {src}')
                raise

    def _collect_inputs_from_most_recent_tensor(self, src):
        """Fetch the most recent tensor at src and copy its input files."""
        tensor_index = leedbase.getMaxTensorIndex(src)
        if not tensor_index:
            raise FileNotFoundError(f'No {DEFAULT_TENSORS} at {src}')
        # Unpack the most recent tensor at self.workdir
        try:
            iotensors.fetch_unpacked_tensor(tensor_index, base_dir=src,
                                            target_dir=self.workdir)
        except (OSError, BadZipFile) as exc:
            _LOGGER.warning(f'Error fetching {DEFAULT_TENSORS}: {exc}')
            raise
        # Finally, pull the input files there into the domain root
        tensor_dir = (
            self.workdir
            / f'{DEFAULT_TENSORS}/{DEFAULT_TENSORS}_{tensor_index:03d}'
            )
        try:
            self._collect_inputs_from_tensor_folder(tensor_dir)
        except FileNotFoundError as exc:
            file = Path(exc.filename).name
            _LOGGER.warning(f'Input file {file} is missing in '
                            f'{DEFAULT_TENSORS} directory. A new '
                            'reference calculation is required.')
            raise
        return tensor_dir

    def _collect_inputs_from_tensor_file(self, src_zip):
        """Fetch input files from a user-given Tensor file."""
        tensor_index = leedbase.getMaxTensorIndex(self.workdir)
        tensor_dir = (
            self.workdir
            / f'{DEFAULT_TENSORS}/{DEFAULT_TENSORS}_{tensor_index + 1:03d}'
            )
        try:
            iotensors.unpack_tensor_file(src_zip, tensor_dir)
        except (OSError, BadZipFile):
            _LOGGER.error(f'Failed to unpack {DEFAULT_TENSORS} '
                          f'for {self} from file {src_zip}')
            raise
        try:
            self._collect_inputs_from_tensor_folder(tensor_dir)
        except FileNotFoundError as exc:
            file = Path(exc.filename).name
            _LOGGER.error(f'Required file {file} for {self} not found '
                          f'in {DEFAULT_TENSORS} directory {tensor_dir}')
            raise
        self.tensorDir = tensor_dir

    def _collect_inputs_from_tensor_folder(self, tensor_dir):
        """Fetch input files from an unpacked tensor_dir."""
        for file in (_DOMAIN_INPUT_FILES + ('IVBEAMS',)):
            shutil.copy2(tensor_dir / file, self.workdir)
