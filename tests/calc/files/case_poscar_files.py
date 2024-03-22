"""Collection of POSCAR files as cases for tests.

These used to be part of poscar_slabs. However, since they are not
slabs, they would mess up other tests (e.g., tests/symmetry).
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__created__ = '2024-01-08'

from pytest_cases import parametrize

from ...helpers import POSCAR_PATH

POSCAR_FILES = tuple(POSCAR_PATH.glob('POSCAR*'))


class CasePOSCARFiles:
    """Collection of POSCAR files. NOT the slab read from them."""

    @staticmethod
    def _string(file_name):
        """Return the string of a POSCAR file name."""
        with (POSCAR_PATH / file_name).open('r', encoding='utf-8') as file:
            return file.read()

    @parametrize(file=POSCAR_FILES)
    def case_poscar_file(self, file):
        """Return a POSCAR file name and its contents."""
        content = self._string(file)
        return file, content
