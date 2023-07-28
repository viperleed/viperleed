"""Test section refcalc

Created on 2023-07-28

@author: Alexander M. Imre
"""

import pytest


@pytest.mark.parametrize('expected_file', (('THEOBEAMS.csv',)))
def test_refcalc_files_present(self, refcalc_files, expected_file):
    assert refcalc_files.expected_file_exists(expected_file)
