"""Test section search

Created on 2023-07-28

@author: Alexander M. Imre
"""

import pytest


class TestSearchAg100():
    def test_exit_code_0(self, search_files_ag100):
        assert search_files_ag100.exit_code == 0

    @pytest.mark.parametrize('expected_file', ('search.steu',))
    def test_search_input_exist(self, search_files_ag100, expected_file):
        assert search_files_ag100.expected_file_exists(expected_file)

    @pytest.mark.parametrize('expected_file', ('SD.TL', 'control.chem'))
    def test_search_raw_files_exist(self, search_files_ag100, expected_file):
        assert search_files_ag100.expected_file_exists(expected_file)

    @pytest.mark.parametrize('expected_file', ('Search-report.pdf', 'Search-progress.pdf'))
    def test_search_pdf_files_exist(self, search_files_ag100, expected_file):
        assert search_files_ag100.expected_file_exists(expected_file)
