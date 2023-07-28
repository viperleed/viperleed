"""Test section refcalc

Created on 2023-07-28

@author: Alexander M. Imre
"""
from pathlib import Path

import pytest


class TestDeltasAg100():
    def test_delta_input_written(self, delta_files_ag100):
        assert delta_files_ag100.expected_file_exists("delta-input")


    def test_exit_code_0(self, delta_files_ag100):
        assert delta_files_ag100.exit_code == 0


    def test_deltas_zip_created(self, delta_files_ag100):
        assert delta_files_ag100.expected_file_exists(Path("Deltas") / "Deltas_001.zip")

