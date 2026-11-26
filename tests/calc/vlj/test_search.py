"""Tests for module viperleed.calc.vlj.search."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-26'
__license__ = 'GPLv3+'

import logging
import subprocess

import pytest

from importlib import import_module
from importlib.util import find_spec

import pytest

from viperleed.calc.files import vibrocc
from viperleed.calc.files.new_displacements.file import DisplacementsFile
from viperleed.calc.symmetry import enforceSymmetry
from viperleed.calc.symmetry import findSymmetry
from viperleed.calc.vlj import search

# The plugin is optional, and importing it swaps the numerical backend
# of the R-factor library: import it only inside the tests that need it
pytestmark = pytest.mark.skipif(find_spec('viperleed_jax') is None,
                                reason='viperleed_jax plugin not installed')

_DISPLACEMENTS = '''\
== SEARCH z

= GEO_DELTA
Ag L(1-4)   z = -0.2 0.2 0.01
'''


@pytest.fixture(name='tree')
def fixture_tree():
    """Return the plugin's module with the parameter-space tree."""
    return import_module('viperleed_jax.transformation_tree.tree')


@pytest.fixture(name='ag100_parameter_space')
def fixture_ag100_parameter_space(ag100, data_path, tmp_path):
    """Return a plugin ParameterSpace for Ag(100) with a search segment."""
    slab, rpars, *_ = ag100
    vibrocc.readVIBROCC(rpars, slab,
                        str(data_path / 'Ag(100)/mergeDisp/VIBROCC'))
    slab.full_update(rpars)
    findSymmetry(slab, rpars)
    enforceSymmetry(slab, rpars)
    from_objects = import_module('viperleed_jax.from_objects')
    parameter_space = from_objects.setup_tl_parameter_space(slab, rpars)
    displacements = tmp_path / 'DISPLACEMENTS'
    displacements.write_text(_DISPLACEMENTS)
    disp_file = DisplacementsFile()
    disp_file.read(displacements)
    parameter_space.apply_search_segment(disp_file.next(0.9))
    return parameter_space


class TestExportParameterSpace:
    """Tests for the graphical export of the plugin's ParameterSpace."""

    def test_export_called(self, ag100_parameter_space, tree, monkeypatch):
        """Check that the plugin's export writes the expected file."""
        written = []
        monkeypatch.setattr(tree.UniqueDotExporter, 'to_picture',
                            lambda _, filename: written.append(filename))
        search._export_parameter_space(ag100_parameter_space, 'z')
        assert written == ['parameter_space_z.pdf']

    _dot_errors = {
        'dot missing': FileNotFoundError('dot'),
        'dot failed': subprocess.CalledProcessError(1, 'dot'),
        }

    @pytest.mark.parametrize('error', _dot_errors.values(), ids=_dot_errors)
    def test_dot_failure_warns(self, error, ag100_parameter_space, tree,
                               monkeypatch, caplog):
        """Check that a failure of the `dot` tool only logs a warning."""
        def _fail(*_):
            raise error
        monkeypatch.setattr(tree.UniqueDotExporter, 'to_picture', _fail)
        with caplog.at_level(logging.WARNING):
            search._export_parameter_space(ag100_parameter_space, 'z')
        assert 'Could not export parameter space tree view' in caplog.text

    def test_interface_mismatch_raises(self):
        """Check that a missing export method is not swallowed."""
        with pytest.raises(AttributeError):
            search._export_parameter_space(object(), 'z')
