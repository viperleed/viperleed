"""Test configuration for tests/calc/sections/deltas."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-07-04'
__license__ = 'GPLv3+'


from pytest_cases import fixture

from viperleed.calc.classes.rparams.rparams import Rparams

from ....helpers import filesystem_from_dict

_MOCK_COMPILER = object()
_MODULE = 'viperleed.calc.sections.deltas'


@fixture(name='tenserleed_tree')
def fixture_tenserleed_tree(tmp_path):
    """Create a tree similar to the one of a TensErLEED source path."""
    mock_root = tmp_path/'tenserleed_root'
    tree = {
        'src': {
            'GLOBAL': 'global contents',
            'deltas.f': 'main deltas source',
            },
        'lib': {
            'lib.tleed.f': 'library of tensor-LEED functions',
            'lib.delta.f': 'library of delta-amplitudes functions',
            },
        }
    filesystem_from_dict(tree, mock_root)
    return mock_root


@fixture(name='rpars')
def fixture_rpars(tenserleed_tree, mocker):
    """Return an Rparams with basic attributes set."""
    rpars = Rparams()
    def _mock_get_fortran_comp():
        rpars.FORTRAN_COMP[0] = _MOCK_COMPILER

    rpars.disp_blocks = ['avoid IndexError']
    rpars.runHistory = [1]
    rpars.N_CORES = 5
    rpars.TENSOR_INDEX = 123456
    mocker.patch.object(rpars, 'updateCores')
    mocker.patch.object(rpars, 'setHaltingLevel')
    mocker.patch.object(rpars,
                        'getFortranComp',
                        side_effect=_mock_get_fortran_comp)
    mocker.patch.object(
        rpars,
        'get_tenserleed_directory',
        return_value=mocker.MagicMock(path=tenserleed_tree),
        )
    rpars.timestamp = '20250101'
    return rpars
