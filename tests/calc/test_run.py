"""Tests for module run of viperleed.calc."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-03-05'
__license__ = 'GPLv3+'

import logging
from pathlib import Path

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.files import parameters
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.run import _get_parent_directory_name
from viperleed.calc.run import _make_rpars_and_slab
from viperleed.calc.run import _interpret_parameters
from viperleed.calc.run import _read_parameters_file
from viperleed.calc.run import _read_poscar_file
from viperleed.calc.run import _set_log_level
from viperleed.calc.run import _set_tensorleed_source
from viperleed.calc.run import _set_system_name
from viperleed.calc.run import run_calc

_MODULE = 'viperleed.calc.run'


class TestRunCalc:
    """Tests for the run_calc function."""

    @staticmethod
    def check_calls(mocks, expected_calls, not_called, n_log_msgs):
        """Ensure mocks have been called or not."""
        for mock_name, call in expected_calls.items():
            assert mocks[mock_name].mock_calls == [call]
        for mock_name in not_called:
            mocks[mock_name].assert_not_called()
        assert len(mocks['logger'].info.mock_calls) == n_log_msgs
        assert len(mocks['logger'].method_calls) == n_log_msgs

    @fixture(name='mocks')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details with mocks."""
        rpars = Rparams()
        slab = mocker.MagicMock()
        mock_manifest = mocker.MagicMock()
        mock_now = mocker.MagicMock(return_value='some_timestamp')
        mock_fmt = mocker.patch(f'{_MODULE}.DateTimeFormat')
        mock_fmt.FILE_SUFFIX.now = mock_now
        mock_fmt.LOG_CONTENTS.now = mock_now
        mocks = {
            'umask': mocker.patch('os.umask'),
            'now': mock_now,
            'logger': mocker.patch(f'{_MODULE}.logger'),
            'prepare_log': mocker.patch(f'{_MODULE}.prepare_calc_logger'),
            'manifest': mocker.patch(f'{_MODULE}.ManifestFile',
                                     return_value=mock_manifest),
            'make_rp_sl': mocker.patch(f'{_MODULE}._make_rpars_and_slab',
                                       return_value=(rpars, slab)),
            'cleanup': mocker.patch(f'{_MODULE}.cleanup'),
            'set_src': mocker.patch(f'{_MODULE}._set_tensorleed_source'),
            'set_name': mocker.patch(f'{_MODULE}._set_system_name'),
            'find_tl_version': mocker.patch.object(rpars,
                                                   'updateDerivedParams'),
            'prerun': mocker.patch(f'{_MODULE}.prerun_clean'),
            'preserve': mocker.patch(f'{_MODULE}.preserve_original_inputs'),
            'section': mocker.patch(f'{_MODULE}.section_loop',
                                    return_value=(0, 'recorder')),
            'handlers': mocker.patch(f'{_MODULE}.close_all_handlers'),
            'shutdown': mocker.patch(f'{_MODULE}.logging.shutdown'),
            }
        return mocks

    _user_args = (
        {'system_name': 'some name'},
        {'console_output': 'log to console'},
        {'slab': 'user-given slab'},
        {'preset_params': 'some parameter presets'},
        {'source': 'path to tensorleed'},
        {'home': 'where calc started'},
        )

    @parametrize(kwargs=_user_args)
    def test_success_calls(self, kwargs, mocks, mocker):
        """Check the expected calls for a successful run."""
        rpars, slab = mocks['make_rp_sl'].return_value
        timestamp = mocks['now']()
        log_name = f'viperleed-calc-{timestamp}.log'
        expect_calls = {
            'umask': mocker.call(0),
            'prepare_log': mocker.call(
                mocks['logger'],
                file_name=log_name,
                with_console=kwargs.get('console_output', True),
                ),
            'manifest': mocker.call('SUPP', 'OUT', log_name),
            'make_rp_sl': mocker.call(
                mocks['manifest'].return_value,
                kwargs.get('preset_params'),
                kwargs.get('slab'),
                kwargs.get('home'),
                ),
            'set_src': mocker.call(rpars, kwargs.get('source')),
            'set_name': mocker.call(rpars, kwargs.get('system_name')),
            'find_tl_version': mocker.call(),
            'prerun': mocker.call(rpars, log_name),
            'preserve': mocker.call(rpars),
            'section': mocker.call(rpars, slab),
            'handlers': mocker.call(mocks['logger']),
            'shutdown': mocker.call(),
            }
        not_called = ('cleanup',)  # Only in error cases

        result = run_calc(**kwargs)
        assert result == mocks['section'].return_value
        assert rpars.timestamp == timestamp
        assert rpars.manifest is mocks['manifest'].return_value

        self.check_calls(mocks, expect_calls, not_called, n_log_msgs=3)

    _raises = (
        Exception,
        FileNotFoundError,
        OSError,
        parameters.errors.ParameterError,
        TypeError,
        )

    @parametrize(exc=_raises)
    def test_read_inputs_fails(self, exc, mocks, mocker):
        """Check behavior when reading inputs fails."""
        mocks['make_rp_sl'].side_effect = exc
        exit_code, states = run_calc()
        assert exit_code
        assert not states

        timestamp = mocks['now']()
        log_name = f'viperleed-calc-{timestamp}.log'
        expect_calls = {
            'umask': mocker.call(0),
            'prepare_log': mocker.call(mocks['logger'],
                                       file_name=log_name,
                                       with_console=True),
            'manifest': mocker.call('SUPP', 'OUT', log_name),
            'make_rp_sl': mocker.call(mocks['manifest'].return_value,
                                      None,   # preset_params
                                      None,   # slab
                                      None),  # home
            'cleanup': mocker.call(mocks['manifest'].return_value),
            }
        not_called = (
            'set_src',
            'set_name',
            'find_tl_version',
            'prerun',
            'preserve',
            'section',
            'handlers',
            'shutdown',
            )
        self.check_calls(mocks, expect_calls, not_called, n_log_msgs=2)

    def test_halt_active(self, mocks, mocker):
        """Check early exit if halting is in effect."""
        rpars, _ = mocks['make_rp_sl'].return_value
        timestamp = mocks['now']()
        rpars.halt = 2  # Should stop right away with default HALTING

        exit_code, states = run_calc()
        assert not exit_code
        assert not states
        assert rpars.timestamp == timestamp
        assert rpars.manifest is mocks['manifest'].return_value

        log_name = f'viperleed-calc-{timestamp}.log'
        expect_calls = {
            'umask': mocker.call(0),
            'prepare_log': mocker.call(mocks['logger'],
                                       file_name=log_name,
                                       with_console=True),
            'manifest': mocker.call('SUPP', 'OUT', log_name),
            'make_rp_sl': mocker.call(mocks['manifest'].return_value,
                                      None,   # preset_params
                                      None,   # slab
                                      None),  # home
            'set_src': mocker.call(rpars, None),
            'set_name': mocker.call(rpars, None),
            'cleanup': mocker.call(rpars),
            }
        not_called = (
            'find_tl_version',
            'prerun',
            'preserve',
            'section',
            'handlers',
            'shutdown',
            )
        self.check_calls(mocks, expect_calls, not_called, n_log_msgs=3)


def test_get_parent_directory_name(mocker, caplog):
    """Test correct result of the _get_parent_directory_name helper."""
    caplog.set_level(5)
    fake_cwd = Path('/home/user/project')
    mocker.patch('pathlib.Path.cwd', return_value=fake_cwd)
    expect_name = 'user'
    assert _get_parent_directory_name() == expect_name
    assert expect_name in caplog.text


class TestMakeRparsAndSlab:
    """Tests for the _make_rpars_and_slab helper."""

    @fixture(name='mocks')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details with mocks."""
        mock_slab = mocker.MagicMock()
        mock_rpars = Rparams()
        return {
            'rpars': mock_rpars,
            'slab': mock_slab,
            'read_param': mocker.patch(f'{_MODULE}._read_parameters_file',
                                       return_value=mock_rpars),
            'read_poscar': mocker.patch(f'{_MODULE}._read_poscar_file',
                                        return_value=mock_slab),
            'interpret': mocker.patch(f'{_MODULE}._interpret_parameters'),
            'warn_slab': mocker.patch(
                f'{_MODULE}.warn_if_slab_has_atoms_in_multiple_c_cells',
                ),
            }

    def test_domains(self, mocks):
        """Check correct outcome for a multi-domain case."""
        manifest = set()
        presets = None
        mocks['rpars'].readParams['DOMAIN'] = 'some domain'
        rpars, slab = _make_rpars_and_slab(manifest,
                                           presets,
                                           slab=None,
                                           home='here')
        assert rpars is mocks['rpars']
        assert slab is None
        assert not rpars.fileLoaded['POSCAR']
        assert rpars.paths.home == Path('here')
        args = {
            'read_param': (presets,),
            'interpret': (rpars, slab, {}),
            }
        for mock_name, mock_args in args.items():
            mocks[mock_name].assert_called_once_with(*mock_args)
        mocks['read_poscar'].assert_not_called()
        mocks['warn_slab'].assert_not_called()
        mocks['slab'].full_update.assert_not_called()

    def test_domains_raises(self, mocks):
        """Check correct outcome for a multi-domain case."""
        manifest = set()
        presets = None
        mocks['rpars'].readParams['DOMAIN'] = 'some domain'
        with pytest.raises(TypeError):
            _make_rpars_and_slab(manifest, presets, slab='a slab', home=None)

        mocks['read_param'].assert_called_once_with(presets)
        not_called = 'interpret', 'read_poscar', 'warn_slab'
        for mock_name in not_called:
            mocks[mock_name].assert_not_called()

    _raises = (
        ('read_param', Exception),
        ('read_param', FileNotFoundError),
        ('read_poscar', Exception),
        ('read_poscar', FileNotFoundError),
        ('read_poscar', OSError),
        ('interpret', parameters.errors.ParameterError),
        )

    @parametrize('mock_name,exc', _raises)
    def test_implementation_raises(self, mock_name, exc, mocks):
        """Check that exceptions raised by the implementation bubble out."""
        mocks[mock_name].side_effect = exc('test exc')
        with pytest.raises(exc, match='test exc'):
            _make_rpars_and_slab({}, None, None, None)

    def test_one_domain_no_slab(self, mocks):
        """Check correct outcome for a single domain and no slab given."""
        manifest = set()
        presets = {'some_preset': 1}
        rpars, slab = _make_rpars_and_slab(manifest,
                                           presets,
                                           slab=None,
                                           home=None)
        assert rpars is mocks['rpars']
        assert slab is mocks['slab']
        assert rpars.fileLoaded['POSCAR']
        assert rpars.paths.home
        args = {
            'read_param': (presets,),
            'read_poscar': (manifest,),
            'interpret': (rpars, slab, presets),
            'warn_slab': (slab, rpars),
            }
        for mock_name, mock_args in args.items():
            mocks[mock_name].assert_called_once_with(*mock_args)
        mocks['slab'].full_update.assert_called_once_with(rpars)

    def test_one_domain_with_slab(self, mocks, mocker):
        """Check correct outcome for a single domain and a user-given slab."""
        manifest = set()
        presets = {'some_preset': 1}
        slab_arg = mocker.MagicMock()
        rpars, slab = _make_rpars_and_slab(manifest,
                                           presets,
                                           slab=slab_arg,
                                           home=None)
        assert slab is slab_arg
        assert rpars.fileLoaded['POSCAR']
        assert rpars.paths.home
        args = {
            'read_param': (presets,),
            'interpret': (rpars, slab, presets),
            'warn_slab': (slab, rpars),
            }
        for mock_name, mock_args in args.items():
            mocks[mock_name].assert_called_once_with(*mock_args)
        mocks['read_poscar'].assert_not_called()
        slab.full_update.assert_called_once_with(rpars)


class TestInterpretParameters:
    """Tests for the _interpret_parameters helper."""

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details with mocks."""
        mock_slab = mocker.MagicMock()
        mock_rpars = Rparams()
        mock_presets = {'LOG_LEVEL': 'DEBUG'}
        mocks = {
            'interpret': mocker.patch(f'{_MODULE}.parameters.interpret'),
            'set_log': mocker.patch(f'{_MODULE}._set_log_level'),
            'update': mocker.patch.object(mock_rpars, 'update'),
            }
        return mock_rpars, mock_slab, mock_presets, mocks

    def test_interpret_and_update(self, mock_implementation, mocker, caplog):
        """Check correct interpretation and update of parameters."""
        caplog.set_level(logging.INFO)
        rpars, slab, presets, mocks = mock_implementation
        _interpret_parameters(rpars, slab, presets)
        calls = {
            'interpret': mocker.call(rpars, slab=slab, silent=False),
            'update': mocker.call(presets),
            'set_log': mocker.call(rpars, presets),
            }
        for mock_name, expect_call in calls.items():
            assert mocks[mock_name].mock_calls == [expect_call]
        assert not caplog.text

    interpret_raises = (
        parameters.errors.ParameterNeedsSlabError,
        parameters.errors.SuperfluousParameterError,
        parameters.errors.ParameterError,
        )

    @parametrize(exc=interpret_raises)
    def test_interpret_raises(self, exc, mock_implementation, caplog):
        """Check bubble-out of exceptions in implementation details."""
        *args, mocks = mock_implementation
        mocks['interpret'].side_effect = exc('test exc')
        with pytest.raises(exc, match='test exc'):
            _interpret_parameters(*args)
        assert caplog.text

    @parametrize(exc=(ValueError, TypeError))
    def test_apply_presets_fails(self, exc, mock_implementation, caplog):
        """Check that logs are emitted when applying presets fails."""
        *args, mocks = mock_implementation
        mocks['update'].side_effect = exc
        _interpret_parameters(*args)  # No exception here
        assert caplog.text


class TestReadParametersFile:
    """Tests for the _read_parameters_file helper."""

    @fixture(name='mocks')
    def fixture_mock(self, mocker):
        """Replace implementation details with mocks."""
        mock_rpars = Rparams()
        mock_read = mocker.patch(f'{_MODULE}.parameters.read',
                                 return_value=mock_rpars)
        return mock_rpars, mock_read

    def test_reads_existing_parameters(self, mocks, caplog):
        """Check correct result when PARAMETERS is found."""
        caplog.set_level(5)  # < DEBUG
        rpars, mock_read = mocks
        result = _read_parameters_file('some unused preset')
        mock_read.assert_called_once()
        assert result is rpars
        assert not caplog.text

    def test_not_found_no_presets(self, mocks, caplog):
        """Check complaints with neither PARAMETERS nor presets."""
        _, mock_read = mocks
        mock_read.side_effect = FileNotFoundError
        with pytest.raises(FileNotFoundError):
            _read_parameters_file(None)
        assert caplog.text

    def test_not_found_with_presets(self, mocks, caplog):
        """Check no complaints with no PARAMETERS but user presets."""
        caplog.set_level(5)  # < DEBUG
        _, mock_read = mocks
        mock_read.side_effect = FileNotFoundError
        result = _read_parameters_file('some preset')
        assert isinstance(result, Rparams)
        assert not caplog.text

    def test_read_fails(self, mocks, caplog):
        """Check exceptions while reading PARAMETERS bubble out."""
        _, mock_read = mocks
        mock_read.side_effect = Exception('test exc')
        with pytest.raises(Exception, match='test exc'):
            _read_parameters_file(None)
        assert caplog.text


class TestReadPoscarFile:
    """Tests for the _read_poscar_file helper."""

    @fixture(name='mocks')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details with mocks."""
        mock_slab = mocker.MagicMock()
        return {
            'slab': mock_slab,
            'read': mocker.patch(f'{_MODULE}.poscar.read',
                                 return_value=mock_slab),
            'copy': mocker.patch('shutil.copy2'),
            }

    def test_success(self, mocks, caplog):
        """Test successful reading of a POSCAR file."""
        caplog.set_level(logging.INFO)
        manifest = set()
        mocks['slab'].preprocessed = True
        slab = _read_poscar_file(manifest)
        expect_log = 'Reading structure'

        assert slab is mocks['slab']
        mocks['read'].assert_called_once_with(filename=Path('POSCAR'))
        mocks['copy'].assert_not_called()
        assert not manifest
        assert expect_log in caplog.text

    def test_poscar_user(self, mocks, caplog):
        """Test successful reading of a "fresh" POSCAR file."""
        caplog.set_level(logging.INFO)
        manifest = set()
        mocks['slab'].preprocessed = False
        slab = _read_poscar_file(manifest)
        expect_log = 'Reading structure'

        assert slab is mocks['slab']
        mocks['read'].assert_called_once_with(filename=Path('POSCAR'))
        mocks['copy'].assert_called_once_with(Path('POSCAR'), 'POSCAR_user')
        assert manifest
        assert expect_log in caplog.text

    _raises = (
        ('read', FileNotFoundError),
        ('read', Exception),
        ('copy', OSError),
        )

    @parametrize('mock_name,exc', _raises)
    def test_raises(self, mock_name, exc, mocks, caplog):
        """Check complaints when the implementation raises exceptions."""
        caplog.set_level(logging.ERROR)
        mocks[mock_name].side_effect = exc('test exc')
        mocks['slab'].preprocessed = False
        with pytest.raises(exc, match='test exc'):
            _read_poscar_file(None)
        assert caplog.text


class TestSetLogLevel:
    """Tests for the _set_log_level helper."""

    @fixture(name='mocks')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details with mocks."""
        rpars = Rparams()
        rpars.LOG_LEVEL = 37
        mock_logger = mocker.patch(f'{_MODULE}.logger')
        return rpars, mock_logger

    def test_preset(self, mocks):
        """Test setting log level from preset parameters."""
        rpars, mock_logger = mocks
        presets = {'LOG_LEVEL'}
        _set_log_level(rpars, presets)
        expect_log = 'Overriding log level to 37.'
        mock_logger.setLevel.assert_called_once_with(rpars.LOG_LEVEL)
        mock_logger.info.assert_called_once_with(expect_log)
        n_log_calls = 2  # setLevel and info
        assert len(mock_logger.method_calls) == n_log_calls

    def test_no_preset(self, mocks):
        """Test setting log level without presets."""
        rpars, mock_logger = mocks
        _set_log_level(rpars, {})
        mock_logger.setLevel.assert_called_once_with(rpars.LOG_LEVEL)
        assert len(mock_logger.method_calls) == 1  # only setLevel


class TestSetTensorleedSource:
    """Tests for the _set_tensorleed_source helper."""

    @fixture(name='mocks')
    def fixture_mock_implementation(self, mocker):
        """Replace internal implementation details with mocks."""
        mock_source = Path('/some/path/to/tensorleed')
        mock_rpars = Rparams()
        mock_get = mocker.patch(f'{_MODULE}.get_tensorleed_path',
                                return_value=mock_source)
        return {
            'src': mock_source.resolve(),
            'rpars': mock_rpars,
            'get': mock_get,
            }

    def test_success(self, mocks, caplog):
        """Check the successful setting of a source path."""
        caplog.set_level(5)  # < DEBUG
        rpars = mocks['rpars']
        src = 'some/source/path'
        _set_tensorleed_source(rpars, src)
        mocks['get'].assert_called_once_with(src)
        assert rpars.paths.tensorleed == mocks['src']
        assert not caplog.text

    _fails = {
        None: Path.cwd(),
        'some/source/path': Path('some/source/path'),
        }

    @parametrize(exc=(ValueError, FileNotFoundError))
    @parametrize('src,expect', _fails.items())
    # pylint: disable-next=too-many-arguments  # 2/6 fixtures
    def test_get_fails(self, exc, src, expect, mocks, caplog):
        """Check the successful setting of a source path."""
        caplog.set_level(logging.WARNING)
        rpars = mocks['rpars']
        mocks['get'].side_effect = exc
        _set_tensorleed_source(rpars, src)
        assert rpars.paths.tensorleed == expect
        assert caplog.text


class TestSetSystemName:
    """Tests for the _set_system_name helper."""

    def test_user_given(self, mocker):
        """Test setting a system name when explicitly provided."""
        rpars, explicit = (mocker.MagicMock() for _ in range(2))
        _set_system_name(rpars, explicit)
        assert rpars.systemName is explicit

    def test_from_parent_directory(self, mocker):
        """Test setting system name from the parent directory."""
        rpars, parent = (mocker.MagicMock() for _ in range(2))
        mock_get_parent = mocker.patch(f'{_MODULE}._get_parent_directory_name',
                                       return_value=parent)
        _set_system_name(rpars, None)
        assert rpars.systemName is parent
        mock_get_parent.assert_called_once()
