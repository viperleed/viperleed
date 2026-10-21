"""Tests for the --rng-seed option of the test-suite configuration."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-15'
__license__ = 'GPLv3+'

import os
from pathlib import Path
import re
import subprocess
import sys
import zlib

import pytest

from .conftest import N_DRAWS

# The tests in this module start pytest subprocesses
pytestmark = pytest.mark.slow

pytest.importorskip('xdist')

_TESTS_PARENT = Path(__file__).resolve().parent.parent
_TEST_SOURCE = """
import os

def test_one(draws):
    with open(os.environ['SEED_LOG'], 'a') as log:
        for draw in draws:
            log.write(f'{draw.index} {draw.seed}\\n')
"""
_BASE_PATTERN = re.compile(r'rng-seed base: (\d+) ')


def _run_seed_logging_pytest(tmp_path, name, *args):
    """Run a random-input test in a subprocess. Return (base, seeds)."""
    folder = tmp_path / name
    folder.mkdir()
    (folder / 'test_small.py').write_text(_TEST_SOURCE)
    log = folder / 'seeds.log'
    python_path = [str(_TESTS_PARENT)]
    if os.environ.get('PYTHONPATH'):
        python_path.append(os.environ['PYTHONPATH'])
    env = dict(os.environ,
               SEED_LOG=str(log),
               PYTHONPATH=os.pathsep.join(python_path))
    env.pop('PYTEST_ADDOPTS', None)  # Could disable xdist
    result = subprocess.run(  # noqa: S603  # Fixed, trusted arguments
        (
            sys.executable,
            '-m',
            'pytest',
            '-p',
            'tests.conftest',
            '--no-cov',
            '-p',
            'no:cacheprovider',
            str(folder),
            *args,
            ),
        cwd=folder,
        env=env,
        capture_output=True,
        text=True,
        check=False,
        )
    assert not result.returncode, result.stdout + result.stderr
    match = _BASE_PATTERN.search(result.stdout)
    assert match
    return int(match.group(1)), sorted(log.read_text().splitlines())


def test_random_base_reproducible_under_xdist(request, tmp_path):
    """Check that the printed random base reproduces a parallel run."""
    if not request.config.pluginmanager.hasplugin('xdist'):
        pytest.skip('pytest-xdist is not active')
    base, seeds = _run_seed_logging_pytest(
        tmp_path,
        'parallel',
        '--rng-seed=random',
        '-n',
        '2',
        )
    base_serial, seeds_serial = _run_seed_logging_pytest(
        tmp_path,
        'serial',
        f'--rng-seed={base}',
        )
    assert base_serial == base
    assert len(seeds) == N_DRAWS
    assert seeds == seeds_serial
    expected = sorted(
        f'{k} {zlib.crc32(f"{base}:test_one#{k}".encode())}'
        for k in range(N_DRAWS)
        )
    assert seeds == expected
