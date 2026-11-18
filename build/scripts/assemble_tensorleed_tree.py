"""Script assemble_tensorleed_tree.

Assembles the TensErLEED tree that the tensorleed test tier runs
against (pytest -m tensorleed, with VIPERLEED_TENSORLEED pointing to
the tree). Used by .github/workflows/tensorleed_tests.yml and by
developers. Needs Internet access, gcc, gfortran, and make (Linux,
macOS) or make.bat (Windows, shipped with the EEASiSSS sources).

Usage: python assemble_tensorleed_tree.py TARGET_DIR

The tree is made from release tags of viperleed/viperleed-tensorleed,
downloaded as archives of the tagged commits (no git needed). The
module constants say which: TAG_COMMITS pins the commit of each tag,
TREE_CONTENTS lists the folders copied from each tag, and the
eeasisss executable is compiled from EEASISSS_TAG. In LEGACY_LIB, the
shipped random_.o and MPIrandom_.o are replaced by a random_.o
compiled on this machine, as the installation instructions of that
TensErLEED version require. It is compiled from RANDOM_SOURCE, the
cross-platform random_.c of the ViPErLEED documentation, which also
replaces the shipped random_.c in the tree.

TARGET_DIR must not exist or be empty: the script refuses to write
into a non-empty directory. Downloads and builds are made in a
temporary subfolder of TARGET_DIR, which is removed at the end. If a
step fails, everything the script wrote is removed, so that it can
be run again on the same TARGET_DIR.

Tested on Linux only.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-17'
__license__ = 'GPLv3+'

import argparse
import contextlib
import ntpath
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import urllib.request

ARCHIVE_URL = ('https://github.com/viperleed/viperleed-tensorleed/'
               'archive/{commit}.tar.gz')

EEASISSS_TAG = 'EEASiSSS-v1.0.2'

# Commits of the release tags. Pinned, so that a moved tag cannot
# change the tree silently.
TAG_COMMITS = {
    'v1.7.5': '71a4cd1f522fe7ab26bf19d09178e57bed544228',
    'v2.0.0': '286211173d560accaef4e046117351a188e3794c',
    'v2.1.0': '22e399826f67837378a95af077eaa3fb048d5cdf',
    EEASISSS_TAG: 'c8e3dbfb45e3bb9c3a9aadf37258346aae80bbaa',
    }

# (tag, path in the tag, path in the tree)
TREE_CONTENTS = (
    ('v1.7.5', 'TensErLEED-v1.61', 'TensErLEED-v1.61'),
    ('v1.7.5', 'TensErLEED-v1.74', 'TensErLEED-v1.74'),
    ('v1.7.5', 'TensErLEED-v1.75', 'TensErLEED-v1.75'),
    ('v2.0.0', 'TensErLEED', 'TensErLEED-v2.0.0'),
    ('v2.1.0', 'TensErLEED', 'TensErLEED-v2.1.0'),
    ('v2.0.0', 'atom_density_files', 'atom_density_files'),
    ('v2.0.0', 'seSernelius', 'seSernelius'),
    )

EEASISSS_MAKE = (
    ('cmd', '/c', 'make.bat', 'gcc') if sys.platform.startswith('win')
    else ('make', 'gcc')
    )
LEGACY_LIB = Path('TensErLEED-v1.61', 'lib')
# The shipped random_.c compiles on UNIX only. This one is used on
# all platforms, as the installation instructions suggest for Windows.
RANDOM_SOURCE = (Path(__file__).resolve().parents[2]
                 / 'doc/_static/tenserleed_windows/random_.c')
WORK_DIR_NAME = '_assemble_tmp'


class TreeAssemblyError(Exception):
    """The TensErLEED tree could not be assembled."""


def assemble_tree(target):
    """Assemble the TensErLEED test tree into a directory.

    Parameters
    ----------
    target : str or Path
        The directory to assemble the tree into. It is created if it
        does not exist.

    Raises
    ------
    TreeAssemblyError
        If target is not an empty directory, if an archive contains
        an unsafe member or more than one root folder, or if a
        compilation step fails.
    OSError
        If a download fails (urllib.error.URLError is an OSError) or
        a file cannot be written.
    tarfile.TarError
        If a downloaded archive cannot be read.
    """
    target = Path(target).resolve()
    if target.exists() and (not target.is_dir() or any(target.iterdir())):
        raise TreeAssemblyError(
            f'{target} exists and is not an empty directory. Remove it '
            'or choose another target.'
            )
    with _removed_on_failure(target):
        _assemble_tree(target)
    print(f'TensErLEED tree assembled in {target}')


def _assemble_tree(target):
    """Assemble the tree into the empty or missing directory target."""
    work = target / WORK_DIR_NAME
    work.mkdir(parents=True)
    sources = {tag: _download_tag(tag, work) for tag in TAG_COMMITS}
    for tag, src_path, tree_path in TREE_CONTENTS:
        src = sources[tag] / src_path
        print(f'Copying {tag}:{src_path} to {tree_path}')
        if src.is_dir():
            shutil.copytree(src, target / tree_path)
        else:
            shutil.copy2(src, target / tree_path)
    _compile_random(target / LEGACY_LIB)
    _compile_eeasisss(sources[EEASISSS_TAG], target)
    shutil.rmtree(work)


def _checked_members(tar, tag):
    """Return the members of tar, refusing unsafe ones.

    For Python versions without tarfile.data_filter. Refuse members
    that are not regular files or folders, or whose name is absolute,
    has a drive (also checked on POSIX, as C:name is unsafe on
    Windows) or contains '..'. Clear the setuid, setgid and sticky
    bits.
    """
    members = tar.getmembers()
    for member in members:
        drive, path = ntpath.splitdrive(member.name)
        parts = path.replace('\\', '/').split('/')
        if (not (member.isfile() or member.isdir())
                or drive or path.startswith(('/', '\\'))
                or os.pardir in parts):
            raise TreeAssemblyError(f'Unsafe member {member.name} in '
                                    f'the archive of {tag}')
        member.mode &= 0o777
    return members


def _compile_eeasisss(src, target):
    """Compile eeasisss in folder src and move it to target."""
    print(f'Compiling eeasisss from {EEASISSS_TAG}')
    _run(EEASISSS_MAKE, src)
    for name in ('eeasisss', 'eeasisss.exe'):
        executable = src / name
        if executable.is_file():
            shutil.copy2(executable, target / name)
            return
    raise TreeAssemblyError(f'Compiling eeasisss produced no executable '
                            f'in {src}')


def _compile_random(lib):
    """Compile RANDOM_SOURCE in place of the shipped randomizer."""
    print(f'Compiling random_.o in {lib} from {RANDOM_SOURCE}')
    for name in ('random_.o', 'MPIrandom_.o'):
        (lib / name).unlink()
    shutil.copyfile(RANDOM_SOURCE, lib / 'random_.c')
    _run(['gcc', '-O3', '-c', 'random_.c', '-o', 'random_.o'], lib)


def _download_tag(tag, work):
    """Download and unpack the archive of tag; return its root."""
    url = ARCHIVE_URL.format(commit=TAG_COMMITS[tag])
    print(f'Downloading {tag} from {url}')
    archive = work / f'{tag}.tar.gz'
    # A fixed https URL, no user input
    response = urllib.request.urlopen(url)  # noqa: S310
    with response, archive.open('wb') as file:
        shutil.copyfileobj(response, file)
    unpacked = work / tag
    with tarfile.open(archive) as tar:
        # Members are checked by the data filter or _checked_members
        if hasattr(tarfile, 'data_filter'):  # Python 3.12, backports
            tar.extractall(unpacked, filter='data')
        else:
            members = _checked_members(tar, tag)
            tar.extractall(unpacked, members=members)  # noqa: S202
    archive.unlink()
    roots = list(unpacked.iterdir())
    if len(roots) != 1:  # GitHub archives have one root folder
        raise TreeAssemblyError(f'The archive of {tag} has {len(roots)} '
                                'root folders instead of one')
    return roots[0]


@contextlib.contextmanager
def _removed_on_failure(target):
    """Remove what was written into the empty target if body fails."""
    target_existed = target.exists()
    try:
        yield
    except BaseException:  # Also on KeyboardInterrupt
        print(f'Removing the partial tree in {target}', file=sys.stderr)
        if target_existed:  # And was empty: keep only the folder
            for child in target.iterdir():
                if child.is_dir() and not child.is_symlink():
                    shutil.rmtree(child, ignore_errors=True)
                else:
                    child.unlink()
        else:
            shutil.rmtree(target, ignore_errors=True)
        raise


def _run(command, cwd):
    """Run command in folder cwd, raise TreeAssemblyError on failure."""
    try:
        subprocess.run(command, cwd=str(cwd), check=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        raise TreeAssemblyError(
            f'{subprocess.list2cmdline(command)} failed in {cwd}: {exc}'
            ) from exc


def main():
    """Assemble the TensErLEED test tree into the directory given."""
    parser = argparse.ArgumentParser(
        description='Assemble the TensErLEED tree for the tensorleed '
                    'test tier into TARGET_DIR (must not exist or be '
                    'empty).',
        )
    parser.add_argument('target_dir', help='directory to assemble into')
    args = parser.parse_args()
    try:
        assemble_tree(args.target_dir)
    except (OSError, TreeAssemblyError, tarfile.TarError) as exc:
        sys.exit(f'Error: {exc}')


if __name__ == '__main__':
    main()
