"""Script assemble_tensorleed_tree.

Assembles the TensErLEED tree that the tensorleed test tier runs
against (pytest -m tensorleed, with VIPERLEED_TENSORLEED pointing to
the tree). Used by .github/workflows/tensorleed_tests.yml and by
developers. Needs Internet access, a C and a Fortran compiler (by
default gcc and gfortran; see --help for the Intel compilers), and
make (Linux, macOS) or make.bat (Windows, shipped with the EEASiSSS
sources).

Usage: python assemble_tensorleed_tree.py [options] TARGET_DIR

The tree is made from release tags of viperleed/viperleed-tensorleed,
downloaded as archives of the tagged commits (no git needed). The
module constants say which: TAG_COMMITS pins the commit of each tag,
TREE_CONTENTS lists the folders copied from each tag, and the
eeasisss executable is compiled from EEASISSS_TAG. In LEGACY_LIB, the
shipped random_.o and MPIrandom_.o are replaced by a random_.o
compiled on this machine, as the installation instructions of that
TensErLEED version require. It is compiled from RANDOM_SOURCE, the
cross-platform random_.c of the ViPErLEED documentation, which also
replaces the shipped random_.c in the tree. At the end, the script
warns if viperleed-tensorleed has release tags newer than the newest
pinned one of their kind (TensErLEED or EEASiSSS), as TAG_COMMITS may
need an update then. If the tags cannot be listed, it says so and
carries on.

TARGET_DIR must not exist or be empty: the script refuses to write
into a non-empty directory. Downloads and builds are made in a
temporary subfolder of TARGET_DIR, which is removed at the end. If a
step fails, everything the script wrote is removed, so that it can
be run again on the same TARGET_DIR.

Tested on Linux only, and in CI with gcc and gfortran only.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-17'
__license__ = 'GPLv3+'

import argparse
import contextlib
import json
import ntpath
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import urllib.request

ARCHIVE_URL = ('https://github.com/viperleed/viperleed-tensorleed/'
               'archive/{commit}.tar.gz')

EEASISSS_TAG = 'EEASiSSS-v1.0.2'
TAGS_URL = ('https://api.github.com/repos/viperleed/viperleed-tensorleed/'
            'tags?per_page=100')
# A release tag: a kind ('' for TensErLEED, 'EEASiSSS-') and a version
VERSION_TAG = re.compile(r'(?P<kind>.*?)v(?P<version>\d+(?:\.\d+)*)')

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
    (EEASISSS_TAG, 'atom_density_files', 'atom_density_files'),
    (EEASISSS_TAG, 'seSernelius', 'seSernelius'),
    )

# Compiler families, named as the make targets of eeasisss, and the
# C compiler each uses for random_.c by default
DEFAULT_C_COMPILERS = {
    'gcc': 'gcc',
    'intel': 'icx',
    }
# Compiler families whose C compiler takes MSVC-style options on Windows
MSVC_STYLE_ON_WINDOWS = {'intel'}
LEGACY_LIB = Path('TensErLEED-v1.61', 'lib')
# The shipped random_.c compiles on UNIX only. This one is used on
# all platforms, as the installation instructions suggest for Windows.
RANDOM_SOURCE = (Path(__file__).resolve().parents[2]
                 / 'doc/_static/tenserleed_windows/random_.c')
WORK_DIR_NAME = '_assemble_tmp'


class TreeAssemblyError(Exception):
    """The TensErLEED tree could not be assembled."""


def assemble_tree(target, compilers='gcc', c_compiler=None,
                  fortran_compiler=None):
    """Assemble the TensErLEED test tree into a directory.

    Parameters
    ----------
    target : str or Path
        The directory to assemble the tree into. It is created if it
        does not exist.
    compilers : {'gcc', 'intel'}, optional
        The compiler family: the make target of eeasisss, and the
        command-line syntax of the C compiler on Windows. Default is
        'gcc'.
    c_compiler : str, optional
        The C compiler that compiles random_.c. Default is the one of
        DEFAULT_C_COMPILERS for compilers.
    fortran_compiler : str, optional
        The Fortran compiler that compiles eeasisss. Default is the
        one of the make target of eeasisss (gfortran or ifort). Not
        supported on Windows, where make.bat fixes it.

    Raises
    ------
    TreeAssemblyError
        If target is not an empty directory, if fortran_compiler is
        given on Windows, if an archive contains an unsafe member or
        more than one root folder, or if a compilation step fails.
    OSError
        If a download fails (urllib.error.URLError is an OSError) or
        a file cannot be written.
    tarfile.TarError
        If a downloaded archive cannot be read.
    """
    random_command = _random_command(compilers, c_compiler)
    make_command = _make_command(compilers, fortran_compiler)
    target = Path(target).resolve()
    if target.exists() and (not target.is_dir() or any(target.iterdir())):
        raise TreeAssemblyError(
            f'{target} exists and is not an empty directory. Remove it '
            'or choose another target.'
            )
    with _removed_on_failure(target):
        _assemble_tree(target, random_command, make_command)
    print(f'TensErLEED tree assembled in {target}')
    _warn_about_newer_tags()


def _assemble_tree(target, random_command, make_command):
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
    _compile_random(target / LEGACY_LIB, random_command)
    _compile_eeasisss(sources[EEASISSS_TAG], target, make_command)
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


def _compile_eeasisss(src, target, make_command):
    """Compile eeasisss in folder src and move it to target."""
    print(f'Compiling eeasisss from {EEASISSS_TAG}')
    _run(make_command, src)
    for name in ('eeasisss', 'eeasisss.exe'):
        executable = src / name
        if executable.is_file():
            shutil.copy2(executable, target / name)
            return
    raise TreeAssemblyError(f'Compiling eeasisss produced no executable '
                            f'in {src}')


def _compile_random(lib, command):
    """Compile RANDOM_SOURCE in place of the shipped randomizer."""
    print(f'Compiling random_.o in {lib} from {RANDOM_SOURCE}')
    for name in ('random_.o', 'MPIrandom_.o'):
        (lib / name).unlink()
    shutil.copyfile(RANDOM_SOURCE, lib / 'random_.c')
    _run(command, lib)


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


def _list_tags():
    """Return the names of the tags of viperleed-tensorleed."""
    with urllib.request.urlopen(TAGS_URL, timeout=30) as response:
        return [tag['name'] for tag in json.load(response)]


def _make_command(compilers, fortran_compiler):
    """Return the command that compiles eeasisss."""
    if _on_windows():
        if fortran_compiler:
            msg = ('A Fortran compiler cannot be chosen on Windows: '
                   'make.bat of eeasisss always uses its own for target '
                   f'{compilers}')
            raise TreeAssemblyError(msg)
        return ['cmd', '/c', 'make.bat', compilers]
    command = ['make', compilers]
    if fortran_compiler:  # Variables gcc_fcomp, intel_fcomp
        command.append(f'{compilers}_fcomp={fortran_compiler}')
    return command


def _make_parser():
    """Return the parser of the command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Assemble the TensErLEED tree for the tensorleed '
                    'test tier into TARGET_DIR (must not exist or be '
                    'empty). Only gcc and gfortran are tested in CI.',
        )
    parser.add_argument('target_dir', help='directory to assemble into')
    parser.add_argument(
        '--compilers',
        choices=sorted(DEFAULT_C_COMPILERS),
        default='gcc',
        help='compiler family: the make target of eeasisss (default: '
             '%(default)s)',
        )
    parser.add_argument(
        '--cc',
        help='C compiler for random_.c (default: gcc for gcc, icx for '
             'intel; e.g., icc for older Intel oneAPI)',
        )
    parser.add_argument(
        '--fc',
        help='Fortran compiler for eeasisss (default: the one of its '
             'make target, gfortran or ifort; e.g., ifx). Not on '
             'Windows, where make.bat fixes it',
        )
    return parser


def _newer_tags(tags):
    """Return the tags newer than the newest pinned of their kind."""
    newest = {}
    for tag in TAG_COMMITS:
        kind, version = _parse_tag(tag)
        newest[kind] = max(version, newest.get(kind, version))
    newer = []
    for tag in tags:
        kind, version = _parse_tag(tag)
        if kind in newest and version > newest[kind]:
            newer.append(tag)
    return sorted(newer)


def _on_windows():
    """Return whether this is Windows."""
    return sys.platform.startswith('win')


def _parse_tag(tag):
    """Return the kind and the version tuple of tag, or (None, ())."""
    match = VERSION_TAG.fullmatch(tag)
    if not match:
        return None, ()
    version = tuple(int(part) for part in match['version'].split('.'))
    return match['kind'], version


def _random_command(compilers, c_compiler):
    """Return the command that compiles random_.c to random_.o."""
    c_compiler = c_compiler or DEFAULT_C_COMPILERS[compilers]
    if compilers in MSVC_STYLE_ON_WINDOWS and _on_windows():
        return [c_compiler, '/c', 'random_.c', '/Fo:random_.o', '/O3']
    return [c_compiler, '-O3', '-c', 'random_.c', '-o', 'random_.o']


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


def _warn_about_newer_tags():
    """Print a warning if newer tags exist than those in TAG_COMMITS."""
    try:
        tags = _list_tags()
    except (OSError, ValueError, KeyError, TypeError) as exc:
        # No connection (URLError is an OSError), or not the expected
        # JSON. Never fail for this: the tree is assembled already
        print('Note: could not check viperleed-tensorleed for newer '
              f'release tags: {exc!r}', file=sys.stderr)
        return
    newer = ', '.join(_newer_tags(tags))
    if newer:
        print('Warning: viperleed-tensorleed has release tags newer than '
              f'those in TAG_COMMITS of {Path(__file__).name}: {newer}. '
              'Consider adding them.',
              file=sys.stderr)


def main(argv=None):
    """Assemble the TensErLEED test tree into the directory given.

    Parameters
    ----------
    argv : list of str, optional
        The command-line arguments. Default is sys.argv[1:].
    """
    args = _make_parser().parse_args(argv)
    try:
        assemble_tree(args.target_dir,
                      compilers=args.compilers,
                      c_compiler=args.cc,
                      fortran_compiler=args.fc)
    except (OSError, TreeAssemblyError, tarfile.TarError) as exc:
        sys.exit(f'Error: {exc}')


if __name__ == '__main__':
    main()
