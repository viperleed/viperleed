"""Supply the license in the form the installed setuptools supports.

PEP 639 license metadata (a SPDX license expression) needs
setuptools>=77, which does not run on Python<=3.8. The older form
(a TOML table) is deprecated in recent setuptools versions. Keep
'license' dynamic in pyproject.toml and fill it in here instead.
Everything else is configured in pyproject.toml. See also #364.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-13'
__license__ = 'GPLv3+'

import setuptools

_LICENSE = 'GPL-3.0-or-later'
_PEP639_SETUPTOOLS = 77

if int(setuptools.__version__.split('.')[0]) >= _PEP639_SETUPTOOLS:
    setuptools.setup(license_expression=_LICENSE)
else:
    setuptools.setup(license=_LICENSE)
