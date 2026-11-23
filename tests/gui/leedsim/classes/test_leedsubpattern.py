"""Tests for module leedsubpattern of viperleed.gui.leedsim.classes."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-23'
__license__ = 'GPLv3+'

import matplotlib
from matplotlib import colors as mpl_colors
import numpy as np
import pytest
from pytest_cases import fixture, parametrize

from viperleed.gui.leedsim.classes.leedpattern import LEEDPattern
from viperleed.gui.leedsim.classes.leedsubpattern import LEEDSubpattern

# A (2x1) superstructure on a square lattice, i.e., with several
# symmetry-induced domains, which are colored using a colormap
_MULTI_DOMAIN = {
    'eMax': 100,
    'SUPERLATTICE': ((2, 0), (0, 1)),
    'surfBasis': ((2.0, 0.0), (0.0, 1.0)),
    'surfGroup': 'p1',
    'bulkGroup': 'p4m',
    }
_UNREGISTERED_CMAP = mpl_colors.ListedColormap(('red', 'green', 'blue'))


def _expected_colors(cmap, n_colors):
    """Return the colors that a LEEDSubpattern picks from cmap."""
    return cmap(np.linspace(0.1, 0.9, n_colors))


class TestLEEDSubpattern:
    """Tests for the LEEDSubpattern class."""

    @fixture(name='leed')
    def fixture_leed(self):
        """Return a LEEDPattern with several domains."""
        return LEEDPattern(_MULTI_DOMAIN)

    @fixture(name='make_subpattern')
    def factory_subpattern(self, leed):
        """Return a LEEDSubpattern for one of the domains of leed."""
        def _make(**kwargs):
            domain_id = leed.domains.domain_ids[0] + '0'
            beams = ((0.5, 0.0), (1.5, 0.0))
            return LEEDSubpattern(leed, beams, domain_id, **kwargs)
        return _make

    def test_cmap_default(self, leed):
        """Check the colors of the domains for the default colormap."""
        subpattern, *_ = leed.pattern['surf']
        colors = np.array(list(subpattern.leed_colors.values()))
        assert len(colors) > 1
        default = matplotlib.colormaps[LEEDSubpattern.defaults['cmap']]
        assert np.allclose(colors, _expected_colors(default, len(colors)))

    _invalid_cmap = {
        'unknown name': 'not-a-colormap',
        'unhashable': ['viridis'],
        }

    @parametrize(cmap=_invalid_cmap.values(), ids=_invalid_cmap)
    def test_cmap_invalid(self, cmap, make_subpattern):
        """Check complaints for an invalid colormap."""
        with pytest.raises(ValueError):
            make_subpattern(cmap=cmap)

    _valid_cmap = {  # cmap, expected colormap
        'name': ('viridis', matplotlib.colormaps['viridis']),
        'colormap': (_UNREGISTERED_CMAP, _UNREGISTERED_CMAP),
        }

    @parametrize('cmap,expect', _valid_cmap.values(), ids=_valid_cmap)
    def test_cmap_valid(self, cmap, expect, make_subpattern):
        """Check the colors of the domains for a given colormap."""
        subpattern = make_subpattern(cmap=cmap)
        colors = np.array(list(subpattern.leed_colors.values()))
        assert len(colors) > 1
        assert np.allclose(colors, _expected_colors(expect, len(colors)))
