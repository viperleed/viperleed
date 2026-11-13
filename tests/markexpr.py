"""Helper for the plugins that exclude a marker from the selection."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-29'
__license__ = 'GPLv3+'


def exclude_marker(config, marker):
    """Add 'not `marker`' to the marker expression of `config`.

    The existing -m expression, if any, is parenthesized, so that
    an 'or' in it only applies to the tests that are not excluded.

    Parameters
    ----------
    config : pytest.Config
        The configuration whose option `markexpr` is modified.
    marker : str
        Name of the marker of the tests to exclude.

    Returns
    -------
    None.
    """
    expression = config.option.markexpr
    config.option.markexpr = (
        f'({expression}) and not {marker}' if expression
        else f'not {marker}'
        )
