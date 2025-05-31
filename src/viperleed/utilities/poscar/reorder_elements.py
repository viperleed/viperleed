"""ViPErLEED utility: Reorder elements."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

from viperleed.calc.lib import periodic_table
from viperleed.utilities.poscar.base import _PoscarStreamCLI


class ReorderElementsCLI(_PoscarStreamCLI, cli_name='reorder_elements'):
    """Change the sort order of the elements in a POSCAR."""

    long_name = 'reorder elements'

    def add_parser_arguments(self, parser):
        """Add ordering flags as arguments to parser."""
        super().add_parser_arguments(parser)
        order = parser.add_mutually_exclusive_group()
        order.add_argument(
            '--ascending',
            help=('sort elements by atomic number in ascending order. '
                  'This is the default if no argument is given.'),
            action='store_const',
            dest='ordering',
            const=order_ascending,
            )
        order.add_argument(
            '-d', '--descending',
            help='sort elements by atomic number in descending order',
            action='store_const',
            dest='ordering',
            const=order_descending,
            )
        order.add_argument(
            '-a', '--alphabetical',
            help='sort elements alphabetically',
            action='store_const',
            dest='ordering',
            const=order_alphabetical,
            )
        order.add_argument(
            '-c', '--custom',
            help=('Sort elements by custom order. Provide a comma-separated '
                  'list of elements in the desired order.'),
            type=str,
            )

    def parse_cli_args(self, args):
        """Set a default if no argument is given."""
        parsed_args = super().parse_cli_args(args)
        if parsed_args.ordering is None and parsed_args.custom is None:
            parsed_args.ordering = order_ascending
        return parsed_args

    def process_slab(self, slab, args):
        """Reorder elements in slab according to args and return slab."""
        try:
            ordering = args.ordering or make_custom_ordering(args.custom, slab)
        except KeyError as exc:
            self.parser.error(str(exc))

        ordered_elements = sorted(slab.elements, key=ordering)
        slab.n_per_elem = {el: slab.n_per_elem[el] for el in ordered_elements}
        slab.sort_by_z()
        return slab


def make_custom_ordering(custom_order, slab):
    """Return a callable that sorts items according to custom_order."""
    custom_order_map = {el.strip().capitalize(): i
                        for i, el in enumerate(custom_order.split(','))}
    missing_elements = set(slab.elements) - set(custom_order_map)
    if missing_elements:
        # Prevents KeyError later on when accessing the slab elements
        raise KeyError('Missing elements in custom order: '
                       + ','.join(sorted(missing_elements)))

    def order_custom(element_name):
        return custom_order_map[element_name]

    return order_custom


def order_alphabetical(element_name):
    """Return a key for sorting element_name alphabetically."""
    return element_name.lower()


def order_ascending(element_name):
    """Return a key for sorting element_name by ascending atomic number."""
    return periodic_table.get_atomic_number(element_name)


def order_descending(element_name):
    """Return a key for sorting element_name by descending atomic number."""
    return -order_ascending(element_name)


if __name__ == '__main__':
    ReorderElementsCLI.run_as_script()
