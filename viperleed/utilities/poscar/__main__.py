#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED poscar utilities.
"""
import argparse

from importlib import import_module
from viperleed.utilities import poscar
from viperleed.utilities.poscar import POSCAR_UTILITIES, add_verbose_option

def add_poscar_parser_arguments(parser):
    subparsers = parser.add_subparsers()

    for utility in POSCAR_UTILITIES:
        _util_parser = subparsers.add_parser(
            utility,
            help=f"call poscar utility {utility}"
        )
        import_module(f"viperleed.utilities.poscar.{utility}")
        sub_module = getattr(poscar, utility)
        add_verbose_option(_util_parser)
        sub_module.add_cli_parser_arguments(_util_parser)
        _util_parser.set_defaults(func=sub_module.main)


def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        add_poscar_parser_arguments(parser)
        args = parser.parse_args()

    args.func(args)

if __name__ == "__main__":
    main()
