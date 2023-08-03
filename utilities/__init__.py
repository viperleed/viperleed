import argparse

def default_cli_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v", "--verbose",
        help=("increase output verbosity"),
        action="store_true"
    )
    return parser