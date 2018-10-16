"""
CLI entry point.
"""

import argparse

from . import usage, version
from .retriever import retrieve


def main():
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser(
        description=usage[0], epilog=usage[1],
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', action='version', version=version(parser.prog))

    parser.add_argument("reference", help="the reference id")

    parser.add_argument("--sizeoff", help="do not consider file size",
                        action="store_true")

    args = parser.parse_args()

    print(retrieve(args.reference, not args.sizeoff))