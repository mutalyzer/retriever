"""
CLI entry point.
"""

import argparse
import json

from . import usage, version
from .retriever import retrieve_model, retrieve_model_from_file, retrieve_raw


def _arg_parser():
    """
    Command line argument parsing.
    """
    parser = argparse.ArgumentParser(
        description=usage[0],
        epilog=usage[1],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("-v", action="version", version=version(parser.prog))

    parser.add_argument("--id", help="the reference id")

    parser.add_argument(
        "-s", "--source", help="retrieval source", choices=["ncbi", "ensembl", "lrg"]
    )

    parser.add_argument(
        "-t",
        "--type",
        help="reference type",
        choices=["gff3", "genbank", "json", "fasta"],
    )

    parser.add_argument(
        "-p", "--parse", help="parse reference content", action="store_true"
    )

    parser.add_argument(
        "-m",
        "--model_type",
        help="include the complete model or parts of it",
        choices=["all", "sequence", "annotations"],
        default="all",
    )

    parser.add_argument("--timeout", help="timeout", type=int)

    parser.add_argument("--indent", help="indentation spaces", default=None)

    parser.add_argument(
        "--sizeoff", help="do not consider file size", action="store_true"
    )

    parser.add_argument("-c", "--configuration", help="configuration file path")

    subparsers = parser.add_subparsers(dest="from_file")

    parser_from_file = subparsers.add_parser("from_file", help="parse files to get the model")

    parser_from_file.add_argument(
        "--paths",
        help="both gff3 and fasta paths or just an lrg",
        nargs="+",
    )
    parser_from_file.add_argument(
        "--is_lrg",
        help="there is one file which is lrg",
        action="store_true",
        default=False,
    )
    return parser


def main():
    """
    Main entry point.
    """
    parser = _arg_parser()

    try:
        args = parser.parse_args()
    except IOError as error:
        parser.error(error)

    if args.indent:
        args.indent = int(args.indent)

    if args.from_file:
        output = retrieve_model_from_file(paths=args.paths, is_lrg=args.is_lrg)
        print(json.dumps(output, indent=args.indent))
    elif args.parse:
        output = retrieve_model(
            reference_id=args.id,
            reference_source=args.source,
            reference_type=args.type,
            model_type=args.model_type,
            size_off=args.sizeoff,
            timeout=args.timeout,
        )
        print(json.dumps(output, indent=args.indent))
    else:
        output = retrieve_raw(
            reference_id=args.id,
            reference_source=args.source,
            reference_type=args.type,
            size_off=args.sizeoff,
            configuration_path=args.configuration,
            timeout=args.timeout,
        )
        print(output[0])
