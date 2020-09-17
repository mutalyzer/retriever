"""
CLI entry point.
"""

import argparse
import json

from . import usage, version
from .retriever import retrieve_raw, retrieve_model, retrieve_model_from_file


def main():
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser(
        description=usage[0],
        epilog=usage[1],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser = argparse.ArgumentParser()

    parser.add_argument("-v", action="version", version=version(parser.prog))

    parser.add_argument("--reference_id", help="the reference id")

    parser.add_argument(
        "--sizeoff", help="do not consider file size", action="store_true"
    )

    parser.add_argument("--parse", help="parse reference content", action="store_true")

    parser.add_argument(
        "--indent", help="indentation spaces", default=None
    )

    parser.add_argument(
        "--source", help="retrieval source", choices=["ncbi", "ensembl", "lrg"]
    )

    parser.add_argument(
        "--type", help="reference type", choices=["gff3", "genbank", "json", "fasta"]
    )

    parser.add_argument(
        "--model_type",
        help="include the complete model or parts of it",
        choices=["all", "sequence", "annotations"],
        default="all",
    )

    parser.add_argument("-c", "--configuration", help="configuration file path")

    subparsers = parser.add_subparsers(help="parse files to get the model", dest="from_file")

    parser_from_file = subparsers.add_parser("from_file")

    parser_from_file.add_argument(
        "--paths",
        help="either gff3 and fasta paths or just one for lrg",
        nargs="+",
    )
    parser_from_file.add_argument(
        "--is_lrg",
        help="there is one file which is lrg",
        action="store_true",
        default=False,
    )

    args = parser.parse_args()

    if args.indent: args.indent = int(args.indent)

    if args.from_file:
        output = retrieve_model_from_file(paths=args.paths, is_lrg=args.is_lrg)
        print(json.dumps(output, indent=args.indent))
    elif args.parse:
        output = retrieve_model(
            reference_id=args.reference_id,
            reference_source=args.source,
            reference_type=args.type,
            model_type=args.model_type,
            size_off=args.sizeoff,
        )
        print(json.dumps(output, indent=args.indent))
    else:
        output = retrieve_raw(
            reference_id=args.reference_id,
            reference_source=args.source,
            reference_type=args.type,
            size_off=args.sizeoff,
            configuration_path=args.configuration,
        )
        print(output[0])
