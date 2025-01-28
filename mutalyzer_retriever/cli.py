"""
CLI entry point.
"""

import argparse
import json
import sys

from . import usage, version
from .related import get_cds_to_mrna, get_related
from .retriever import retrieve_model, retrieve_model_from_file, retrieve_raw
from .sources.ncbi_assemblies import annotations_summary, retrieve_assemblies


def _args_parser():
    """
    Command line argument parser.
    """
    parser = argparse.ArgumentParser(
        description=usage[0],
        epilog=usage[1],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("-v", action="version", version=version(parser.prog))

    parser.add_argument("--id", help="the reference id")

    parser.add_argument(
        "-s",
        "--source",
        help="retrieval source",
        choices=["ncbi", "ensembl", "ensembl_tark", "ensembl_rest", "lrg"],
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

    parser.add_argument(
        "-r", "--related", help="retrieve related reference ids", action="store_true"
    )

    parser.add_argument(
        "--mrna_id", help="retrieve the mrna_id from the cds_id", action="store_true"
    )

    parser.add_argument("--timeout", help="timeout", type=int)

    parser.add_argument("--indent", help="indentation spaces", default=None, type=int)

    parser.add_argument(
        "--sizeoff", help="do not consider file size", action="store_true"
    )

    parser.add_argument("--output", help="directory output path")

    parser.add_argument("--split", action="store_true")

    parser.add_argument("--only_annotations", action="store_true")

    subparsers = parser.add_subparsers(dest="command")

    parser_from_file = subparsers.add_parser(
        "from_file", help="parse files to get the model"
    )

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

    parser_ncbi_assemblies = subparsers.add_parser(
        "ncbi_assemblies", help="retrieve NCBI genomic annotations (including history)"
    )

    parser_ncbi_assemblies.add_argument(
        "--input", default="./ncbi_annotation_releases", help="input (downloaded) directory path"
    )
    parser_ncbi_assemblies.add_argument(
        "--output", default="./ncbi_annotation_models", help="output (models) directory path"
    )
    parser_ncbi_assemblies.add_argument(
        "--downloaded", help="already downloaded", action="store_true"
    )
    parser_ncbi_assemblies.add_argument(
        "--assembly_id_start", help="assembly id should start with"
    )
    parser_ncbi_assemblies.add_argument(
        "--ref_id_start", help="reference id should start with"
    )

    parser_ncbi_assemblies.add_argument(
        "--include_sequence", help="download also the sequence", action="store_true"
    )

    parser_assemblies_summary = subparsers.add_parser(
        "summary", help="gather references summary"
    )

    parser_assemblies_summary.add_argument("--directory", help="models directory path")
    parser_assemblies_summary.add_argument("--ref_id_start")

    return parser


def parse_args(args=None):
    return _args_parser().parse_args(args)


def _write_model(model, args):
    if args.split:
        if model.get("annotations"):
            with open(f"{args.output}/{args.id}.annotations", "w") as f:
                f.write(json.dumps(model["annotations"], indent=args.indent))
        if model.get("sequence"):
            with open(f"{args.output}/{args.id}.sequence", "w") as f:
                f.write(model["sequence"]["seq"])
    else:
        with open(f"{args.output}/{args.id}", "w") as f:
            f.write(json.dumps(model, indent=args.indent))


def _output_model(model, args):
    if args.output:
        _write_model(model, args)
    else:
        print(json.dumps(model, indent=args.indent))


def _from_file(args):
    model = retrieve_model_from_file(paths=args.paths, is_lrg=args.is_lrg)
    _output_model(model, args)


def _retrieve_assemblies(args):
    retrieve_assemblies(
        directory_input=args.input,
        directory_output=args.output,
        assembly_id_start=args.assembly_id_start,
        ref_id_start=args.ref_id_start,
        downloaded=args.downloaded,
        include_sequence=args.include_sequence,
    )


def _retrieve_model(args):
    model = retrieve_model(
        reference_id=args.id,
        reference_source=args.source,
        reference_type=args.type,
        model_type=args.model_type,
        size_off=args.sizeoff,
        timeout=args.timeout,
    )
    _output_model(model, args)


def _related(args):
    output = get_related(
        reference_id=args.id,
        timeout=args.timeout,
    )
    print(json.dumps(output, indent=args.indent))


def _cds_to_mrna(args):
    output = get_cds_to_mrna(
        cds_id=args.id,
        timeout=args.timeout,
    )
    print(json.dumps(output, indent=args.indent))


def _retrieve_raw(args):
    output = retrieve_raw(
        reference_id=args.id,
        reference_source=args.source,
        reference_type=args.type,
        size_off=args.sizeoff,
        timeout=args.timeout,
    )
    print(output[0])


def _assemblies_summary(args):
    annotations_summary(args.directory, args.ref_id_start)


def _endpoint(args):
    if args.command == "from_file":
        return _from_file
    elif args.command == "ncbi_assemblies":
        return _retrieve_assemblies
    elif args.command == "summary":
        return _assemblies_summary
    elif args.parse:
        return _retrieve_model
    elif args.related:
        return _related
    elif args.mrna_id:
        return _cds_to_mrna
    else:
        return _retrieve_raw


def main():
    """
    Main entry point.
    """
    args_parser = _args_parser()

    if len(sys.argv) == 1:
        args_parser.print_help(sys.stderr)
        sys.exit(1)

    args = args_parser.parse_args()
    _endpoint(args)(args)
