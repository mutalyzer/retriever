from .parsers import fasta, gff3, json_ensembl, lrg


def _get_reference_type(content):
    if content.startswith("<?xml version="):
        return "lrg"
    elif content.startswith("LOCUS"):
        return "genbank_ncbi"


def parse(reference_content, reference_type=None, reference_source=None):

    if reference_type is None:
        reference_type = _get_reference_type(reference_content)

    if reference_type == "lrg":
        model = lrg.parse(reference_content)
    elif reference_type == "gff3":
        model = gff3.parse(reference_content, reference_source)
    elif reference_type == "fasta":
        model = fasta.parse(reference_content)
    elif reference_type == "json":
        if reference_content["results"]:
            model = json_ensembl.parse(reference_content["results"])
    else:
        return None

    return model
