from . import configuration, parser
from .sources import ensembl, lrg, ncbi


class NoReferenceRetrieved(Exception):
    pass


class NoReferenceError(Exception):
    def __init__(self, status):
        self.no_reference = []
        message = ""
        for source in status.keys():
            source_errors = []
            for error in status[source]["errors"]:
                if isinstance(error, ValueError):
                    source_errors.append("ValueError")
                elif isinstance(error, NameError):
                    source_errors.append("NameError")
                elif isinstance(error, ConnectionError):
                    source_errors.append("ConnectionError")
            if source_errors:
                message += "\n{}: {}.".format(source, ", ".join(source_errors))

        self.message = message

    def __str__(self):
        return self.message


def raise_error(status):
    no_reference_retrieved = []
    for source in status.keys():
        if len(status[source]["errors"]) == 1 and isinstance(
            status[source]["errors"][0], NameError
        ):
            no_reference_retrieved.append(True)
        else:
            no_reference_retrieved.append(False)
    if all(no_reference_retrieved):
        raise NoReferenceRetrieved

    raise NoReferenceError(status)


def fetch_unknown_source(reference_id, reference_type, size_off=True):

    status = {"lrg": {"errors": []}, "ncbi": {"errors": []}, "ensembl": {"errors": []}}

    # LRG
    if reference_type in [None, "lrg"]:
        try:
            reference_content = lrg.fetch_lrg(reference_id)
        except (NameError, ConnectionError) as e:
            status["lrg"]["errors"].append(e)
        else:
            return reference_content, "lrg", "lrg"
    else:
        status["lrg"]["errors"].append(
            ValueError(
                "Lrg fetch does not support '{}' reference type.".format(reference_type)
            )
        )

    # NCBI
    try:
        reference_content, reference_type = ncbi.fetch(
            reference_id, reference_type, size_off
        )
    except (NameError, ConnectionError, ValueError) as e:
        status["ncbi"]["errors"].append(e)
    else:
        return reference_content, reference_type, "ncbi"

    # Ensembl
    try:
        reference_content, reference_type = ensembl.fetch(reference_id, reference_type)
    except (NameError, ConnectionError, ValueError) as e:
        status["ensembl"]["errors"].append(e)
    else:
        return reference_content, reference_type, "ensembl"

    raise_error(status)


def retrieve_raw(
    reference_id,
    reference_source=None,
    reference_type=None,
    size_off=True,
    configuration_path=None,
):
    configuration.settings = configuration.setup_settings(configuration_path)

    reference_content = None

    if reference_source is None:
        reference_content, reference_type, reference_source = fetch_unknown_source(
            reference_id, reference_type, size_off
        )
    elif reference_source == "ncbi":
        reference_content, reference_type = ncbi.fetch(reference_id, reference_type)
    elif reference_source == "ensembl":
        reference_content, reference_type = ensembl.fetch(reference_id, reference_type)
    elif reference_source == "lrg":
        reference_content = lrg.fetch_lrg(reference_id)
        if reference_content:
            reference_type = "lrg"

    return reference_content, reference_type, reference_source


def retrieve_model(
    reference_id,
    reference_source=None,
    reference_type=None,
    size_off=True,
    model_type="all",
    configuration_path=None,
):
    configuration.settings = configuration.setup_settings(configuration_path)

    reference_content, reference_type, reference_source = retrieve_raw(
        reference_id, reference_source, reference_type, size_off
    )

    if reference_type == "lrg":
        model = parser.parse(reference_content, reference_type, reference_source)
        if model_type == "all":
            return model
        elif model_type == "sequence":
            return model["sequence"]
        elif model_type == "annotations":
            return model["annotations"]
    elif reference_type == "gff3":
        if model_type == "all":
            fasta = retrieve_raw(reference_id, reference_source, "fasta", size_off)
            return {
                "annotations": parser.parse(
                    reference_content, reference_type, reference_source
                ),
                "sequence": parser.parse(fasta[0], "fasta"),
            }
        elif model_type == "sequence":
            fasta = retrieve_raw(reference_id, "fasta", size_off)
            return {"sequence": parser.parse(fasta, "fasta")}
        elif model_type == "annotations":
            return parser.parse(
                reference_content, reference_source, "fasta", reference_source
            )
    elif reference_type == "fasta":
        return {
            "sequence": parser.parse(reference_content, "fasta"),
        }


def retrieve_model_from_file(paths=[], is_lrg=False):
    if is_lrg:
        with open(paths[0]) as f:
            content = f.read()
            model = parser.parse(content, "lrg")
            return model

    gff3 = paths[0]
    fasta = paths[1]

    model = {}
    with open(gff3) as f:
        annotations = f.read()
        model["annotations"] = parser.parse(annotations, "gff3")

    with open(fasta) as f:
        sequence = f.read()
        model["sequence"] = parser.parse(sequence, "fasta")

    return model
