from . import parser
from .sources import ensembl, lrg, ncbi
from .configuration import cache_dir, cache_url
import json
from pathlib import Path
import requests
from functools import lru_cache


class NoReferenceRetrieved(Exception):
    pass


class NoReferenceError(Exception):
    def __init__(self, status, uncertain_sources):
        self.uncertain_sources = uncertain_sources
        message = ""
        if uncertain_sources is not []:
            message = f"\n\nUncertain sources: {', '.join(uncertain_sources)}\n"

        for source in status.keys():
            source_errors = []
            message += f"\n{source}:"
            for error in status[source]["errors"]:
                if isinstance(error, ValueError):
                    detail = {"type": "ValueError", "details": str(error)}
                elif isinstance(error, NameError):
                    detail = {"type": "NameError", "details": str(error)}
                elif isinstance(error, ConnectionError):
                    detail = {"type": "ConnectionError", "details": str(error)}
                else:
                    detail = {"type": "Unknown", "details": str(error)}
                source_errors.append(detail)
                message += f"\n {detail['type']}: {detail['details']}"

        self.message = message

    def __str__(self):
        return self.message


def _raise_error(status):
    uncertain_sources = []
    for source in status.keys():
        if not (
            len(status[source]["errors"]) == 1
            and isinstance(status[source]["errors"][0], NameError)
        ):
            uncertain_sources.append(source)
    if uncertain_sources is []:
        raise NoReferenceRetrieved
    raise NoReferenceError(status, uncertain_sources)


def _fetch_unknown_source(reference_id, reference_type, size_off=True, timeout=1):

    status = {"lrg": {"errors": []}, "ncbi": {"errors": []}, "ensembl": {"errors": []}}

    # LRG
    if reference_type in [None, "lrg"]:
        try:
            reference_content = lrg.fetch_lrg(reference_id, timeout=timeout)
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
            reference_id, reference_type, size_off, timeout
        )
    except (NameError, ConnectionError, ValueError) as e:
        status["ncbi"]["errors"].append(e)
    else:
        return reference_content, reference_type, "ncbi"

    # Ensembl
    try:
        reference_content, reference_type = ensembl.fetch(
            reference_id, reference_type, timeout
        )
    except (NameError, ConnectionError, ValueError) as e:
        status["ensembl"]["errors"].append(e)
    else:
        return reference_content, reference_type, "ensembl"

    _raise_error(status)


def retrieve_raw(
    reference_id,
    reference_source=None,
    reference_type=None,
    size_off=True,
    timeout=1,
):
    """
    Retrieve a reference based on the provided id.

    :arg str reference_id: The id of the reference to retrieve.
    :arg str reference_source: A dedicated retrieval source.
    :arg str reference_type: A dedicated retrieval type.
    :arg bool size_off: Download large files.
    :arg float timeout: Timeout.
    :returns: Reference content.
    :rtype: str
    """
    reference_content = None

    if reference_source is None:
        reference_content, reference_type, reference_source = _fetch_unknown_source(
            reference_id, reference_type, size_off, timeout
        )
    elif reference_source == "ncbi":
        reference_content, reference_type = ncbi.fetch(
            reference_id, reference_type, timeout
        )
    elif reference_source == "ensembl":
        reference_content, reference_type = ensembl.fetch(
            reference_id, reference_type, timeout
        )
    elif reference_source == "lrg":
        reference_content = lrg.fetch_lrg(reference_id, timeout=timeout)
        if reference_content:
            reference_type = "lrg"

    return reference_content, reference_type, reference_source


def retrieve_model(
    reference_id,
    reference_source=None,
    reference_type=None,
    size_off=True,
    model_type="all",
    timeout=1,
):
    """
    Obtain the model of the provided reference id.

    :arg str reference_id: The id of the reference to retrieve.
    :arg str reference_source: A dedicated retrieval source.
    :arg str reference_type: A dedicated retrieval type.
    :arg bool size_off: Download large files.
    :arg float timeout: Timeout.
    :returns: Reference model.
    :rtype: dict
    """
    reference_content, reference_type, reference_source = retrieve_raw(
        reference_id, reference_source, reference_type, size_off, timeout=timeout
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
            fasta = retrieve_raw(
                reference_id, reference_source, "fasta", size_off, timeout=timeout
            )
            return {
                "annotations": parser.parse(
                    reference_content, reference_type, reference_source
                ),
                "sequence": parser.parse(fasta[0], "fasta"),
            }
        elif model_type == "sequence":
            fasta = retrieve_raw(reference_id, "fasta", size_off, timeout=timeout)
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
    """

    :arg list paths: Path towards the gff3, fasta, or lrg files.
    :arg bool is_lrg: If there is only one file path of an lrg.
    :returns: Reference model.
    :rtype: dict
    """
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


@lru_cache(maxsize=16)
def _get_sequence(file_path):
    with open(file_path) as f:
        return f.read()


def get_from_api_cache(r_id, s_id):
    api_url = cache_url()
    cache_path = cache_dir()
    if api_url:
        url = api_url + "/reference/" + r_id
        if s_id:
            url += f"?selector_id={s_id}"
        try:
            annotations = requests.get(url).text
        except Exception:
            return
        annotations = json.loads(annotations)

        file_path = Path(cache_path) / (r_id + ".sequence")
        if cache_path:
            if file_path.is_file():
                return {"annotations": annotations, "sequence": {"seq": _get_sequence(file_path)}}


def get_from_file_cache(r_id):
    cache_path = cache_dir()
    if cache_path and (Path(cache_path) / r_id).is_file():
        with open(Path(cache_path) / r_id) as json_file:
            return json.load(json_file)


def get_overlap_models(r_id, l_min, l_max):
    api_url = cache_url()
    if api_url:
        url = f"{api_url}/overlap/{r_id}?min={l_min}&max={l_max}"
        try:
            annotations = requests.get(url).text
        except Exception:
            return
        annotations = json.loads(annotations)
        return annotations


def get_reference_model(r_id, s_id=None):
    model = get_from_api_cache(r_id, s_id)
    if model:
        return model
    model = get_from_file_cache(r_id)
    if model:
        return model
    return retrieve_model(r_id, timeout=10)
