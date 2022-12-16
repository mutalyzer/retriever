import json
from copy import deepcopy
from functools import lru_cache
from pathlib import Path

import requests

from . import parser
from .configuration import cache_add, cache_dir, cache_url, lru_cache_maxsize
from .sources import ensembl, lrg, ncbi


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
    except Exception as e:
        raise e
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


@lru_cache(maxsize=lru_cache_maxsize())
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
            annotations = parser.parse(
                reference_content, reference_type, reference_source
            )
            fasta = retrieve_raw(
                reference_id, reference_source, "fasta", size_off, timeout=timeout
            )
            return {
                "annotations": annotations,
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


@lru_cache(maxsize=lru_cache_maxsize())
def get_annotations_from_file_cache(r_id):
    cache_path = cache_dir()
    if cache_path and (Path(cache_path) / (r_id + ".annotations")).is_file():
        with open(Path(cache_path) / (r_id + ".annotations")) as json_file:
            return json.load(json_file)


@lru_cache(maxsize=lru_cache_maxsize())
def get_sequence_from_file_cache(r_id):
    cache_path = cache_dir()
    if cache_path and (Path(cache_path) / (r_id + ".sequence")).is_file():
        with open(Path(cache_path) / (r_id + ".sequence")) as seq_file:
            return {"seq": seq_file.read()}


def get_from_api_cache(r_id, s_id):
    api_url = cache_url()
    if api_url:
        url = api_url + "/reference/" + r_id
        if s_id:
            url += f"?selector_id={s_id}"
        try:
            annotations = requests.get(url).text
            annotations = json.loads(annotations)
            sequence = get_sequence_from_file_cache(r_id)
        except Exception:
            return

        if annotations and sequence:
            return {
                "annotations": annotations,
                "sequence": get_sequence_from_file_cache(r_id),
            }


def get_from_file_cache(r_id):
    annotations = get_annotations_from_file_cache(r_id)
    sequence = get_sequence_from_file_cache(r_id)
    if annotations and sequence:
        return {"annotations": annotations, "sequence": sequence}


def get_overlap_models(r_id, l_min, l_max):
    api_url = cache_url()
    if api_url:
        url = f"{api_url}/overlap/{r_id}?min={l_min}&max={l_max}"
        try:
            # print("- get overlap models from api cache")
            annotations = requests.get(url).text
            annotations = json.loads(annotations)
        except Exception:
            return
        return annotations
    model = get_from_file_cache(r_id)
    if model:
        return model


def get_reference_model(r_id, s_id=None):

    model = get_from_api_cache(r_id, s_id)
    if model:
        model["annotations"]["source"] = "api_cache"
        return model
    model = get_from_file_cache(r_id)
    if model:
        return model
    model = retrieve_model(r_id, timeout=10)

    cache_path = cache_dir()
    if cache_add() and cache_path:
        if (
            model.get("annotations")
            and model.get("sequence")
            and model["sequence"].get("seq")
        ):
            with open(Path(cache_path) / (r_id + ".annotations"), "w") as f:
                f.write(json.dumps(model["annotations"]))
            with open(Path(cache_path) / (r_id + ".sequence"), "w") as f:
                f.write(model["sequence"]["seq"])
    return model


def get_reference_model_segmented(
    reference_id, feature_id=None, siblings=False, ancestors=True, descendants=True
):
    def _get_from_api_cache():
        api_url = cache_url()
        if api_url:
            if feature_id is not None:
                args = f"&feature_ud={descendants}"
            else:
                args = ""
            args += (
                f"&siblings={siblings}&ancestors={ancestors}&descendants={descendants}"
            )
            url = f"{api_url}/reference_model_segmented/{reference_id}{args}"
            try:
                annotations = requests.get(url).text
                return json.loads(annotations)
            except Exception:
                return

    model = _get_from_api_cache()
    if model:
        return model
    model = get_from_file_cache(reference_id)
    if model is None:
        model = retrieve_model(reference_id, timeout=10)

    if feature_id is not None:
        return extract_feature_model(
            model["annotations"],
            feature_id,
            siblings,
            ancestors,
            descendants,
        )[0]
    return model


def get_chromosome_from_selector(assembly_id, selector_id):
    api_url = cache_url()
    if api_url:
        url = f"{api_url}/chromosome_from_selector/{selector_id}?assembly_id={assembly_id}"
        try:
            response = requests.get(url).text
            return json.loads(response)["id"]
        except Exception:
            return


def extract_feature_model(
    feature, feature_id, siblings=False, ancestors=True, descendants=True
):
    output_model = None
    just_found = False
    if feature.get("id") is not None and feature_id == feature["id"]:
        output_model = deepcopy(feature)
        if not descendants:
            if output_model.get("features"):
                output_model.pop("features")
            return output_model, True, True
        return output_model, True, False
    elif feature.get("features"):
        for f in feature["features"]:
            output_model, just_found, propagate = extract_feature_model(
                f, feature_id, siblings, ancestors, descendants
            )
            if output_model:
                break
        if output_model and just_found:
            if siblings:
                output_model = deepcopy(feature["features"])
            if not ancestors:
                return output_model, False, True
        elif propagate:
            return output_model, False, True
    if output_model is not None:
        if isinstance(output_model, dict):
            output_model = [output_model]
        return (
            {
                **{
                    k: deepcopy(feature[k])
                    for k in list(set(feature.keys()) - {"features"})
                },
                **{"features": output_model},
            },
            False,
            False,
        )
    return None, False, False
