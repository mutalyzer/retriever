import json

from ..configuration import settings
from ..request import Http400, RequestErrors, request
from ..util import f_e


def fetch_json(feature_id, api_base, timeout=1):
    url = f"{api_base}/lookup/id/{feature_id}"
    params = {"feature": ["gene", "transcript", "cds"], "expand": 1}
    headers = {"Content-Type": "application/json"}
    try:
        response = request(url, params, headers, timeout=timeout)
    except RequestErrors as e:
        raise ConnectionError(f"(json) {str(e)}")
    except Http400 as e:
        response_json = e.response.json()
        if response_json and response_json.get("error") == "ID '{}' not found".format(
            feature_id
        ):
            raise NameError(f"(json) {str(e)}")
        else:
            raise e
    else:
        return response


def fetch_fasta(feature_id, api_base, timeout=1):
    url = f"{api_base}/sequence/id/{feature_id}"
    params = {"format": "fasta", "type": "genomic"}
    headers = {"Content-Type": "text/x-fasta"}

    try:
        response = request(url, params, headers, timeout=timeout)
    except RequestErrors as e:
        raise ConnectionError(f_e("gff3", e))
    except Http400 as e:
        response_json = e.response.json()
        if response_json and response_json.get("error") == "ID '{}' not found".format(
            feature_id
        ):
            raise NameError(f_e("fasta", e, response_json.get("error")))
        else:
            raise e
    else:
        return response


def fetch_gff3(feature_id, api_base, timeout=1):
    url = f"{api_base}/overlap/id/{feature_id}"
    params = {"feature": ["gene", "transcript", "cds", "exon"]}
    headers = {"Content-Type": "text/x-gff3"}

    try:
        response = request(url, params, headers, timeout=timeout)
    except RequestErrors as e:
        raise ConnectionError(f_e("gff3", e))
    except Http400 as e:
        response_json = e.response.json()
        if response_json and response_json.get("error") == "ID '{}' not found".format(
            feature_id
        ):
            raise NameError(f_e("gff3", e, response_json.get("error")))
        else:
            raise e
    else:
        return response


def _get_most_recent_version(reference_id, api_base, timeout=1):
    return int(_get_reference_information(reference_id, api_base, timeout)["version"])


def _get_reference_information(reference_id, api_base, timeout=1):
    url = f"{api_base}/lookup/id/{reference_id}"
    headers = {"Content-Type": "application/json"}
    return json.loads(request(url, headers=headers, timeout=timeout))


def _get_id_and_version(reference_id):
    r_id = None
    r_version = None
    if reference_id.startswith("ENS"):
        if (
            "." in reference_id
            and len(reference_id.split(".")) == 2
            and reference_id.split(".")[1].isdigit()
        ):
            r_id, r_version = reference_id.split(".")
            r_version = int(r_version)
        else:
            r_id = reference_id
    return r_id, r_version


def _in_grch37(r_id, r_version, r_info, timeout):
    api_base = settings.get("ENSEMBL_API_GRCH37")
    if r_info["species"] == "homo_sapiens" and int(r_info["version"]) > r_version:
        grch37_version = _get_most_recent_version(r_id, api_base, timeout)
        if grch37_version and grch37_version == r_version:
            return True
    return False


def fetch(reference_id, reference_type=None, timeout=1):
    api_base = settings.get("ENSEMBL_API")
    r_id, r_version = _get_id_and_version(reference_id)

    if r_id is None:
        raise NameError
    elif r_version is not None:
        r_info = _get_reference_information(r_id, api_base, timeout)
        if int(r_info["version"]) > r_version:
            if _in_grch37(r_id, r_version, r_info, timeout):
                api_base = settings.get("ENSEMBL_API_GRCH37")
            else:
                raise NameError

    if reference_type in [None, "gff3"]:
        return fetch_gff3(r_id, api_base, timeout), "gff3"
    elif reference_type == "fasta":
        return fetch_fasta(r_id, api_base, timeout), "fasta"
    elif reference_type == "json":
        return fetch_json(r_id, api_base, timeout), "json"
    elif reference_type == "genbank":
        return None, "genbank"

    raise ValueError(
        "Ensembl fetch does not support '{}' reference type.".format(reference_type)
    )
