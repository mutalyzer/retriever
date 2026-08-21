import json

import requests

from ..configuration import DEFAULT_TIMEOUT, settings
from ..request import Http400, RequestErrors, request
from ..util import f_e


def fetch_fasta(feature_id, api_base, timeout=DEFAULT_TIMEOUT):
    url = f"{api_base}/sequence/id/{feature_id}"
    params = {"format": "fasta", "type": "genomic"}
    headers = {"Content-Type": "text/x-fasta"}

    try:
        response = request(url, params, headers, timeout=timeout)
    except RequestErrors as e:
        raise ConnectionError(f_e("gff3", e))
    except Http400 as e:
        response_json = e.response.json()
        if response_json and response_json.get("error") == f"ID '{feature_id}' not found":
            raise NameError(f_e("fasta", e, response_json.get("error")))
        raise
    return response


def fetch_gff3(feature_id, api_base, timeout=DEFAULT_TIMEOUT):
    url = f"{api_base}/overlap/id/{feature_id}"
    params = {"feature": ["gene", "transcript", "cds", "exon"]}
    headers = {"Content-Type": "text/x-gff3"}

    try:
        response = request(url, params, headers, timeout=timeout)
    except RequestErrors as e:
        raise ConnectionError(f_e("gff3", e))
    except Http400 as e:
        response_json = e.response.json()
        if response_json and response_json.get("error") == f"ID '{feature_id}' not found":
            raise NameError(f_e("gff3", e, response_json.get("error")))
        raise
    return response


def _get_tark_versions(reference_id, api_base, timeout=DEFAULT_TIMEOUT):
    endpoint = "transcript"
    params = {"stable_id": reference_id}
    tark_req = json.loads(
        request(url=f"{api_base}/{endpoint}", params=params, timeout=timeout)
    )
    tark_versions_38 = []
    tark_versions_37 = []
    if tark_req["results"]:
        for r in tark_req["results"]:
            if r["assembly"] == "GRCh37":
                tark_versions_37.append(int(r["stable_id_version"]))
            elif r["assembly"] == "GRCh38":
                tark_versions_38.append(int(r["stable_id_version"]))

    return tark_versions_38, tark_versions_37


def _get_most_recent_version(reference_id, api_base, timeout=DEFAULT_TIMEOUT):
    return int(_get_reference_information(reference_id, api_base, timeout)["version"])


def _get_reference_information(reference_id, api_base, timeout=DEFAULT_TIMEOUT):
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


def fetch_json(
    reference_id,
    reference_version,
    api_base,
    assembly="GRCh38",
    timeout=DEFAULT_TIMEOUT,
):
    endpoint = "transcript"
    params = {
        "stable_id": reference_id,
        "assembly_name": assembly,
        "stable_id_version": reference_version,
        "expand": "translations, genes, exons",
    }
    req = requests.request(
        method="get", url=f"{api_base}/{endpoint}", params=params, timeout=timeout
    )
    return req.json()


def get_rest_api_base(r_id, r_version, timeout=DEFAULT_TIMEOUT):
    rest_version_38 = _get_most_recent_version(
        r_id, settings.get("ENSEMBL_API"), timeout
    )
    if r_version in [None, rest_version_38]:
        return settings.get("ENSEMBL_API"), "GRCh38"
    if r_version == _get_most_recent_version(
        r_id, settings.get("ENSEMBL_API_GRCH37"), timeout
    ):
        return settings.get("ENSEMBL_API_GRCH37"), "GRCh37"
    raise NameError(f"Cannot fetch {r_id}.{r_version} from Ensembl REST")


def get_transcript_api_base(r_id, r_version, r_source, timeout=DEFAULT_TIMEOUT):
    if r_source == "ensembl_rest":
        return get_rest_api_base(r_id, r_version, timeout)

    tark_versions_38, tark_versions_37 = _get_tark_versions(
        r_id, settings.get("ENSEMBL_TARK_API"), timeout
    )
    if r_version is None or r_version in tark_versions_38:
        return settings.get("ENSEMBL_TARK_API"), "GRCh38"
    if r_version in tark_versions_37:
        return settings.get("ENSEMBL_TARK_API"), "GRCh37"
    raise NameError(f"Cannot fetch {r_id} from Ensembl Tark")


def _fetch_default(r_id, r_version, reference_source, timeout):
    # An explicit "ensembl_tark" source goes straight to Tark, same as
    # reference_type="json" would. gff3 isn't served there anyway.
    if reference_source == "ensembl_tark":
        api_base, assembly = get_transcript_api_base(
            r_id, r_version, reference_source, timeout
        )
        return fetch_json(r_id, r_version, api_base, assembly, timeout), "json"

    # gff3 only exists on REST. Resolving or fetching it can fail as
    # NameError, ConnectionError, or a raw RequestErrors, all of which
    # fall back to Tark below.
    try:
        api_base, assembly = get_rest_api_base(r_id, r_version, timeout)
        return fetch_gff3(r_id, api_base, timeout), "gff3"
    except (NameError, ConnectionError, RequestErrors):
        if "ENST" not in r_id:
            raise
        api_base, assembly = get_transcript_api_base(
            r_id, r_version, reference_source, timeout
        )
        return fetch_json(r_id, r_version, api_base, assembly, timeout), "json"


def fetch(
    reference_id,
    reference_type=None,
    reference_source=None,
    timeout=DEFAULT_TIMEOUT,
):
    r_id, r_version = _get_id_and_version(reference_id)
    if r_id is None:
        raise NameError

    if reference_type is None:
        return _fetch_default(r_id, r_version, reference_source, timeout)

    # gff3/fasta only exist on REST, never on Tark.
    if reference_type in ["gff3", "fasta"]:
        api_base, assembly = get_rest_api_base(r_id, r_version, timeout)
    elif "ENST" in r_id:
        api_base, assembly = get_transcript_api_base(
            r_id, r_version, reference_source, timeout
        )
    else:
        api_base, assembly = get_rest_api_base(r_id, r_version, timeout)

    if reference_type == "gff3":
        return fetch_gff3(r_id, api_base, timeout), "gff3"
    elif reference_type == "fasta":
        return fetch_fasta(r_id, api_base, timeout), "fasta"
    elif reference_type == "json":
        if reference_source in [None, "ensembl", "ensembl_tark"]:
            return fetch_json(r_id, r_version, api_base, assembly, timeout), "json"

    elif reference_type == "genbank":
        return None, "genbank"

    raise ValueError(
        f"{reference_source} fetch does not support {reference_type} reference type."
    )
