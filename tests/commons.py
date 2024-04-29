import json
from pathlib import Path
import pytest
from mutalyzer_retriever.retriever import NoReferenceError


@pytest.fixture(autouse=True)
def patch_retriever(monkeypatch):
    """retrieve all monkeypath"""
    from .test_fetch import (
        _fetch_gff3,
        _fetch_json,
        _get_reference_information,
        _get_tark_versions,
    )

    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.fetch_gff3", _fetch_gff3)
    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl._get_reference_information",
        _get_reference_information,
    )
    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl._get_tark_versions", _get_tark_versions
    )
    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.fetch_json", _fetch_json)
    monkeypatch.setattr("mutalyzer_retriever.retriever.retrieve_raw", _retrieve_raw)


def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    try:
        with open(str(data_file), "r") as file:
            content = file.read()
    except FileNotFoundError as exc:
        raise NoReferenceError({}, []) from exc
    return content


def _retrieve_raw(
    r_id,
    r_source=None,
    r_type=None,
    size_off=True,
    configuration_path=None,
    timeout=1,
):
    if r_type == "fasta":
        return _get_content("data/" + r_id + ".fasta"), "fasta", "ncbi"
    elif r_id.startswith("LRG_"):
        return _get_content("data/" + r_id + ".lrg"), "lrg", "lrg"
    elif r_type == "json":
        return (
            json.loads(_get_content("data/" + r_id + ".tark_raw.json")),
            "json",
            "ensembl_tark",
        )
    else:
        return _get_content("data/" + r_id + ".gff3"), "gff3", "ncbi"


references = {
    "ncbi": {
        "gff3": [
            "NM_078467.2",
            "NM_152263.2",
            "NM_152263.3",
            "NM_000077.4",
            "NM_002001.2",
            "NG_012337.1",
            "NR_002196.2",
            "L41870.1",
            "NG_007485.1",
            "NC_012920.1",
            "NG_009930.1",
            "AA010203.1",
            "NP_060665.3",
            "D64137.1",
            "AB006684.1",
            "NM_004152.3",
            "7",
            "M65131.1",
            "XR_948219.2",
            "NR_023343.1",
        ]
    },
    "ensembl_rest": {
        "gff3": [
            "ENSG00000147889",
            "ENST00000383925",
            "ENST00000304494",
            "ENSG00000198899",
        ]
    },
    "ensembl_tark": {
        "json": [
            "ENST00000383925.1",
            "ENST00000383925",
            "ENST00000304494",
            "ENST00000304494.10",
        ]
    },
    "lrg": {"lrg": ["LRG_11", "LRG_417", "LRG_857"]},
}
