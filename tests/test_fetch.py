from pathlib import Path

import pytest
import json

from mutalyzer_retriever.sources.ensembl import fetch
from mutalyzer_retriever.configuration import settings

from .commons import _get_content

API_BASE = settings["ENSEMBL_API"]
API_BASE_GRCH37 = settings["ENSEMBL_API_GRCH37"]
TARK_API_BASE = settings["ENSEMBL_TARK_API"]


API_BASE_MAP = {
    "ENSG00000147889": {"version": 18, "species": "homo_sapiens"},
    "ENSMUSG00000022346": {"version": 18, "species": "mus_musculus"},
    "ENST00000304494": {"version": 10, "species": "homo_sapiens"},
}
API_BASE_GRCH37_MAP = {
    "ENSG00000147889": {"version": 12, "species": "homo_sapiens"},
    "ENST00000304494": {"version": 5, "species": "homo_sapiens"},
}

TARK_API_BASE_MAP = {
    "ENST00000304494": {"GRCH38_version": [10,9,8,7,6], "GRCH37_version":[5]}   
}



@pytest.fixture(autouse=True)
def patch_retriever(monkeypatch):
    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.fetch_gff3", _fetch_gff3)
    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl._get_reference_information",
        _get_reference_information,
    )
    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl._get_tark_versions",_get_tark_versions)
    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.fetch_tark",_fetch_tark)
    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.fetch_json",_fetch_json)


def _fetch_tark(reference_id, reference_version, api_base, assembly):
    if api_base == TARK_API_BASE:
        print(f"data/{reference_id}.{reference_version}.tark_raw.model.json")
        return _get_content(f"data/{reference_id}.{reference_version}.tark_raw.model.json")

def _get_tark_versions(reference_id, api_base, timeout=1):
    if api_base == TARK_API_BASE and reference_id in TARK_API_BASE_MAP.keys():
        return TARK_API_BASE_MAP[reference_id]["GRCH38_version"],TARK_API_BASE_MAP[reference_id]["GRCH37_version"]


def _fetch_json(reference_id, api_base, timeout=1):
    if api_base == API_BASE:
        return _get_content(f"data/{reference_id}.rest_raw.json")

def _fetch_gff3(feature_id, api_base, timeout=1):
    if api_base == API_BASE_GRCH37:
        return _get_content(
            f"data/{feature_id}.{API_BASE_GRCH37_MAP[feature_id]['version']}.gff3"
        )
    return _get_content(f"data/{feature_id}.gff3")


def _get_reference_information(reference_id, api_base, timeout=1):
    if api_base == API_BASE and reference_id in API_BASE_MAP.keys():
        return API_BASE_MAP[reference_id]
    if api_base == API_BASE_GRCH37 and reference_id in API_BASE_GRCH37_MAP.keys():
        return API_BASE_GRCH37_MAP[reference_id]


@pytest.mark.parametrize("reference_id", [("ENSG00000147889")])
def test_ensembl_fetch_no_version(reference_id):
    assert fetch(reference_id)[0] == _get_content(f"data/{reference_id}.gff3")


@pytest.mark.parametrize("reference_id", [("ENSG00000147889.18")])
def test_ensembl_fetch_version_newest(reference_id):
    assert fetch(reference_id)[0] == _get_content(f"data/{reference_id}.gff3")


@pytest.mark.parametrize("reference_id", [("ENST00000304494")])
def test_ensembl_fetch_transcript_no_version(reference_id):
    assert fetch(reference_id)[0] == _get_content(f"data/{reference_id}.gff3")


@pytest.mark.parametrize("reference_id, reference_type", [("ENST00000304494", 'json')])
def test_ensembl_fetch_transcript_rest_38(reference_id, reference_type):
    assert fetch(reference_id,reference_type)[0] == _get_content(f"data/{reference_id}.rest_raw.json")


@pytest.mark.parametrize("reference_id, reference_type", [("ENST00000304494.5", 'json')])
def test_ensembl_fetch_transcript_rest_37(reference_id, reference_type):
    ## for transcripts with version, go to tark api
    pass


@pytest.mark.parametrize("reference_id, reference_type", [("ENST00000304494.7", 'json')])
def test_ensembl_fetch_transcript_tark_38(reference_id, reference_type):
    assert fetch(reference_id, reference_type)[0] == _get_content(f"data/{reference_id}.tark_raw.model.json")


def test_ensembl_fetch_transcript_tark_37():
    ## have not found a transcript 37 that is unique in tark
    pass


@pytest.mark.parametrize("reference_id", [("ENSG00000147889.12")])
def test_ensembl_fetch_version_grch37(reference_id):
    assert fetch(reference_id)[0] == _get_content(f"data/{reference_id}.gff3")


@pytest.mark.parametrize("reference_id", [("ENSG00000147889.15")])
def test_ensembl_fetch_other_version(reference_id):
    with pytest.raises(ValueError):
        fetch(reference_id)[0]


@pytest.mark.parametrize("reference_id", [("ENSMUSG00000022346.18")])
def test_ensembl_fetch_no_version_mouse(reference_id):
    assert fetch(reference_id)[0] == _get_content(f"data/{reference_id}.gff3")


@pytest.mark.parametrize("reference_id", [("ENSMUSG00000022346")])
def test_ensembl_fetch_version_newest_mouse(reference_id):
    assert fetch(reference_id)[0] == _get_content(f"data/{reference_id}.gff3")
