from pathlib import Path

import pytest

from mutalyzer_retriever.sources.ensembl import fetch
from mutalyzer_retriever.configuration import settings

API_BASE = settings["ENSEMBL_API"]
API_BASE_GRCH37 = settings["ENSEMBL_API_GRCH37"]

API_BASE_MAP = {
    "ENSG00000147889": {"version": 18, "species": "homo_sapiens"},
    "ENSMUSG00000022346": {"version": 18, "species": "mus_musculus"},
}
API_BASE_GRCH37_MAP = {"ENSG00000147889": {"version": 12, "species": "homo_sapiens"}}


@pytest.fixture(autouse=True)
def patch_retriever(monkeypatch):
    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.fetch_gff3", _fetch_gff3)
    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl._get_reference_information",
        _get_reference_information,
    )


def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    with open(str(data_file), "r") as file:
        content = file.read()
    return content


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


@pytest.mark.parametrize("reference_id", [("ENSG00000147889.12")])
def test_ensembl_fetch_version_grch37(reference_id):
    assert fetch(reference_id)[0] == _get_content(f"data/{reference_id}.gff3")


@pytest.mark.parametrize("reference_id", [("ENSG00000147889.15")])
def test_ensembl_fetch_other_version(reference_id):
    with pytest.raises(NameError):
        fetch(reference_id)[0]


@pytest.mark.parametrize("reference_id", [("ENSMUSG00000022346.18")])
def test_ensembl_fetch_no_version_mouse(reference_id):
    assert fetch(reference_id)[0] == _get_content(f"data/{reference_id}.gff3")


@pytest.mark.parametrize("reference_id", [("ENSMUSG00000022346")])
def test_ensembl_fetch_version_newest_mouse(reference_id):
    assert fetch(reference_id)[0] == _get_content(f"data/{reference_id}.gff3")
