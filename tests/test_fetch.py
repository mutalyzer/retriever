import pytest
from mutalyzer_retriever.configuration import settings
from mutalyzer_retriever.sources.ensembl import fetch
from .commons import _get_content, patch_retriever

API_BASE = settings["ENSEMBL_API"]
API_BASE_GRCH37 = settings["ENSEMBL_API_GRCH37"]
TARK_API_BASE = settings["ENSEMBL_TARK_API"]


API_BASE_MAP = {
    "ENSG00000147889": {"version": 18, "species": "homo_sapiens"},
    "ENSMUSG00000022346": {"version": 18, "species": "mus_musculus"},
    "ENST00000304494": {"version": 10, "species": "homo_sapiens"},
    "ENST00000000000": {"version": 20, "species": "homo_sapiens"},
}
API_BASE_GRCH37_MAP = {
    "ENSG00000147889": {"version": 12, "species": "homo_sapiens"},
    "ENST00000304494": {"version": 5, "species": "homo_sapiens"},
    "ENST00000000000": {"version": 6, "species": "homo_sapiens"},
}

TARK_API_BASE_MAP = {
    "ENST00000304494": {"GRCH38_version": [10, 9, 8, 7, 6], "GRCH37_version": [5]},
    "ENST00000000000": {"GRCH38_version": [20, 19, 18], "GRCH37_version": [6, 5]},
}


def _fetch_json(r_id, r_version, api_base, assembly, timeout):
    if api_base == TARK_API_BASE:
        return _get_content(
            f"data/{r_id}.{r_version}.tark_raw.model.json"
        )


def _get_tark_versions(r_id, api_base, timeout=1):
    if api_base == TARK_API_BASE and r_id in TARK_API_BASE_MAP.keys():
        return (
            TARK_API_BASE_MAP[r_id]["GRCH38_version"],
            TARK_API_BASE_MAP[r_id]["GRCH37_version"],
        )


def _fetch_gff3(feature_id, api_base, timeout=1):
    if api_base == API_BASE_GRCH37:
        return _get_content(
            f"data/{feature_id}.{API_BASE_GRCH37_MAP[feature_id]['version']}.gff3"
        )
    return _get_content(f"data/{feature_id}.gff3")


def _get_reference_information(r_id, api_base, timeout=1):
    if api_base == API_BASE and r_id in API_BASE_MAP.keys():
        return API_BASE_MAP[r_id]
    if api_base == API_BASE_GRCH37 and r_id in API_BASE_GRCH37_MAP.keys():
        return API_BASE_GRCH37_MAP[r_id]


@pytest.mark.parametrize("r_id", [("ENSG00000147889")])
def test_ensembl_fetch_no_version(r_id):
    assert fetch(r_id)[0] == _get_content(f"data/{r_id}.gff3")


@pytest.mark.parametrize("r_id", [("ENSG00000147889.18")])
def test_ensembl_fetch_version_newest(r_id):
    assert fetch(r_id)[0] == _get_content(f"data/{r_id}.gff3")


@pytest.mark.parametrize("r_id", [("ENST00000304494")])
def test_ensembl_fetch_transcript_no_version(r_id):
    assert fetch(r_id)[0] == _get_content(f"data/{r_id}.gff3")


@pytest.mark.parametrize("r_id", [("ENST00000304494")])
def test_ensembl_fetch_transcript_rest_38(r_id):
    assert fetch(r_id)[0] == _get_content(f"data/{r_id}.gff3")


@pytest.mark.parametrize("r_id, r_type, r_source", [("ENST00000304494.5", "json", "ensembl_rest")])
def test_ensembl_fetch_transcript_rest_37(r_id, r_type, r_source):
    with pytest.raises(ValueError):
        fetch(r_id, r_type, r_source)


@pytest.mark.parametrize("r_id, r_type", [("ENST00000304494.7", "json")])
def test_ensembl_fetch_transcript_tark_38(r_id, r_type):
    assert fetch(r_id, r_type)[0] == _get_content(
        f"data/{r_id}.tark_raw.model.json"
    )


@pytest.mark.parametrize("r_id, r_type", [("ENST00000000000.5", "json")])
def test_ensembl_fetch_transcript_tark_37(r_id, r_type):
    assert fetch(r_id, r_type)[0] == _get_content(
        f"data/{r_id}.tark_raw.model.json"
    )


@pytest.mark.parametrize("r_id", [("ENSG00000147889.12")])
def test_ensembl_fetch_version_grch37(r_id):
    assert fetch(r_id)[0] == _get_content(f"data/{r_id}.gff3")


@pytest.mark.parametrize("r_id", [("ENSG00000147889.15")])
def test_ensembl_fetch_other_version(r_id):
    with pytest.raises(NameError):
        assert fetch(r_id)[0] == None


@pytest.mark.parametrize("r_id", [("ENSMUSG00000022346.18")])
def test_ensembl_fetch_no_version_mouse(r_id):
    assert fetch(r_id)[0] == _get_content(f"data/{r_id}.gff3")


@pytest.mark.parametrize("r_id", [("ENSMUSG00000022346")])
def test_ensembl_fetch_version_newest_mouse(r_id):
    assert fetch(r_id)[0] == _get_content(f"data/{r_id}.gff3")
