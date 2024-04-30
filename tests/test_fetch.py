import pytest

from mutalyzer_retriever.sources.ensembl import fetch

from .commons import _get_content, patch_retriever


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


@pytest.mark.parametrize(
    "r_id, r_type, r_source", [("ENST00000304494.5", "json", "ensembl_rest")]
)
def test_ensembl_fetch_transcript_rest_37(r_id, r_type, r_source):
    with pytest.raises(ValueError):
        fetch(r_id, r_type, r_source)


@pytest.mark.parametrize("r_id, r_type", [("ENST00000304494.7", "json")])
def test_ensembl_fetch_transcript_tark_38(r_id, r_type):
    assert fetch(r_id, r_type)[0] == _get_content(f"data/{r_id}.tark_raw.model.json")


@pytest.mark.parametrize("r_id, r_type", [("ENST00000000000.5", "json")])
def test_ensembl_fetch_transcript_tark_37(r_id, r_type):
    assert fetch(r_id, r_type)[0] == _get_content(f"data/{r_id}.tark_raw.model.json")


@pytest.mark.parametrize("r_id", [("ENSG00000147889.12")])
def test_ensembl_fetch_version_grch37(r_id):
    assert fetch(r_id)[0] == _get_content(f"data/{r_id}.gff3")


@pytest.mark.parametrize("r_id", [("ENSG00000147889.15")])
def test_ensembl_fetch_other_version(r_id):
    with pytest.raises(NameError):
        fetch(r_id)[0] is None


@pytest.mark.parametrize("r_id", [("ENSMUSG00000022346.18")])
def test_ensembl_fetch_no_version_mouse(r_id):
    assert fetch(r_id)[0] == _get_content(f"data/{r_id}.gff3")


@pytest.mark.parametrize("r_id", [("ENSMUSG00000022346")])
def test_ensembl_fetch_version_newest_mouse(r_id):
    assert fetch(r_id)[0] == _get_content(f"data/{r_id}.gff3")
