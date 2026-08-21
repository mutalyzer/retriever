import pytest

from mutalyzer_retriever.request import RequestErrors
from mutalyzer_retriever.sources import ensembl
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
        fetch(r_id)


@pytest.mark.parametrize("r_id", [("ENSMUSG00000022346.18")])
def test_ensembl_fetch_no_version_mouse(r_id):
    assert fetch(r_id)[0] == _get_content(f"data/{r_id}.gff3")


@pytest.mark.parametrize("r_id", [("ENSMUSG00000022346")])
def test_ensembl_fetch_version_newest_mouse(r_id):
    assert fetch(r_id)[0] == _get_content(f"data/{r_id}.gff3")


def test_ensembl_fetch_default_tries_gff3_via_rest_not_tark(monkeypatch):
    def get_reference_information(reference_id, api_base, timeout):
        return {"version": "10"}

    def get_tark_versions(reference_id, api_base, timeout):
        raise AssertionError("should not resolve via Tark when gff3 succeeds")

    calls = []

    def fetch_gff3(feature_id, api_base, timeout):
        calls.append(api_base)
        return "##gff-version 3\n"

    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl._get_reference_information",
        get_reference_information,
    )
    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl._get_tark_versions", get_tark_versions
    )
    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.fetch_gff3", fetch_gff3)

    _, reference_type = fetch("ENST00000304494", timeout=1)

    assert reference_type == "gff3"
    assert calls == [ensembl.settings.get("ENSEMBL_API")]


def test_ensembl_fetch_explicit_gff3_type_uses_rest_not_tark(monkeypatch):
    def get_reference_information(reference_id, api_base, timeout):
        return {"version": "10"}

    def get_tark_versions(reference_id, api_base, timeout):
        raise AssertionError("gff3 must not resolve via Tark")

    calls = []

    def fetch_gff3(feature_id, api_base, timeout):
        calls.append(api_base)
        return "##gff-version 3\n"

    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl._get_reference_information",
        get_reference_information,
    )
    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl._get_tark_versions", get_tark_versions
    )
    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.fetch_gff3", fetch_gff3)

    _, reference_type = fetch("ENST00000304494", reference_type="gff3", timeout=1)

    assert reference_type == "gff3"
    assert calls == [ensembl.settings.get("ENSEMBL_API")]


def test_ensembl_fetch_default_falls_back_to_tark_when_rest_lacks_version(monkeypatch):
    # e.g. an old ENST version REST no longer serves gff3/fasta for, but
    # that still exists in Tark's history.
    def get_rest_api_base(r_id, r_version, timeout):
        raise NameError(f"Cannot fetch {r_id}.{r_version} from Ensembl REST")

    calls = []

    def get_transcript_api_base(r_id, r_version, r_source, timeout):
        calls.append((r_id, r_version))
        return ensembl.settings.get("ENSEMBL_TARK_API"), "GRCh38"

    def fetch_json(r_id, r_version, api_base, assembly, timeout):
        return {"stable_id": r_id, "stable_id_version": r_version}

    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl.get_rest_api_base", get_rest_api_base
    )
    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl.get_transcript_api_base",
        get_transcript_api_base,
    )
    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.fetch_json", fetch_json)

    _, reference_type = fetch("ENST00000304494.7", timeout=1)

    assert reference_type == "json"
    assert calls == [("ENST00000304494", 7)]


def test_ensembl_fetch_default_falls_back_to_tark_on_rest_request_error(monkeypatch):
    def get_rest_api_base(r_id, r_version, timeout):
        raise RequestErrors(["503 Service Unavailable"])

    def get_transcript_api_base(r_id, r_version, r_source, timeout):
        return ensembl.settings.get("ENSEMBL_TARK_API"), "GRCh38"

    def fetch_json(r_id, r_version, api_base, assembly, timeout):
        return {"stable_id": r_id}

    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl.get_rest_api_base", get_rest_api_base
    )
    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl.get_transcript_api_base",
        get_transcript_api_base,
    )
    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.fetch_json", fetch_json)

    _, reference_type = fetch("ENST00000304494.7", timeout=1)

    assert reference_type == "json"


def test_ensembl_fetch_default_non_enst_reraises_when_rest_fails(monkeypatch):
    def get_rest_api_base(r_id, r_version, timeout):
        raise NameError("not found")

    def fail_transcript(*args, **kwargs):
        raise AssertionError("should not attempt Tark for a non-ENST id")

    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl.get_rest_api_base", get_rest_api_base
    )
    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl.get_transcript_api_base", fail_transcript
    )

    with pytest.raises(NameError):
        fetch("ENSG00000147889.99", timeout=1)


def test_ensembl_fetch_default_explicit_tark_source_uses_tark_not_rest(monkeypatch):
    def fail_rest(*args, **kwargs):
        raise AssertionError("should not try REST when ensembl_tark is explicit")

    calls = []

    def get_transcript_api_base(r_id, r_version, r_source, timeout):
        calls.append((r_id, r_version, r_source))
        return ensembl.settings.get("ENSEMBL_TARK_API"), "GRCh38"

    def fetch_json(r_id, r_version, api_base, assembly, timeout):
        return {"stable_id": r_id}

    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl.get_rest_api_base", fail_rest
    )
    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.fetch_gff3", fail_rest)
    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl.get_transcript_api_base",
        get_transcript_api_base,
    )
    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.fetch_json", fetch_json)

    _, reference_type = fetch(
        "ENST00000304494", reference_source="ensembl_tark", timeout=1
    )

    assert reference_type == "json"
    assert calls == [("ENST00000304494", None, "ensembl_tark")]
