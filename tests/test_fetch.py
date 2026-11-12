import json
from urllib.error import HTTPError

import pytest
import requests

from mutalyzer_retriever.request import Http400, RequestErrors
from mutalyzer_retriever.sources import ensembl
from mutalyzer_retriever.sources.ensembl import (
    _get_reference_information,
    _get_tark_versions,
    fetch,
    fetch_gff3,
    fetch_json,
)
from mutalyzer_retriever.sources.lrg import fetch_lrg

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
        raise ConnectionError("503 Service Unavailable")

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


def test_lrg_fetch_non_ascii_id(monkeypatch):
    # Seen in the production logs: a full-width colon in the description.
    def urlopen(url, timeout):
        url.encode("ascii")  # what http.client does with the request line
        raise HTTPError(url, 404, "Not Found", None, None)

    monkeypatch.setattr("mutalyzer_retriever.sources.lrg.urlopen", urlopen)

    with pytest.raises(NameError):
        fetch_lrg("NM_145239.3：c.1013-29del")


def _http400(error):
    response = requests.Response()
    response.status_code = 400
    response._content = json.dumps({"error": error}).encode()
    return Http400(requests.exceptions.HTTPError(response=response))


@pytest.mark.parametrize(
    "function, r_id, error",
    [
        # Overlap does not support protein ids, the error is just the id.
        (fetch_gff3, "ENSP00000462054", "ENSP00000462054"),
        (
            _get_reference_information,
            "ENST00000000000",
            "ID 'ENST00000000000' not found",
        ),
    ],
)
def test_ensembl_http400_is_name_error(function, r_id, error, monkeypatch):
    def request(*args, **kwargs):
        raise _http400(error)

    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.request", request)

    with pytest.raises(NameError):
        function(r_id, "https://rest.ensembl.org", timeout=1)


@pytest.mark.parametrize("function", [_get_reference_information, _get_tark_versions])
def test_ensembl_request_errors_are_connection_errors(function, monkeypatch):
    def request(*args, **kwargs):
        raise RequestErrors([requests.exceptions.ReadTimeout("Read timed out.")])

    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.request", request)

    with pytest.raises(ConnectionError):
        function("ENST00000257430", "https://rest.ensembl.org", timeout=1)


@pytest.mark.parametrize(
    "status_code, content, error",
    [
        (200, {"results": []}, NameError),
        (503, {"detail": "Service temporarily unavailable"}, ConnectionError),
    ],
)
def test_ensembl_fetch_json_no_results(status_code, content, error, monkeypatch):
    response = requests.Response()
    response.status_code = status_code
    response._content = json.dumps(content).encode()

    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl.requests.request",
        lambda **kwargs: response,
    )

    with pytest.raises(error):
        fetch_json("ENST00000257430.4:c.835-8A>G", None, "https://tark.ensembl.org/api")


def test_ensembl_fetch_json_timeout(monkeypatch):
    def request(**kwargs):
        raise requests.exceptions.ReadTimeout("Read timed out. (read timeout=10)")

    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.requests.request", request)

    with pytest.raises(ConnectionError):
        fetch_json("ENST00000257430", 4, "https://tark.ensembl.org/api")
