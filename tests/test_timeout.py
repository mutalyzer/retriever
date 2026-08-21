from mutalyzer_retriever.retriever import retrieve_raw
from mutalyzer_retriever.sources import ensembl


def test_retrieve_raw_ncbi_forwards_size_and_timeout(monkeypatch):
    calls = {}

    def fetch(reference_id, reference_type, size_on, timeout):
        calls["reference_id"] = reference_id
        calls["reference_type"] = reference_type
        calls["size_on"] = size_on
        calls["timeout"] = timeout
        return "raw", "gff3"

    retrieve_raw.cache_clear()
    monkeypatch.setattr("mutalyzer_retriever.retriever.ncbi.fetch", fetch)

    assert retrieve_raw("NM_000077.4", "ncbi", "gff3", False, timeout=42) == (
        "raw",
        "gff3",
        "ncbi",
    )
    assert calls == {
        "reference_id": "NM_000077.4",
        "reference_type": "gff3",
        "size_on": False,
        "timeout": 42,
    }

    retrieve_raw.cache_clear()


def test_get_rest_api_base_forwards_timeout(monkeypatch):
    calls = []

    def get_reference_information(reference_id, api_base, timeout):
        calls.append((reference_id, api_base, timeout))
        return {"version": "18"}

    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl._get_reference_information",
        get_reference_information,
    )

    ensembl.get_rest_api_base("ENSG00000147889", None, timeout=42)

    assert calls == [("ENSG00000147889", ensembl.settings.get("ENSEMBL_API"), 42)]


def test_get_transcript_api_base_forwards_timeout(monkeypatch):
    calls = []

    def get_tark_versions(reference_id, api_base, timeout):
        calls.append((reference_id, api_base, timeout))
        return [10], [5]

    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl._get_tark_versions", get_tark_versions
    )

    ensembl.get_transcript_api_base("ENST00000304494", 10, "ensembl_tark", timeout=42)

    assert calls == [("ENST00000304494", ensembl.settings.get("ENSEMBL_TARK_API"), 42)]
