import pytest

from mutalyzer_retriever.retriever import (
    NoReferenceError,
    NoReferenceRetrieved,
    _fetch_unknown_source,
    _raise_error,
    extract_feature_model,
)


def test_raise_error_all_name_errors_means_not_retrieved():
    status = {
        "lrg": {"errors": [NameError("not lrg")]},
        "ncbi": {"errors": [NameError("not found")]},
        "ensembl": {"errors": [NameError("not found")]},
    }
    with pytest.raises(NoReferenceRetrieved):
        _raise_error(status)


def test_raise_error_non_name_error_means_uncertain():
    status = {
        "lrg": {"errors": [NameError("not lrg")]},
        "ncbi": {"errors": [ConnectionError("network down")]},
        "ensembl": {"errors": [NameError("not found")]},
    }
    with pytest.raises(NoReferenceError) as exc_info:
        _raise_error(status)
    assert exc_info.value.uncertain_sources == ["ncbi"]


def test_fetch_unknown_source_lrg_succeeds_first(monkeypatch):
    def fetch_lrg(reference_id, timeout):
        return "lrg content"

    def fail(*args, **kwargs):
        raise AssertionError("should not be reached once lrg succeeds")

    monkeypatch.setattr("mutalyzer_retriever.retriever.lrg.fetch_lrg", fetch_lrg)
    monkeypatch.setattr("mutalyzer_retriever.retriever.ncbi.fetch", fail)
    monkeypatch.setattr("mutalyzer_retriever.retriever.ensembl.fetch", fail)

    assert _fetch_unknown_source("LRG_1", None, None) == ("lrg content", "lrg", "lrg")


def test_fetch_unknown_source_falls_back_lrg_to_ncbi(monkeypatch):
    def fail_lrg(reference_id, timeout):
        raise NameError("not an lrg id")

    def fetch_ncbi(reference_id, reference_type, size_off, timeout):
        return "ncbi content", "gff3"

    def fail_ensembl(*args, **kwargs):
        raise AssertionError("should not be reached once ncbi succeeds")

    monkeypatch.setattr("mutalyzer_retriever.retriever.lrg.fetch_lrg", fail_lrg)
    monkeypatch.setattr("mutalyzer_retriever.retriever.ncbi.fetch", fetch_ncbi)
    monkeypatch.setattr("mutalyzer_retriever.retriever.ensembl.fetch", fail_ensembl)

    assert _fetch_unknown_source("NM_1.1", None, None) == (
        "ncbi content",
        "gff3",
        "ncbi",
    )


def test_fetch_unknown_source_falls_back_ncbi_to_ensembl(monkeypatch):
    def fail(*args, **kwargs):
        raise NameError("not found")

    def fetch_ensembl(reference_id, reference_type, reference_source, timeout):
        return "ensembl content", "gff3"

    monkeypatch.setattr("mutalyzer_retriever.retriever.lrg.fetch_lrg", fail)
    monkeypatch.setattr("mutalyzer_retriever.retriever.ncbi.fetch", fail)
    monkeypatch.setattr("mutalyzer_retriever.retriever.ensembl.fetch", fetch_ensembl)

    assert _fetch_unknown_source("ENST00000000000", None, None) == (
        "ensembl content",
        "gff3",
        "ensembl",
    )


def test_fetch_unknown_source_all_fail_raises_not_retrieved(monkeypatch):
    def fail(*args, **kwargs):
        raise NameError("not found")

    monkeypatch.setattr("mutalyzer_retriever.retriever.lrg.fetch_lrg", fail)
    monkeypatch.setattr("mutalyzer_retriever.retriever.ncbi.fetch", fail)
    monkeypatch.setattr("mutalyzer_retriever.retriever.ensembl.fetch", fail)

    with pytest.raises(NoReferenceRetrieved):
        _fetch_unknown_source("BOGUS_ID", None, None)


def test_fetch_unknown_source_unexpected_ncbi_error_propagates(monkeypatch):
    def fail_lrg(reference_id, timeout):
        raise NameError("not an lrg id")

    def broken_ncbi(*args, **kwargs):
        raise RuntimeError("something genuinely broke")

    def fail_ensembl(*args, **kwargs):
        raise AssertionError("should not be reached, ncbi error must propagate")

    monkeypatch.setattr("mutalyzer_retriever.retriever.lrg.fetch_lrg", fail_lrg)
    monkeypatch.setattr("mutalyzer_retriever.retriever.ncbi.fetch", broken_ncbi)
    monkeypatch.setattr("mutalyzer_retriever.retriever.ensembl.fetch", fail_ensembl)

    with pytest.raises(RuntimeError, match="something genuinely broke"):
        _fetch_unknown_source("NM_1.1", None, None)


def _sample_annotations():
    return {
        "id": "NC_000001.1",
        "type": "record",
        "features": [
            {
                "id": "gene1",
                "type": "gene",
                "features": [
                    {
                        "id": "mrna1",
                        "type": "mRNA",
                        "features": [
                            {"id": "exon1", "type": "exon"},
                            {"id": "cds1", "type": "CDS"},
                        ],
                    },
                    {
                        "id": "mrna2",
                        "type": "mRNA",
                        "features": [{"id": "exon2", "type": "exon"}],
                    },
                ],
            }
        ],
    }


def test_extract_feature_model_not_found():
    assert extract_feature_model(_sample_annotations(), "nonexistent") == (
        None,
        False,
        False,
    )


def test_extract_feature_model_keeps_ancestors_and_descendants_by_default():
    model, _, _ = extract_feature_model(_sample_annotations(), "mrna1")
    # ancestor chain (record -> gene) is kept, mrna2 is pruned, mrna1 keeps its exons
    gene = model["features"][0]
    assert [f["id"] for f in gene["features"]] == ["mrna1"]
    assert [f["id"] for f in gene["features"][0]["features"]] == ["exon1", "cds1"]


def test_extract_feature_model_descendants_false_drops_matched_features():
    model, _, _ = extract_feature_model(
        _sample_annotations(), "mrna1", descendants=False
    )
    mrna = model["features"][0]["features"][0]
    assert mrna["id"] == "mrna1"
    assert "features" not in mrna


def test_extract_feature_model_ancestors_false_returns_match_only():
    model, _, _ = extract_feature_model(_sample_annotations(), "mrna1", ancestors=False)
    assert model["id"] == "mrna1"
    assert [f["id"] for f in model["features"]] == ["exon1", "cds1"]


def test_extract_feature_model_siblings_true_keeps_unmatched_siblings():
    model, _, _ = extract_feature_model(_sample_annotations(), "mrna1", siblings=True)
    gene = model["features"][0]
    assert [f["id"] for f in gene["features"]] == ["mrna1", "mrna2"]
