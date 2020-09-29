import json
from pathlib import Path

import pytest
from mutalyzer_retriever import retrieve_model


def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    with open(str(data_file), "r") as file:
        content = file.read()
    return content


def ordered(obj):
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return sorted(ordered(x) for x in obj)
    else:
        return obj


def get_tests(references):

    tests = []

    for r_source in references.keys():
        for r_type in references[r_source].keys():
            for r_id in references[r_source][r_type]:
                p = Path(Path(__file__).parent) / "data" / str(r_id + ".model.json")
                with p.open() as f:
                    r_model = json.loads(f.read())
                tests.append((r_id, r_source, r_type, r_model))

    return tests


@pytest.mark.parametrize(
    "r_id, r_source, r_type, expected_model",
    get_tests(
        {
            "ncbi": {
                "gff3": [
                    "NM_078467.2",
                    "NM_152263.2",
                    "NM_152263.3",
                    "NM_000077.4",
                    "NG_012337.1",
                    "NR_002196.2",
                    "L41870.1",
                    "NG_007485.1",
                ]
            },
            "ensembl": {"gff3": ["ENSG00000147889"]},
            "lrg": {"lrg": ["LRG_11"]},
        }
    ),
)
def test_model(r_id, r_source, r_type, expected_model, monkeypatch):
    def mock_fetch_gff3(reference_id):
        return _get_content("data/" + reference_id + ".gff3")

    def mock_fetch_fasta(reference_id):
        return _get_content("data/" + reference_id + ".fasta")

    def mock_fetch_lrg(reference_id):
        return _get_content("data/" + reference_id + ".lrg")

    def mock_fetch_name_error(reference_id):
        raise NameError

    if r_source == 'ncbi':
        monkeypatch.setattr(
            "mutalyzer_retriever.sources.ncbi.fetch_gff3", mock_fetch_gff3)
        monkeypatch.setattr(
            "mutalyzer_retriever.sources.ncbi.fetch_fasta", mock_fetch_fasta )
        monkeypatch.setattr(
            "mutalyzer_retriever.sources.ensembl.fetch_gff3", mock_fetch_name_error)
        monkeypatch.setattr(
            "mutalyzer_retriever.sources.lrg.fetch_lrg", mock_fetch_name_error)
    elif r_source == 'ensembl':
        monkeypatch.setattr(
            "mutalyzer_retriever.sources.ensembl.fetch_gff3", mock_fetch_gff3)
        monkeypatch.setattr(
            "mutalyzer_retriever.sources.ensembl.fetch_fasta", mock_fetch_fasta)
        monkeypatch.setattr(
            "mutalyzer_retriever.sources.ncbi.fetch_gff3", mock_fetch_name_error)
        monkeypatch.setattr(
            "mutalyzer_retriever.sources.lrg.fetch_lrg", mock_fetch_name_error)
    elif r_source == 'lrg':
        monkeypatch.setattr(
            "mutalyzer_retriever.sources.lrg.fetch_lrg", mock_fetch_lrg)
        monkeypatch.setattr(
            "mutalyzer_retriever.sources.ensembl.fetch_gff3", mock_fetch_name_error)
        monkeypatch.setattr(
            "mutalyzer_retriever.sources.ncbi.fetch_gff3", mock_fetch_name_error )

    assert ordered(retrieve_model(r_id, r_source)) == ordered(expected_model)