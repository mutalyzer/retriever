import json
from pathlib import Path

import pytest

from mutalyzer_retriever.related import get_related


def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    with open(str(data_file), "r") as file:
        content = file.read()
    return content


def _fetch_1(
    reference_id,
    timeout=1,
):
    return _get_content("data/" + reference_id + ".ncbi_datasets_gene_accession.json")


def _fetch_2(
    gene_id,
    timeout=1,
):
    return _get_content("data/gene_summary_" + gene_id + ".json")


def get_tests(references):

    tests = []

    for r_id in references:
        p = Path(Path(__file__).parent) / "data" / str(r_id + ".related.json")
        with p.open() as f:
            r_model = json.loads(f.read())
        tests.append(
            pytest.param(
                r_id,
                r_model,
                id="{}".format(r_id),
            )
        )

    return tests


@pytest.mark.parametrize(
    "r_id, expected_model",
    get_tests(["NM_021803.4", "NM_003002.2"]),
)
def test_model(r_id, expected_model, monkeypatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.related._fetch_ncbi_datasets_gene_accession", _fetch_1
    )
    monkeypatch.setattr(
        "mutalyzer_retriever.related._fetch_ncbi_entrez_eutils_esummary", _fetch_2
    )
    related = get_related(r_id)
    assert related.keys() == expected_model.keys()
    for k in related:
        assert set(related[k]) == set(expected_model[k])
