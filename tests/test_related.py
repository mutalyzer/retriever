import hashlib
import json
from pathlib import Path

import pytest

from mutalyzer_retriever.related import get_related


def _get_hash(uids):
    uids.sort()
    s_i = ",".join(uids)
    return hashlib.md5(s_i.encode("utf-8")).hexdigest()


def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    with open(str(data_file), "r") as file:
        content = file.read()
    return content


def _fetch_ncbi_esummary(db, query_id, timeout=1):
    if "," in query_id:
        query_id = _get_hash(query_id.split(","))
    return json.loads(_get_content(f"data/esummary_{db}_{query_id}.json"))


def _fetch_ncbi_elink(db, dbfrom, query_id, timeout=1):
    return json.loads(_get_content(f"data/elink_{db}_{dbfrom}_{query_id}.json"))


def _fetch_ncbi_datasets_gene_accession(accession_id, timeout=1):
    return json.loads(
        _get_content(f"data/{accession_id}.ncbi_datasets_gene_accession.json")
    )


def _fetch_ensembl_xrefs(query_id, timeout=1):
    if query_id.startswith("LRG") or query_id.startswith("ENS"):
        return json.loads(_get_content(f"data/ensembl_xrefs_{query_id}.json"))
    else:
        return []


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
    get_tests(
        [
            "NM_021803.4",
            "NM_003002.2",
            "NM_003002.4",
            "NR_002196.2",
            "NG_012337.1",
            "NG_012337.2",
            "NG_012337.3",
            "LRG_24",
            "NC_000022.10",
            "NC_000022.11",
            "NC_000022",
            "CYP2D6",
            "ENSG00000159189",
        ]
    ),
)
def test_model(r_id, expected_model, monkeypatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.related._fetch_ncbi_esummary", _fetch_ncbi_esummary
    )
    monkeypatch.setattr(
        "mutalyzer_retriever.related._fetch_ncbi_elink", _fetch_ncbi_elink
    )
    monkeypatch.setattr(
        "mutalyzer_retriever.related._fetch_ncbi_datasets_gene_accession",
        _fetch_ncbi_datasets_gene_accession,
    )
    monkeypatch.setattr(
        "mutalyzer_retriever.related._fetch_ensembl_xrefs",
        _fetch_ensembl_xrefs,
    )

    assert get_related(r_id) == expected_model
