import json
from pathlib import Path

import pytest

from mutalyzer_retriever import retrieve_model


def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    with open(str(data_file), "r") as file:
        content = file.read()
    return content


def _retrieve_raw(
    reference_id,
    reference_source=None,
    reference_type=None,
    size_off=True,
    configuration_path=None,
    timeout=1,
):
    if reference_type == "fasta":
        return _get_content("data/" + reference_id + ".fasta"), "fasta", "ncbi"
    elif reference_id.startswith("LRG_"):
        return _get_content("data/" + reference_id + ".lrg"), "lrg", "lrg"
    else:
        return _get_content("data/" + reference_id + ".gff3"), "gff3", "ncbi"


def get_tests(references):

    tests = []

    for r_source in references.keys():
        for r_type in references[r_source].keys():
            for r_id in references[r_source][r_type]:
                p = Path(Path(__file__).parent) / "data" / str(r_id + ".model.json")
                with p.open() as f:
                    r_model = json.loads(f.read())
                tests.append(
                    pytest.param(
                        r_id,
                        r_source,
                        r_type,
                        r_model,
                        id="{}-{}-{}".format(r_source, r_type, r_id),
                    )
                )

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
                    "NM_002001.2",
                    "NG_012337.1",
                    "NR_002196.2",
                    "L41870.1",
                    "NG_007485.1",
                    "NC_012920.1",
                    "NG_009930.1",
                    "AA010203.1",
                    "NP_060665.3",
                    "D64137.1",
                    "AB006684.1",
                    "NM_004152.3",
                    "7",
                    "M65131.1",
                    "XR_948219.2",
                    "NR_023343.1",
                ]
            },
            "ensembl": {
                "gff3": [
                    "ENSG00000147889",
                    "ENST00000383925",
                    "ENST00000304494",
                    "ENSG00000198899",
                ]
            },
            "lrg": {"lrg": ["LRG_11", "LRG_417", "LRG_857"]},
        }
    ),
)
def test_model(r_id, r_source, r_type, expected_model, monkeypatch):
    monkeypatch.setattr("mutalyzer_retriever.retriever.retrieve_raw", _retrieve_raw)

    assert retrieve_model(r_id, r_source) == expected_model
