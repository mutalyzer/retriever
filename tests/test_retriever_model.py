import json
from pathlib import Path

import pytest

from mutalyzer_retriever import retrieve_model
from mutalyzer_retriever.retriever import retrieve_model_from_file

from .commons import _get_content, patch_retriever, references


def get_tests(references):

    tests = []

    for r_source in references:
        for r_type in references[r_source]:
            for r_id in references[r_source][r_type]:
                if r_type == "json":
                    p = Path(__file__).parent / "data" / str(r_id + ".tark.model.json")
                else:
                    p = Path(__file__).parent / "data" / str(r_id + ".model.json")
                with p.open() as f:
                    r_model = json.loads(f.read())
                tests.append(
                    pytest.param(
                        r_id,
                        r_source,
                        r_type,
                        r_model,
                        id=f"{r_source}-{r_type}-{r_id}",
                    )
                )

    return tests


def _seq_from_rest(r_id):
    return _get_content("data/" + str(r_id) + ".sequence")


def load_expected_model(r_id):
    with (Path(__file__).parent / "data" / f"{r_id}.model.json").open() as f:
        return json.load(f)


@pytest.mark.parametrize(
    "r_id, r_source, r_type, expected_model", get_tests(references)
)
def test_model(r_id, r_source, r_type, expected_model, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.parsers.json_ensembl._seq_from_rest",
        lambda _0, _1, _2, _3, _4, _5: _seq_from_rest(r_id),
    )
    assert retrieve_model(r_id, r_source, r_type) == expected_model


def test_model_type_sequence():
    # model_type only changes what retrieve_model() returns, not per-reference
    # behaviour, so a single fixture is enough to cover the "gff3" branch.
    r_id = "NM_000077.4"
    expected_model = load_expected_model(r_id)

    assert retrieve_model(r_id, "ncbi", "gff3", model_type="sequence") == {
        "sequence": expected_model["sequence"]
    }


def test_model_type_annotations():
    r_id = "NM_000077.4"
    expected_model = load_expected_model(r_id)

    assert (
        retrieve_model(r_id, "ncbi", "gff3", model_type="annotations")
        == expected_model["annotations"]
    )


def test_retrieve_model_from_file_gff3_and_fasta():
    data_dir = Path(__file__).parent / "data"
    model = retrieve_model_from_file(
        paths=[str(data_dir / "AA010203.1.gff3"), str(data_dir / "AA010203.1.fasta")]
    )
    with (data_dir / "AA010203.1.model.json").open() as f:
        expected_model = json.load(f)
    assert model == expected_model


def test_retrieve_model_from_file_lrg():
    data_dir = Path(__file__).parent / "data"
    model = retrieve_model_from_file(paths=[str(data_dir / "LRG_11.lrg")], is_lrg=True)
    with (data_dir / "LRG_11.model.json").open() as f:
        expected_model = json.load(f)
    assert model == expected_model
