import json
from pathlib import Path

import pytest

from mutalyzer_retriever import retrieve_model

from .commons import _get_content, patch_retriever, references


def get_tests(references):

    tests = []

    for r_source in references.keys():
        for r_type in references[r_source].keys():
            for r_id in references[r_source][r_type]:
                if r_type == "json":
                    p = (
                        Path(Path(__file__).parent)
                        / "data"
                        / str(r_id + ".tark.model.json")
                    )
                else:
                    p = Path(Path(__file__).parent) / "data" / str(r_id + ".model.json")
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


@pytest.mark.parametrize(
    "r_id, r_source, r_type, expected_model", get_tests(references)
)
def test_model(r_id, r_source, r_type, expected_model, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.parsers.json_ensembl._seq_from_rest",
        lambda _0, _1, _2, _3, _4: _seq_from_rest(r_id),
    )
    assert retrieve_model(r_id, r_source, r_type) == expected_model
