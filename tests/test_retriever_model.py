import json
from pathlib import Path

import pytest
from .commons import references, _retrieve_raw

from mutalyzer_retriever import retrieve_model


def get_tests(references):

    tests = []

    for r_source in references.keys():
        for r_type in references[r_source].keys():
            for r_id in references[r_source][r_type]:
                if r_type == "json":
                    p = Path(Path(__file__).parent) / "data" / str(r_id + ".tark.model.json")
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
                        id="{}-{}-{}".format(r_source, r_type, r_id),
                    )
                )

    return tests


@pytest.mark.parametrize("r_id, r_source, r_type, expected_model",get_tests(references))
def test_model(r_id, r_source, r_type, expected_model, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr("mutalyzer_retriever.retriever.retrieve_raw", _retrieve_raw)
    assert retrieve_model(r_id, r_source,r_type) == expected_model

