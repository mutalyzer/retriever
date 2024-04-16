import json
from pathlib import Path

import pytest

from mutalyzer_retriever import parser
from mutalyzer_retriever.schema_validation import validate
from .test_retriever_model import _seq_from_rest

from .commons import references, patch_retriever



def get_references_content(references):
    '''read raw response from tests data folder'''
    r_contents = []
    for r_source in references.keys():
        for r_type in references[r_source]:
            for r_id in references[r_source][r_type]:
                if r_type == "json":
                    path_gb = (
                        Path(Path(__file__).parent)
                        / "data"
                        / f"{r_id}.tark.model.{r_type}"
                    )
                    r_content = json.loads(path_gb.open().read())
                else:
                    path_gb = (
                        Path(Path(__file__).parent)
                        / "data"
                        / f"{r_id}.{r_type}"
                    )
                    with path_gb.open() as f:
                        r_content = f.read()
                r_contents.append(
                    pytest.param(
                        r_source,
                        r_type,
                        r_content,
                        r_id,
                        id=f"{r_source}-{r_type}-{r_id}",
                    )
                )
    return r_contents


@pytest.mark.parametrize(
    "r_source, r_type, r_content, r_id",
    get_references_content(references),
)
def test_schema_validation(r_source, r_type, r_content, r_id, monkeypatch: pytest.MonkeyPatch):
    '''parse raw response and check its output schema'''
    monkeypatch.setattr("mutalyzer_retriever.parsers.json_ensembl._seq_from_rest", 
                        lambda _1, _2, _3, _4: _seq_from_rest(r_id))
    r_model = parser.parse(
        reference_content=r_content,
        reference_type=r_type,
        reference_source=r_source,
    )
    if r_source == "lrg":
        assert validate(r_model["annotations"]) is None
    else:
        assert validate(r_model) is None
