import json
from pathlib import Path

import pytest

from mutalyzer_retriever import parser
from mutalyzer_retriever.schema_validation import validate

from .commons import references



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
                        / f"{r_id}.tark_raw.{r_type}"
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
                        id=f"{r_source}-{r_type}-{r_id}",
                    )
                )
    return r_contents


@pytest.mark.parametrize(
    "r_source, r_type, r_content",
    get_references_content(references),
)
def test_schema_validation(r_source, r_type, r_content):
    '''parse raw response and check its output schema'''
    r_model = parser.parse(
        reference_content=r_content,
        reference_type=r_type,
        reference_source=r_source,
    )
    if r_source == "lrg":
        assert validate(r_model["annotations"]) is None
    else:
        assert validate(r_model) is None
