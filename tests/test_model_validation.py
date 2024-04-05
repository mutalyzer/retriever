from pathlib import Path
import json
import pytest
from .commons import references
from mutalyzer_retriever import parser
from mutalyzer_retriever.schema_validation import validate

#check if schema is right



def get_references_content(references):
    references_content = []
    for reference_source in references.keys():
        for reference_type in references[reference_source]:
            for reference_id in references[reference_source][reference_type]:
                if reference_type == "json":
                    path_gb = (
                        Path(Path(__file__).parent)
                        / "data"
                        / "{}.tark_raw.{}".format(reference_id, reference_type)
                    )
                    reference_content = json.loads(path_gb.open().read())
                else:
                    path_gb = (
                                        Path(Path(__file__).parent)
                                        / "data"
                                        / "{}.{}".format(reference_id, reference_type)
                                    )
                    with path_gb.open() as f:
                        reference_content = f.read()
                references_content.append(
                    pytest.param(
                        reference_source,
                        reference_type,
                        reference_content,
                        reference_id,
                        id="{}-{}-{}".format(
                            reference_source, reference_type, reference_id
                        ),
                    )
                )
    return references_content

@pytest.mark.parametrize("reference_source, reference_type, reference_content, reference_id",get_references_content(references))
def test_schema_validation(reference_source, reference_type, reference_content, reference_id):
    reference_model = parser.parse(
        reference_content,
        reference_type=reference_type,
        reference_source=reference_source,
    )
    if reference_source == "lrg":
        assert validate(reference_model["annotations"]) is None
    else:
        assert validate(reference_model) is None
