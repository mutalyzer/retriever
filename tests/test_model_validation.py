from pathlib import Path

import pytest

from mutalyzer_retriever import parser
from mutalyzer_retriever.schema_validation import validate


def get_references_content(references):
    references_content = []
    for reference_source in references.keys():
        for reference_type in references[reference_source]:
            for reference_id in references[reference_source][reference_type]:
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
                        id="{}-{}-{}".format(
                            reference_source, reference_type, reference_id
                        ),
                    )
                )
    return references_content


@pytest.mark.parametrize(
    "reference_source, reference_type, reference_content",
    get_references_content(
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
def test_schema_validation(reference_source, reference_type, reference_content):
    reference_model = parser.parse(
        reference_content,
        reference_type=reference_type,
        reference_source=reference_source,
    )
    if reference_source == "lrg":
        assert validate(reference_model["annotations"]) is None
    else:
        assert validate(reference_model) is None
