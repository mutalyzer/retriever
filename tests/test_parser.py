import pytest
from pathlib import Path

from retriever import parser
from retriever.schema_validation import validate


def get_references():
    references = {'gff3': ['NM_078467.2',
                           'NM_152263.2',
                           'NM_152263.3',
                           'NG_012337.1',
                           'NR_002196.2',
                           'L41870.1',
                           'NG_007485.1',
                           'ENSG00000147889'],
                  'lrg': ['LRG_11']}
    references_content = []
    for reference_type in references.keys():
        for reference_id in references[reference_type]:
            path_gb = Path(Path(__file__).parent) / 'data' / '{}.{}'.format(
                reference_id, reference_type)
            with path_gb.open() as f:
                reference_content = f.read()
            references_content.append(
                (reference_id, reference_type, reference_content))
    return references_content


@pytest.mark.parametrize('reference_id, reference_type, reference_content',
                         get_references())
def test_schema_validation(reference_id, reference_type, reference_content):
    reference_model = parser.parse(reference_content, reference_type)
    if reference_type == 'lrg':
        assert validate(reference_model['model']) is None
    else:
        assert validate(reference_model) is None


