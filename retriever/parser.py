from .parsers import genbank, lrg, gff3


def get_reference_type(content):
    if content.startswith('<?xml version='):
        return 'lrg'
    elif content.startswith('LOCUS'):
        return 'genbank_ncbi'


def parse(reference_content, reference_type=None):

    if reference_type is None:
        reference_type = get_reference_type(reference_content)

    if reference_type == 'lrg':
        model = lrg.parse(reference_content)
    elif reference_type == 'genbank':
        model = genbank.parse(reference_content)
    elif reference_type == 'gff3':
        model = gff3.parse(reference_content)
    else:
        return None

    return model
