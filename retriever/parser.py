from .parsers import genbank, lrg, gff


def get_reference_type(content):
    if content.startswith('<?xml version='):
        return 'lrg'
    elif content.startswith('LOCUS'):
        return 'genbank_ncbi'


def parse(reference, reference_type=None):

    if reference_type is None:
        reference_type = get_reference_type(reference)

    if reference_type == 'lrg':
        model = lrg.parse(reference)
    elif reference_type == 'genbank_ncbi':
        model = genbank.parse(reference)
    elif reference_type in ['gff_ncbi', 'gff_ensembl']:
        model = gff.get_raw_record(reference)
    else:
        return None

    return model
