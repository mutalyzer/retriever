from .sources import ncbi, ensembl, lrg
from retriever import parser


def fetch_annotations(reference_id, reference_type=None):
    annotations, reference_type = \
        ncbi.fetch_annotations(reference_id, reference_type)
    if annotations is not None:
        return annotations, reference_type, 'ncbi'
    annotations = lrg.fetch_lrg(reference_id)
    if annotations is not None:
        return annotations, 'lrg', 'lrg'
    annotations, reference_type = \
        ensembl.get_annotations(reference_id, reference_type)
    return annotations, reference_type, 'ensembl'


def fetch_sequence(reference_id, reference_source=None):
    if reference_source is None:
        sequence = ncbi.get_sequence(reference_id)
        if sequence is not None:
            return sequence
        sequence = ensembl.get_sequence(reference_id)
        return sequence
    else:
        if reference_source == 'ncbi':
            return ncbi.get_sequence(reference_id)
        elif reference_source == 'ensembl':
            return ensembl.get_sequence(reference_id)


def retrieve(reference_id, reference_source=None, reference_type=None,
             size_off=True, parse=False):
    """
    Main retriever entry point. Identifies and calls the appropriate specific
    retriever methods.

    :param reference_type: The type of the reference: gff, genbank, or lrg.
    :param reference_source: The source of the reference, e.g., ncbi, ensembl.
    :arg bool parse: Flag for parsing or not the reference.
    :arg str reference_id: The id of the reference.
    :arg bool size_off: Flag for the maximum sequence length.
    :return: The reference raw content and its type.
    :rtype: tuple
    """
    annotations = model = sequence = None
    if reference_source is None:
        annotations, reference_type, reference_source = \
            fetch_annotations(reference_id, reference_type)
    elif reference_source == 'ncbi':
        if reference_type is None or reference_type == 'gff':
            annotations = ncbi.get_gff(reference_id)
            reference_type = 'gff'
        elif reference_type == 'genbank':
            annotations = ncbi.get_genbank(reference_id, not size_off)
            reference_type = 'genbank'
        elif reference_type == 'sequence':
            return fetch_sequence(reference_id, reference_source)
    elif reference_source == 'ensembl':
        if reference_type is None or reference_type == 'gff':
            annotations = ensembl.get_gff(reference_id)
            reference_type = 'gff'
        elif reference_type == 'json':
            annotations = ensembl.get_json(reference_id)
    elif reference_source == 'lrg':
        annotations = lrg.fetch_lrg(reference_id)

    if parse:
        model = parser.parse(annotations, reference_type)
        if reference_type is 'gff':
            sequence = fetch_sequence(reference_id, reference_source)
            return {'model': model,
                    'sequence': sequence,
                    'source': reference_source}
        else:
            model.update({'source': reference_source})
            return model
    else:
        return annotations
