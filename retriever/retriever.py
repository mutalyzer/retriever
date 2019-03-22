from .sources import ncbi, ensembl, lrg
from retriever import parser


def fetch_annotations(reference_id):
    annotations = ncbi.get_gff(reference_id)
    if annotations is None:
        annotations = ensembl.get_gff(reference_id)
    else:
        return annotations, 'gff_ncbi'
    if annotations is None:
        annotations = lrg.fetch_lrg(reference_id)


def retrieve(reference_id, reference_source=None, reference_type=None,
             size_on=True, parse=False):
    """
    Main retriever entry point. Identifies and calls the appropriate specific
    retriever methods (lrg or ncbi genbank so far).

    :param reference_type: The type of the reference, e.g., gff, genbank.
    :param reference_source: The source of the reference, e.g., ncbi, ensembl.
    :arg bool parse: Flag for parsing or not the reference.
    :arg str reference_id: The id of the reference.
    :arg bool size_on: Flag for the maximum sequence length.
    :return: The reference raw content and its type.
    :rtype: tuple
    """
    output = None
    if reference_source is None and reference_type is None:
        fetch_annotations(reference_id)

    if reference_source == 'ncbi':
        if reference_type == 'gff':
            annotations = ncbi.get_gff(reference_id)
            sequence = ncbi.get_sequence(reference_id)
            parser_type = 'gff_ncbi'
            output = {'annotations': annotations,
                      'sequence': sequence}
    elif reference_source == 'ensembl':
        if reference_type == 'gff':
            annotations = ensembl.get_gff(reference_id)
            sequence = ensembl.get_sequence(reference_id)
            parser_type = 'gff_ensembl'
        output = {'annotations': annotations,
                  'sequence': sequence}
        if reference_type == 'json':
            annotations = ensembl.get_json(reference_id)
            sequence = ensembl.get_sequence(reference_id)
            parser_type = 'json_ensembl'
            output = {'annotations': annotations,
                      'sequence': sequence}
    elif reference_source == 'lrg':
        annotations = lrg.fetch_lrg(reference_id)
        output = {'model': parser.parse(annotations, 'lrg')}
        parser_type = 'lrg'

    if parse:
        model = parser.parse(annotations, parser_type)
        output['model'] = model

    return output
