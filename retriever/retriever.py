from .sources import ncbi, ensembl
from retriever import parser


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
    if reference_source == 'ncbi':
        if reference_type == 'gff':
            annotations = ncbi.get_gff(reference_id)
            sequence = ncbi.get_sequence(reference_id)
            output = {'annotations': annotations,
                      'sequence': sequence}
            if parse:
                model = parser.parse(annotations, 'gff_ncbi')
                output['model'] = model
    elif reference_source == 'ensembl':
        if reference_type == 'gff':
            annotations = ensembl.get_gff(reference_id)
            # sequence = ncbi.get_sequence(reference_id)
            output = {'annotations': annotations},
                      # 'sequence': sequence}
            if parse:
                model = parser.parse(annotations, 'gff_ensembl')
                output['model'] = model

    return output
