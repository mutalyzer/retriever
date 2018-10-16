from .get_lrg import fetch_lrg
from .get_ncbi import fetch_ncbi


def retrieve(reference_id, size_on=True):
    """
    Main retriever entry point. Identifies and calls the appropriate specific
    retriever methods (lrg or ncbi genbank so far).

    :param reference_id: the id of the reference
    :param size_on: flag for the maximum sequence length
    :return: reference content
    """
    if 'LRG' in reference_id:
        return fetch_lrg(reference_id, size_on)
    else:
        return fetch_ncbi(reference_id, size_on)