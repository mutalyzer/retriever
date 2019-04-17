from retriever.util import make_request
from Bio import SeqIO
from io import StringIO
from urllib.error import HTTPError


def get_json(feature_id):
    url = 'https://rest.ensembl.org/lookup/id/{}'.format(feature_id)
    params = {'feature': ['gene', 'transcript', 'cds'], 'expand': 1}
    headers = {'Content-Type': 'application/json'}
    return make_request(url, params, headers)


def get_gff(feature_id):
    url = 'https://rest.ensembl.org/overlap/id/{}'.format(feature_id)
    params = {'feature': ['gene', 'transcript', 'cds', 'exon']}
    headers = {'Content-Type': 'text/x-gff3'}
    try:
        return make_request(url, params, headers)
    except HTTPError as err:
        print("HTTP error")


def get_sequence(feature_id):
    url = 'https://rest.ensembl.org/sequence/id/{}'.format(feature_id)
    headers = {'Content-Type': 'text/x-fasta'}
    handle = StringIO(make_request(url, headers=headers))

    records = []
    for record in SeqIO.parse(handle, "fasta"):
        records.append(str(record.seq))
    handle.close()
    if len(records) == 1:
        return records[0]


def get_annotations(reference_id, reference_type):
    if reference_type in [None, 'gff']:
        return get_gff(reference_id), 'gff'
    if reference_type == 'json':
        return get_json(reference_id), 'json'
    if reference_type == 'genbank':
        return None, 'genbank'
    return None, None
