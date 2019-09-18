from retriever.util import make_request
from Bio import SeqIO
from io import StringIO
from urllib.error import HTTPError
import json


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


def get_sequence_details(feature_id):
    url = 'https://rest.ensembl.org/lookup/id/{}'.format(feature_id)
    headers = {'Content-Type': 'application/json'}
    response = json.loads(make_request(url, headers=headers))
    return response['start'], response['end'], response['species'], \
           response['seq_region_name']


def get_sequence(feature_id):
    start, end, species, seq_region_name = get_sequence_details(feature_id)
    url = 'https://rest.ensembl.org/sequence/region/{}/{}:{}..{}'.format(
        species, seq_region_name, start, end)
    headers = {'Content-Type': 'application/json'}
    return json.loads(make_request(url, headers=headers))


def get_annotations(reference_id, reference_type):
    if reference_type in [None, 'gff3']:
        return get_gff(reference_id), 'gff3'
    if reference_type == 'json':
        return get_json(reference_id), 'json'
    if reference_type == 'genbank':
        return None, 'genbank'
    return None, None
