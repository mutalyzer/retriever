from retriever.util import make_request
from Bio import SeqIO
from io import StringIO


def get_json(feature_id):
    url = 'https://rest.ensembl.org/lookup/id/{}'.format(feature_id)
    params = {'expand': 1}
    headers = {'Content-Type': 'application/json'}
    return make_request(url, params, headers)


def get_gff(feature_id):
    url = 'http://rest.ensembl.org/overlap/id/{}'.format(feature_id)
    params = {'feature': ['gene', 'transcript', 'cds']}
    headers = {'Content-Type': 'text/x-gff3'}
    return make_request(url, params, headers)


def get_sequence(feature_id):
    url = 'http://rest.ensembl.org/sequence/id/{}'.format(feature_id)
    headers = {'Content-Type': 'text/x-fasta'}
    handle = StringIO(make_request(url, headers=headers))

    records = []
    for record in SeqIO.parse(handle, "fasta"):
        records.append(str(record.seq))
    handle.close()
    if len(records) == 1:
        return records[0]
