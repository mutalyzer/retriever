"""
Module for gff files parsing.

GFF3 specifications:
- Official:
  - https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
- NCBI:
  - ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt

Notes:
    - GFF files have 9 columns, not explicitly mentioned in the file, tab
    delimited, with the following order: seqid, source, type, start, end,
    score, strand, phase, and attributes.
    - According to the official specifications, there can be multiple parents
    for one entry, e.g., for exons. However, it seems like NCBI does not adhere
    to this practice.
    - Multiple entries can have the same parent.
    - There are entries with no parents.
    - '#' is used for comments.
    - mRNA and gene ID fields in the attributes column are unique.
    - CDS ID fields in the attributes column are not unique. However, the CDS
    entries with the same ID are part of the same protein. They are split like
    this in the same manner as the exons are.
"""
import argparse
from BCBio.GFF import GFFParser
import io
import requests
import json
from Bio import Entrez
from Bio import SeqIO
from retriever.validation import validate


def get_qualifiers(feature):
    qualifiers = {}
    for qualifier in feature.qualifiers:
        if len(feature.qualifiers[qualifier]) == 1:
            qualifiers[qualifier] = feature.qualifiers[qualifier][0]
        else:
            qualifiers[qualifier] = feature.qualifiers[qualifier]
    return qualifiers


def extract_raw_features(feature):
    model = {'id': feature.id,
             'type': feature.type,
             'qualifiers': get_qualifiers(feature),
             'start': {'position': feature.location.start},
             'end': {'position': feature.location.end},
             'orientation': feature.location.strand}
    if feature.sub_features:
        model['sub_features'] = []
        for sub_feature in feature.sub_features:
            model['sub_features'].append(extract_raw_features(sub_feature))
    return model


def get_raw_record(feature_id):
    URL = 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi'
    params = {'db': 'nuccore',
              'report': 'gff3',
              'id': feature_id}

    try:
        response = requests.get(URL, params=params)
        response.raise_for_status()
    except requests.exceptions.HTTPError as errh:
        print("HTTP Error:", errh)
    except requests.exceptions.ConnectionError as errc:
        print("Connection Error:", errc)
    except requests.exceptions.Timeout as errt:
        print("Timeout Error:", errt)
    except requests.exceptions.RequestException as err:
        print("Some other Error", err)
    else:
        gff_parser = GFFParser()
        gff = gff_parser.parse(io.StringIO(response.text))

        model = []
        for rec in gff:
            for feature in rec.features:
                model.append(extract_raw_features(feature))
        return model


def get_id(feature):
    if feature.type == 'CDS':
        return feature.qualifiers.get('protein_id')[0]
    elif feature.type in ['mRNA']:
        return feature.qualifiers.get('transcript_id')[0]
    elif feature.type == 'gene':
        for dbxref in feature.qualifiers['Dbxref']:
            if dbxref.startswith('HGNC'):
                return dbxref.split(':')[-1]
    else:
        return feature.id


def get_sequence(feature_id):

    Entrez.email = 'm.lefter@lumc.nl'

    handle = Entrez.efetch(db='nucleotide', id=feature_id, rettype='fasta')
    records = []
    for record in SeqIO.parse(handle, "fasta"):
        records.append(str(record.seq))
    handle.close()
    if len(records) == 1:
        return records[0]


def main():
    """
    Command-line interface to gff parser.
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(dest='id',
                        help='NCBI feature ID.')

    args = parser.parse_args()

    get_raw_record(args.id)


if __name__ == '__main__':
    main()
