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


def extract_model(features, indent):
    for feature in features:
        print('{}id = {}'.format(indent, feature.id))
        print('{}type = {}'.format(indent, feature.type))
        if feature.sub_features:
            extract_model(feature.sub_features, indent + ' ')



def parse(feature_id):
    """
    Entry point for the GFF parser.
    """
    URL = 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi'
    payload = {'db': 'nuccore',
               'report': 'gff3',
               'id': feature_id}
    response = requests.get(URL, params=payload)

    gff_parser = GFFParser()
    gff = gff_parser.parse(io.StringIO(response.text))

    for rec in gff:
        extract_model(rec.features, '')


def main():
    """
    Command-line interface to gff parser.
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(dest='id',
                        help='NCBI feature ID.')

    args = parser.parse_args()

    parse(args.id)


if __name__ == '__main__':
    main()
