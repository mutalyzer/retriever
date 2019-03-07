"""
Module for gff files parsing.

GFF3 specifications:
- https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
- ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt

Notes:
    - GFF files have 9 columns, not explicitly mentioned in the file, tab
    delimited, with the following order: seqid, source, type, start, end,
    score, strand, phase, and attributes.
    - There can be multiple parents for one entry, e.g., for exons.
    - Multiple entries can have the same parent.
    - There are entries with no parents.
    - Multiple entries have the same,
    '#' is used for comments.

    - mRNA and gene ID fields in the attributes column are unique.
    - CDS ID fields in the attributes column are not unique. However, the CDS
    entries with the same ID are part of the same protein. They are split like
    this in the same manner as the exons are.
"""
import argparse
import json

# Column names. (Used as keys in the corresponding dictionary).
gff_columns = [
    'seqid',
    'source',
    'type',
    'start',
    'end',
    'score',
    'strand',
    'phase',
    'attributes'
]

keys = {
    'region': 'seqid',
    'gene': 'HGNC',
    'mRNA': 'transcript_id',
    'CDS': 'protein_id',
    'transcript': 'transcript_id',
}

# Separators
# ----------
# Used between gff columns.
col_sep = '\t'
# Used between different fields in the attributes column.
attr_sep = ';'
# Used to separate keys from values in an attribute field.
attr_key_value_sep = '='


def extract_attributes(attributes):
    attr_dict = dict((k.strip(), v.strip()) for k, v in
                     (item.split(attr_key_value_sep) for item in
                      attributes.strip().split(attr_sep)))
    return attr_dict


def read_gff_raw_records(gff_content):
    """
    Converts a gff line into a raw record dictionary.

    :param gff_content:
    """
    for line in gff_content:
        if line.startswith(b'#'):
            continue
        entry = dict(zip(gff_columns, line.decode().strip().split(col_sep)))

        # Convert positions to ints.
        entry['location'] = [int(entry['start']), int(entry['end'])]

        # Extract the attributes.
        entry['attributes'] = extract_attributes(entry['attributes'])
        entry['attributes'].update({
            'seqid': entry['seqid'],
            'source': entry['source'],
            'score': entry['strand'],
            'strand': entry['strand'],
            'phase': [entry['phase']],
        })

        yield {
            'type': entry['type'],
            'attributes': entry['attributes'],
            'location': entry['location'],
        }


def add_to_composed(records, entry_id, new):
    """
    Some loci are composed of multiple entries. This function deals with
    composing the final loci object.
    """
    existing = records[entry_id]

    if existing.get('parts'):
        existing['parts'].append({'attributes': new['attributes'],
                                  'location': new['location']})
    else:
        existing = {'parts': [{'attributes': existing['attributes'],
                               'location': existing['location']},
                              {'attributes': new['attributes'],
                               'location': new['location']}],
                    'type': existing['type']}
    records[entry_id] = existing


def add_children_to_parents(records):
    """

    """
    for record_id in records:
        if records[record_id].get('parts'):
            parent_id = records[record_id]['parts'][0]['attributes'].get('Parent')
        else:
            parent_id = records[record_id]['attributes'].get('Parent')
        if parent_id:
            if records[parent_id].get('children') is None:
                records[parent_id]['children'] = [record_id]
            else:
                records[parent_id]['children'].append(record_id)


def retrieve_gene(records, gene_name):
    """
    Retrieves a gene object from an GFF record object.
    """
    for record_id in records:
        if records[record_id].get('type') == 'gene':
            if records[record_id]['attributes'].get('gene') == gene_name:
                return records[record_id]


def retrieve_locus(records, locus_name, locus_type):
    """
    Retrieves a gene object from an GFF record object.
    """
    for record_id in records:
        if records[record_id].get('type') == locus_type:
            if records[record_id]['attributes'].get('Name') == locus_name:
                return records[record_id]


def parse(gff_file):
    """
    Entry point for the GFF parser.
    """
    records = {}
    with open(gff_file, 'rb') as gff_content:
        for raw_record in read_gff_raw_records(gff_content):
            if raw_record['attributes']['ID'] in records:
                add_to_composed(records, raw_record['attributes']['ID'],
                                raw_record)
            else:
                records[raw_record['attributes']['ID']] = raw_record
    add_children_to_parents(records)
    print(json.dumps(records, indent=2))
    # print(json.dumps(retrieve_gene(records, 'PIGR'), indent=2))
    # print(json.dumps(records['rna13638'], indent=2))
    # print(json.dumps(records['rna13637'], indent=2))


def main():
    """
    Command-line interface to gff parser.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(dest='gff_file',
                        help='Path towards gff file.')

    args = parser.parse_args()

    parse(args.gff_file)


if __name__ == '__main__':
    main()
