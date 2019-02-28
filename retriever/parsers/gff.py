"""
Module for gff files parsing.

We consider that:
    - gff files have 9 columns.
    - '#' is used to comment lines in the gff file, so we ignore those lines.
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


def add_to_existing(existing, new):
    if existing['location'] != new['location']:
        existing['location'].extend(new['location'])
    for k in new['attributes']:
        if existing['attributes'].get(k) is None:
            existing['attributes'][k] = new['attributes'][k]
            # print('strange thing here since key: {} not in existing.'.format(k))
            # print(json.dumps(existing, indent=2))
            # print(json.dumps(new, indent=2))
        else:
            if k == 'phase':
                existing['attributes'][k].extend(new['attributes'][k])
            elif existing['attributes'][k] != new['attributes'][k]:
                print('different value for key: {}'.format(k))
                print(' {} in existing'.format(existing['attributes'][k]))
                print(' {} in new'.format(new['attributes'][k]))


def parse(gff_file):
    """
    Import transcript mappings from an GFF file.
    """
    records = {}
    already_there = 0
    with open(gff_file, 'rb') as gff_content:
        for raw_record in read_gff_raw_records(gff_content):
            if raw_record['attributes']['ID'] in records:
                add_to_existing(records[raw_record['attributes']['ID']],
                                raw_record)
                already_there += 1
            else:
                records[raw_record['attributes']['ID']] = raw_record
    print("done")
    print("Total records", len(records))
    print("Already there", already_there)
    print(json.dumps(records, indent=2))
    # print(json.dumps(chromosomes, indent=2))


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
