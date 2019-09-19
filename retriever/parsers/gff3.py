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
from BCBio.GFF import GFFParser
import io
from ..util import make_location


def _get_qualifiers(feature):
    qualifiers = {}
    for qualifier in feature.qualifiers:
        if len(feature.qualifiers[qualifier]) == 1:
            qualifiers[qualifier] = feature.qualifiers[qualifier][0]
        else:
            qualifiers[qualifier] = feature.qualifiers[qualifier]
    return qualifiers


def _get_features(feature):
    model = {'id': feature.id,
             'type': feature.type,
             'qualifiers': _get_qualifiers(feature),
             'location': make_location(
                 feature.location.start,
                 feature.location.end,
                 feature.location.strand)}
    if feature.sub_features:
        model['features'] = []
        for sub_feature in feature.sub_features:
            model['features'].append(_get_features(sub_feature))
    return model


def parse(gff_content):
    gff_parser = GFFParser()
    gff = gff_parser.parse(io.StringIO(gff_content))

    records = []
    for rec in gff:
        record = {'type': 'top',
                  'id': rec.id,
                  'qualifiers': {'annotations': rec.annotations},
                  'location': make_location(
                      rec.annotations['sequence-region'][0][2],
                      rec.annotations['sequence-region'][0][1])}
        features = []
        for feature in rec.features:
            features.append(_get_features(feature))
        if features:
            record['features'] = features
        records.append(record)
    if len(records) >= 1:
        return records[0]
    # TODO: Decide what to do when there are multiple records.
