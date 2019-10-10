"""
Module for gff files parsing.

GFF3 specifications:
- Official:
  - [1] https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
- NCBI:
  - [2]: ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt
  - https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
- Ensembl:
  - ftp://ftp.ensembl.org/pub/release-98/gff3/homo_sapiens/README

Download ftp:
- ftp://ftp.ensembl.org/pub/current_gff3/homo_sapiens/
- ftp://ftp.ensembl.org/pub/release-98/gff3/homo_sapiens/

Sequence Ontology Feature Annotation (SOFA):
- http://www.sequenceontology.org/so_wiki/index.php/FAQ

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
from Bio.SeqFeature import SeqFeature
import io
from ..util import make_location

CONSIDERED_TYPES = ['region', 'gene', 'mRNA', 'exon', 'CDS', 'lnc_RNA']
QUALIFIERS = {'gene': {'Name': 'name',
                       'gene_synonym': 'synonym'},
              'region': {'mol_type': 'mol_type',
                         'Is_circular': 'is_circular',
                         'transl_table': 'transl_table'}}
SO_IDS = {'gene': 'SO:0000704',
          'mRNA': 'SO:0000234',
          'ncRNA': 'SO:0000655',
          'exon': 'SO:0000147',
          'CDS': 'SO:0000316'}


def _get_feature_id(feature):
    if feature.type == 'gene':
        if feature.qualifiers.get('gene_id'):
            return feature.qualifiers['gene_id'][0]
        return feature.qualifiers['Name'][0]
    elif feature.type == 'mRNA':
        return feature.qualifiers['transcript_id'][0]
    elif feature.type == 'lnc_RNA':
        return feature.qualifiers['transcript_id'][0]
    elif feature.type == 'CDS':
        return feature.qualifiers['protein_id'][0]
    elif feature.type == 'exon':
        if feature.qualifiers.get('exon_id'):
            return feature.qualifiers['exon_id'][0]
        elif feature.id:
            return feature.id


def _combine_cdses(mrna):
    """
    Combine all the cds features into a single one.
    """
    positions = []
    exons = []
    for feature in mrna['features']:
        if feature['type'] == 'CDS':
            positions.append(feature['location']['start']['position'])
            positions.append(feature['location']['end']['position'])
        elif feature['type'] == 'exon':
            exons.append(feature)
    positions = sorted(positions)
    for feature in mrna['features']:
        if feature['type'] == 'CDS':
            feature['location']['start']['position'] = positions[0]
            feature['location']['end']['position'] = positions[-1]
            mrna['features'] = exons + [feature]
            return


def _get_qualifiers(feature):
    q = feature.qualifiers
    t = feature.type
    if feature.type in QUALIFIERS.keys():
        qs = {QUALIFIERS[t][k]: q[k][0] if len(q[k]) == 1 else q[k]
              for k in q.keys() if k in QUALIFIERS[t].keys()}
        if t == 'gene':
            if q.get('Dbxref'):
                for dbxref_entry in q['Dbxref']:
                    if 'HGNC' in dbxref_entry:
                        qs['HGNC'] = dbxref_entry.split(':')[-1]
        return qs


def _get_feature_type(feature):
    if feature.type in ['gene']:
        return 'gene'
    elif feature.type in ['mRNA']:
        return 'mRNA'
    elif feature.type in ['lnc_RNA']:
        return 'ncRNA'
    elif feature.type in ['exon']:
        return 'exon'
    elif feature.type in ['CDS']:
        return 'CDS'
    else:
        return feature.type


def _get_feature_model(feature, parent=None, skip=None):
    """
    Recursively get the model for a particular feature. If some sub features
    do not need to be included, specify them in the `skip` dictionary.

    The method to combine CDSes into a single feature is also employed.
    """
    if skip and parent and isinstance(parent, SeqFeature):
        if parent in skip.keys() and skip[parent] == feature.type:
            return
    if feature.type in CONSIDERED_TYPES:
        model = {'id': _get_feature_id(feature),
                 'type': _get_feature_type(feature),
                 'location': make_location(
                     feature.location.start,
                     feature.location.end,
                     feature.location.strand)}
        qualifiers = _get_qualifiers(feature)
        if qualifiers:
            model['qualifiers'] = qualifiers
        if feature.sub_features:
            model['features'] = []
            for sub_feature in feature.sub_features:
                sub_feature_model = _get_feature_model(sub_feature, feature,
                                                       skip)
                if sub_feature_model:
                    model['features'].append(sub_feature_model)
        if feature.type == 'mRNA':
            _combine_cdses(model)
        return model


def _get_region_model(features):
    """
    Multiple `region` features can be present in the file. According to
    the NCBI [2], the one that corresponds to the `source` feature that
    appears in a GenBank flatfile format can be identified by the
    `gbkey=Src` attribute and is the first feature row for every seqid.
    """
    for feature in features:
        if feature.type == 'region' and feature.qualifiers.get('gbkey'):
            if feature.qualifiers['gbkey'] == 'Src':
                return _get_feature_model(feature)


def _create_mrna_features(features, mrna_id):
    mrna_model = {'id': mrna_id}
    exon_positions = []
    for sub_feature in features[0]['features']:
        if sub_feature['type'] == 'exon':
            exon_positions.append(sub_feature['location']['start']['position'])
            exon_positions.append(sub_feature['location']['end']['position'])
    if exon_positions:
        mrna_model['location'] = make_location(sorted(exon_positions)[0, -1])
    print(json.dumps(mrna_model, indent=2))
    return mrna_model


def _create_mrna_model(record):
    features = []
    for feature in record.features:
        feature_model = _get_feature_model(feature, record, {'gene': 'exon'})
        if feature_model:
            features.append(feature_model)


def _create_record_model(record, source=None):
    """
    Our model follows the gene-mRNA-CDS/exon and gene-ncRNA-exon conventions.
    Annotations in GFF3 files also conform to this, with some exceptions:
    - `mol_type=mRNA` references (e.g., NM_/XM, NR_/XR), for which the RNA is
       missing: gene-(CDS)/exon. In this case we create the RNA.
    - There may be some floating exons attached directly to a gene. We do not
      add them to our model.
    """

    features = []
    # Consider first if the record is mRNA.
    region_model = _get_region_model(record.features)
    if region_model:
        if region_model.get('qualifiers'):
            if region_model['qualifiers']['mol_type'] == 'mRNA':
               features = _create_mrna_model(record)
    else:
        for feature in record.features:
            feature_model = _get_feature_model(feature, record,
                                               {'gene': 'exon'})
            if feature_model:
                features.append(feature_model)

    model = {'id': record.id,
             'location': make_location(
                 record.annotations['sequence-region'][0][2],
                 record.annotations['sequence-region'][0][1])}

    if features:
        model['features'] = features

    return model


def parse(gff_content, source=None):
    gff_parser = GFFParser()
    gff = gff_parser.parse(io.StringIO(gff_content))

    records = []
    for record in gff:
        records.append(_create_record_model(record, source))
    if len(records) >= 1:
        return records[0]
    # TODO: Decide what to do when there are multiple records.