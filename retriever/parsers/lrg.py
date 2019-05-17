"""
LRG file parser.

An LRG file is an XML formatted file and consists of a fixed and updatable
section. The fixed section contains a DNA sequence and for that sequence a
number of transcripts. The updatable region could contain all sorts of
annotation for the sequence and transcripts. It can also contain additional
(partial) transcripts and mapping information.

More information on LRG files:
    https://www.lrg-sequence.org/faq/
    http://ftp.ebi.ac.uk/pub/databases/lrgex/docs/LRG_XML_schema_documentation_1_9.pdf
    http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG.rnc
    http://ftp.ebi.ac.uk/pub/databases/lrgex/docs/LRG.pdf

This module is based on the result of the minidom xml parser.

NOTE: A strong alternative to the minidom parser would be lxml, which was
already employed by Mutalyzer in other circumstances.
"""

import xml.dom.minidom
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from ..reference import Position, Locus, Reference


def _get_content(data, refname):
    """
    Return string-content of an XML text node.

    :arg data: A minidom object.
    :arg refname: The name of a member of the minidom object.

    :return: The content of the textnode or an emtpy string.
    :rtype: string
    """
    temp = data.getElementsByTagName(refname)
    if temp:
        return temp[0].lastChild.data
    return ""


def _attr2dict(attr):
    """
    Create a dictionary from the attributes of an XML node

    :arg attr: A minidom node.

    :return: A dictionary with pairing of node-attribute names and values.
    Integer string values are converted to integers.
    :rtype: dictionary
    """
    ret = {}
    for key, value in attr.items():
        if value.isdigit():
            value = int(value)
        ret[key] = value
    return ret


def _get_coordinates(data, system=None):
    """
    Get attributes from descendent <coordinates> element as a dictionary. If
    more than one <coordinates> element is found, we have a preference for the
    one with 'coord_system' attribute equal to the `system` argument, if
    defined.
    """
    result = None
    coordinates = data.getElementsByTagName('coordinates')
    for coordinate in coordinates:
        attributes = _attr2dict(coordinate.attributes)
        if result and system and attributes.get('coord_system') != system:
            continue
        result = attributes
    return {'start': int(result['start']),
            'end': int(result['end']),
            'strand': int(result['strand'])}


def _get_gene_name(section):
    """
    Extract the gene name from the LRG record updatable section.

    NOTE: It is necessary to use the updatable section since there is no
    other way to identify the main gene directly from the LRG file.
    A possibility would be to use the HGNC id with some external service.
    Another way would be to make use of the special file with genes to LRG:
    http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_transcripts_GRCh38.txt

    :arg section: (updatable) section of the LRG file
    :return: gene name present under the lrg annotation set
    """
    gene_name = ''
    annotation_nodes = section.getElementsByTagName('annotation_set')
    for anno in annotation_nodes:
        if anno.getAttribute('type') == 'lrg':
            gene_name = _get_content(anno, 'lrg_locus')
    return gene_name


def _get_gene(fixed, updatable):
    """
    Create the gene locus.

    :param fixed: The fixed section of the LRG XML file.
    :param updatable: The updatable section of the LRG XML file.
    :return: Corresponding gene locus.
    """

    return {'id': _get_content(fixed, 'hgnc_id'),
            'type': 'gene'}


def _get_transcripts(section):
    """
    Extracts the transcripts present in the (fixed) section of the LRG file.

    :param section: (fixed) section of the LRG file
    :return: list of transcripts (GenRecord.Locus)
    """
    lrg_id = _get_content(section, 'id')

    transcripts = []
    for tdata in section.getElementsByTagName('transcript'):
        transcript = {'id': tdata.getAttribute('name')}
        transcript.update(_get_coordinates(tdata, lrg_id))

        # Get the transcript exons.
        exons = []
        for exon_data in tdata.getElementsByTagName('exon'):
            exon = {'type': 'exon'}
            exon.update(_get_coordinates(exon_data, lrg_id))
            exons.append(exon)

        transcript['sub_features'] = exons

        # Get the CDSes of the transcript.
        for cds_id, source_cds in enumerate(
                tdata.getElementsByTagName("coding_region")):
            if cds_id > 0:
                # Todo: For now, we only support one CDS per transcript and
                #   ignore all others. This should be discussed.
                continue
            cds = {'type': 'cds',
                   'id': source_cds.getElementsByTagName(
                       'translation')[0].getAttribute('name')}
            cds.update(_get_coordinates(source_cds, lrg_id))
            transcript['sub_features'].append(cds)

        transcripts.append(transcript)
    return transcripts


def _get_loci(fixed, updatable):
    """
    Construct the loci reference model.

    :param fixed: The fixed section of the LRG XML file.
    :param updatable: The updatable section of the LRG XML file.
    :return: Corresponding loci reference model.
    """
    gene = _get_gene(fixed, updatable)

    gene['sub_features'] = _get_transcripts(fixed)

    return gene


def parse(content):
    """
    Parses an LRG <xml> formatted string and calls the appropriate methods to
    create and return the defined reference model.

    :arg bytes content: Content of LRG file

    :return: Corresponding reference instance.
    :rtype: Reference
    """

    # Extract the fixed and updatable section.
    data = xml.dom.minidom.parseString(content)
    fixed = data.getElementsByTagName('fixed_annotation')[0]
    updatable = data.getElementsByTagName('updatable_annotation')[0]

    info = {
        'organism': _get_content(data, 'organism'),
        'sequence_source': _get_content(data, 'sequence_source'),
        'creation_date': _get_content(data, 'creation_date'),
        'molType': 'g'}

    # Get the sequence from the fixed section
    sequence = Seq(_get_content(fixed, 'sequence'), IUPAC.unambiguous_dna)

    features = _get_loci(fixed, updatable)

    return {'model':
                {'info': info,
                 'features': features},
            'sequence': str(sequence)}
