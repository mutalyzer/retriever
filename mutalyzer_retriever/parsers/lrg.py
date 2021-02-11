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

from xml.dom import minidom

from Bio.Seq import Seq

from ..util import make_location


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
    return None


def _attr2dict(attr):
    """
    Create a dictionary from the attributes of an XML node

    :arg attr: A minidom node.

    :return: A dictionary with pairing of node-attribute names and values.
    Integer string values are converted to integers.
    :rtype: dictionary
    """
    attr_dict = {}
    for key, value in attr.items():
        if value.isdigit():
            value = int(value)
        attr_dict[key] = value
    return attr_dict


def _get_location(data, coord_system=None, recursive=False):
    """
    Get attributes from descendent <coordinates> element as a dictionary. If
    more than one <coordinates> element is found, we have a preference for the
    one with 'coord_system' attribute equal to the `coord_system` argument, if
    defined.
    """
    result = None
    coordinates = data.getElementsByTagName("coordinates")
    for coordinate in coordinates:
        attributes = _attr2dict(coordinate.attributes)
        if result and coord_system and attributes.get("coord_system") != coord_system:
            continue
        result = attributes
        if not recursive:
            break
    return {
        "type": "range",
        "start": {"type": "point", "position": int(result["start"]) - 1},
        "end": {"type": "point", "position": int(result["end"])},
        "strand": int(result["strand"]),
    }


def _get_translation_exception(cds):
    output = []
    for t_e in cds.getElementsByTagName("translation_exception"):
        output.append(
            {
                "location": make_location(t_e.getAttribute("codon")),
                "amino_acid": _get_content(t_e, "sequence"),
            }
        )
    if output:
        return {"translation_exception": {"exceptions": output}, "coordinate_system": "p"}


def _get_transcripts(section):
    """
    Extracts the transcripts present in the (fixed) section of the LRG file.

    :param section: (fixed) section of the LRG file
    :return: list of transcripts (GenRecord.Locus)
    """
    lrg_id = _get_content(section, "id")

    transcripts = []
    for transcript_data in section.getElementsByTagName("transcript"):
        transcript = {
            "id": transcript_data.getAttribute("name"),
            "location": _get_location(transcript_data, lrg_id),
        }

        # Get the exons.
        exons = []
        for exon_data in transcript_data.getElementsByTagName("exon"):
            exons.append(
                {
                    "type": "exon",
                    "id": exon_data.getAttribute("label"),
                    "location": _get_location(exon_data, lrg_id),
                }
            )
        transcript["features"] = exons

        # Get the CDS.
        transcript_type = "ncRNA"
        for cds_id, source_cds in enumerate(
            transcript_data.getElementsByTagName("coding_region")
        ):
            if cds_id > 0:
                # Todo: For now, we only support one CDS per transcript and
                #   ignore all others. This should be discussed.
                continue
            feature = {
                "type": "CDS",
                "id": source_cds.getElementsByTagName("translation")[0].getAttribute(
                    "name"
                ),
                "location": _get_location(source_cds, lrg_id),
            }
            translation_exception = _get_translation_exception(source_cds)
            if translation_exception:
                feature["qualifiers"] = translation_exception

            transcript["features"].append(feature)
            transcript_type = "mRNA"

        transcript["type"] = transcript_type

        transcripts.append(transcript)
    return transcripts


def _get_gene(fixed):
    """
    Construct the gene reference model.

    :param fixed: The fixed section of the LRG XML file.
    :return: Corresponding loci reference model.
    """
    gene = {
        "type": "gene",
        "id": _get_content(fixed, "hgnc_id"),
        "features": _get_transcripts(fixed),
    }
    return gene


def parse(content):
    """
    Parses an LRG <xml> formatted string and calls the appropriate methods to
    create and return the defined reference model.

    :arg bytes content: LRG file content.
    :return: Corresponding dictionary model.
    """

    # Extract the fixed section.
    data = minidom.parseString(content)
    fixed = data.getElementsByTagName("fixed_annotation")[0]

    # Get the sequence from the fixed section
    sequence = Seq(_get_content(fixed, "sequence"))

    annotations = {
        "type": "record",
        "id": _get_content(data, "id"),
        "location": {
            "type": "range",
            "start": {"type": "point", "position": 0},
            "end": {"type": "point", "position": len(sequence)},
        },
        "qualifiers": {
            "organism": _get_content(data, "organism"),
            "sequence_source": _get_content(data, "sequence_source"),
            "creation_date": _get_content(data, "creation_date"),
            "hgnc_id": _get_content(data, "hgnc_id"),
            "mol_type": _get_content(data, "mol_type"),
        },
        "features": [_get_gene(fixed)],
    }

    return {"annotations": annotations, "sequence": {"seq": str(sequence)}}
