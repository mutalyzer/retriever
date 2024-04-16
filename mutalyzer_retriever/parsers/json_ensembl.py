import sys

import requests


def point(position: int):
    return {"type": "point", "position": position}


def location(start: int, end: int, strand=None):
    loc_dict = {
        "type": "range",
        "start": point(start),
        "end": point(end),
        "strand": strand,
    }
    return {k: v for k, v in loc_dict.items() if v is not None}


def _feature(raw_dict):
    """Converting a general tark sub-dictionary into internal model, only id and location info"""
    return {
        "id": raw_dict["stable_id"],
        "location": location(
            raw_dict["loc_start"] - 1, raw_dict["loc_end"], raw_dict["loc_strand"]
        ),
    }


def annotations(id, location, features):
    annotations_dict = {
        "id": id,
        "type": "record",
        "location": location,
        "features": features,
    }
    return annotations_dict


def _exons(tark_exons):
    """converting exons info from tark list into internal exon list"""
    exons = []
    for tark_exon in tark_exons:
        exon = _feature(tark_exon)
        exon["type"] = "exon"
        exons.append(exon)
    return exons


def _translation(tark_translations):
    """Converting translations per transcript from tark list into internal translation list"""
    translations = []
    for tark_translation in tark_translations:
        translation = _feature(tark_translation)
        translation["type"] = "CDS"
        translations.append(translation)
    return translations


def _transcript(tark_transcript, exon_features, translation_feature):
    """Converting transcript from tark list into internal transcript list"""
    transcript = {}
    transcript = _feature(tark_transcript)
    transcript["type"] = tark_transcript["biotype"]
    if transcript["type"] == "protein_coding":
        transcript["type"] = "mRNA"
    transcript["qualifiers"] = {
        "assembly_name": tark_transcript["assembly"],
        "version": str(tark_transcript["stable_id_version"]),
        "tag": "basic",
    }
    transcript["features"] = exon_features + translation_feature
    return [transcript]


def _gene(tark_gene, gene_feature):
    gene = {}
    gene = _feature(tark_gene)
    gene["type"] = "gene"
    gene["qualifiers"] = {
        "assembly_name": tark_gene["assembly"],
        "version": str(tark_gene["stable_id_version"]),
        "name": tark_gene["name"],
    }
    gene["features"] = gene_feature
    return [gene]


def _seq_from_rest(chr_No, strand, loc_start, loc_end):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/region/human/{chr_No}:{loc_start}..{loc_end}:{strand}?coord_system_version=GRCh38"
    r = requests.get(server + ext, headers={"Content-Type": "text/plain"})
    if not r.ok:
        raise NameError
    return r.text


def _sequence(tark_result):
    return {
        "seq": _seq_from_rest(
            tark_result["loc_region"],
            tark_result["loc_strand"],
            tark_result["loc_start"],
            tark_result["loc_end"],
        ),
        "description": " ".join(
            [
                f"{tark_result['stable_id']}.{str(tark_result['stable_id_version'])}",
                ":".join(
                    [
                        "chromosome",
                        tark_result["assembly"],
                        str(tark_result["loc_region"]),
                        str(tark_result["loc_start"]),
                        str(tark_result["loc_end"]),
                        str(tark_result["loc_strand"]),
                    ]
                ),
            ]
        ),
    }


def parse(tark_result):
    """convert the tark json response of one transcript into the retriever model json output"""

    exon_features = _exons(tark_result["exons"])

    # TODO: find examples of null or 2
    # one translations per transcript, somtimes null or 2 (strand -1/1)
    translation_features = _translation(tark_result["translations"])

    transcript_features = _transcript(tark_result, exon_features, translation_features)

    # sometimes multiple gene info with different version or different content (eg. name=null/sth)
    # sort the genes based on stable_id_version and take the highest version
    # In case of same stable_id_version, the one with "name" field is considered higher.
    genes = sorted(
        tark_result["genes"],
        key=lambda g: (g["stable_id_version"], 0 if g["name"] is None else 1),
    )
    tark_gene = genes[-1]
    gene_feature = _gene(tark_gene, gene_feature=transcript_features)

    return {
        "annotations": annotations(
            tark_result["loc_region"],
            location(tark_result["loc_start"] - 1, tark_result["loc_end"]),
            gene_feature,
        ),
        "sequence": _sequence(tark_result),
    }
