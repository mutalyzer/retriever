import requests
import sys


def filter_none(dicts):
    dicts = {k:v for k,v in dicts.items() if v is not None}
    return dicts


def point(position: int):
    return {
        "type": "point",
        "position": position
    }

def location(start: int, end: int, strand=None):
    loc_dict = {
        "type": "range",
        "start": point(start),
        "end": point(end),
        "strand": strand
        }
    loc_dict = filter_none(loc_dict)
    return loc_dict

def feature(id: str, feature_type: str, location: dict, qualifiers=None, children_features=None):
    if feature_type == "protein_coding":
        feature_type = "mRNA"
    feature_dict = {
        "id": id,
        "type": feature_type,
        "location": location,
        "qualifiers": qualifiers,
        "features": children_features,
        }
    feature_dict = filter_none(feature_dict)
    return feature_dict

def annotations(id, location, features):
    annotations_dict = {
        "id": id,
        "type": "record",
        "location": location,
        "features": features,
        }
    annotations_dict = filter_none(annotations_dict)
    return annotations_dict

   
def parse(tark_results):
    tark_results = sorted(tark_results, key=lambda g: (g["stable_id_version"]))
    tark_result = tark_results[-1]
    tark_exons = tark_result["exons"]
    features = []
    for tark_exon in tark_exons:
        exon = feature(
            tark_exon['stable_id'], 
            "exon", 
            location(tark_exon["loc_start"]-1, tark_exon["loc_end"],tark_exon['loc_strand']),
        )
        features.append(exon)

    #one translations per transcript, somtimes null or 2 (strand -1/1)
    tark_translations = tark_result["translations"]
    if tark_translations:
        for tark_translation in tark_translations:
            translations = feature(
                tark_translation["stable_id"],
                "CDS",
                location(tark_translation['loc_start']-1, tark_translation['loc_end'],tark_translation['loc_strand']),    
            )

            features.append(translations)

    rna = feature(
        tark_result["stable_id"],
        tark_result["biotype"],
        location(tark_result["loc_start"]-1, tark_result["loc_end"],tark_result["loc_strand"]),
        qualifiers={
            "assembly_name": tark_result["assembly"],
            "version": str(tark_result["stable_id_version"]),
            "tag": "basic"
        },
        children_features = features
    )

    #sometimes multiple gene info with different version or different content (eg. name=null/sth)
    # sort the genes based on stable_id_version and take the highest version
    # In case of same stable_id_version, the one with "name" field is considered higher.
    genes = tark_result["genes"]
    genes = sorted(genes, key=lambda g: (g["stable_id_version"], 0 if g["name"] is None else 1))
    tark_gene = genes[-1]    

    gene = feature(
        tark_gene["stable_id"],
        "gene", 
        location(tark_gene['loc_start']-1, tark_gene['loc_end'],tark_gene['loc_strand']),
        qualifiers={
            "assembly_name": tark_gene["assembly"],
            "version": str(tark_gene["stable_id_version"]),
            "name": tark_gene["name"],
        },
        children_features=[rna]
    )
    return  annotations(
                tark_result["loc_region"], 
                location(tark_result["loc_start"]-1, tark_result["loc_end"]),
                [gene]
            )

    