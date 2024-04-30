import requests
from ..util import make_location, f_e


def _feature(raw_dict):
    """Convert a general tark sub-dictionary into our internal model.
       - only id and location info;
       - Tark locations are 1-based, our model is 0-based.
    """
    return {
        "id": raw_dict["stable_id"],
        "location": make_location(
            raw_dict["loc_start"] - 1, raw_dict["loc_end"], raw_dict.get("loc_strand")
        ),
    }


def _annotations(ref_id, location, features):
    return {
        "id": ref_id,
        "type": "record",
        "location": location,
        "features": features,
    }


def _exons(tark_exons):
    """Convert exons info from tark response list into internal exon list."""
    exons = []
    for tark_exon in tark_exons:
        exon = _feature(tark_exon)
        exon["type"] = "exon"
        exons.append(exon)
    return exons


def _translation(tark_translations):
    """Convert translations per transcript from tark list into internal translation list.
       - null for non-coding RNA in input, return an empty list;
       - one value for coding RNA in input, return a list of one item;
       - rarely multiple values for coding RNA in input with different versions,
         return a list of multiple items.
    """
    translations = []
    for tark_translation in tark_translations:
        translation = _feature(tark_translation)
        translation["type"] = "CDS"
        translations.append(translation)
    return translations


def _transcript(tark_transcript, exon_features, translation_feature):
    """Convert transcript from tark list into internal transcript list.
       - Tark has RNA type as protein_coding, change to internal RNA type mRNA.
    """
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
    """Convert gene info from tark list into internal gene list."""
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


def _seq_from_rest(assembly, chr_idx, strand, loc_start, loc_end, timeout=1):
    """Retrieve sequence from ensembl Rest API."""
    if assembly == "GRCh38":
        server = "https://rest.ensembl.org"
    elif assembly == "GRCh37":
        server = "https://grch37.rest.ensembl.org"
    else:
        raise NameError("Unsupported assembly {assembly}")
    ext = f"/sequence/region/human/{chr_idx}:{loc_start}..{loc_end}:{strand}?"
    r = requests.get(
        server + ext, headers={"Content-Type": "text/plain"}, timeout=timeout
    )
    if not r.ok:
        raise NameError
    return r.text


def _sequence(tark_result):
    return {
        "seq": _seq_from_rest(
            tark_result["assembly"],
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


def parse(tark_results):
    """Convert the Tark json response into the retriever model json output.
       - take the latest version from Tark response if no specific version required;
       - for genes, take the latest version with "name" field in case of same stable ID
    """
    tark_results = tark_results.get("results")
    if tark_results:
        tark_result = tark_results[-1]
    else:
        raise NameError(f_e("ensembl tark", e=None, extra="returns no results"))

    exon_features = _exons(tark_result["exons"])

    translation_features = _translation(tark_result["translations"])

    transcript_features = _transcript(tark_result, exon_features, translation_features)

    genes = sorted(
        tark_result["genes"],
        key=lambda g: (g["stable_id_version"], 0 if g["name"] is None else 1),
    )
    tark_gene = genes[-1]
    gene_feature = _gene(tark_gene, gene_feature=transcript_features)

    return {
        "annotations": _annotations(
            tark_result["loc_region"],
            make_location(tark_result["loc_start"] - 1, tark_result["loc_end"]),
            gene_feature,
        ),
        "sequence": _sequence(tark_result),
    }
