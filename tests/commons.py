import json
from pathlib import Path

import pytest

from mutalyzer_retriever.configuration import settings
from mutalyzer_retriever.retriever import NoReferenceError

API_BASE = settings["ENSEMBL_API"]
API_BASE_GRCH37 = settings["ENSEMBL_API_GRCH37"]
TARK_API_BASE = settings["ENSEMBL_TARK_API"]


API_BASE_MAP = {
    "ENSG00000147889": {"version": 18, "species": "homo_sapiens"},
    "ENSMUSG00000022346": {"version": 18, "species": "mus_musculus"},
    "ENST00000304494": {"version": 10, "species": "homo_sapiens"},
    "ENST00000000000": {"version": 20, "species": "homo_sapiens"},
    "ENSG00000284190": {"version": 1, "species": "homo_sapiens"},
    "ENSG00000141510": {"version": 21, "species": "homo_sapiens"},
}
API_BASE_GRCH37_MAP = {
    "ENSG00000147889": {"version": 12, "species": "homo_sapiens"},
    "ENST00000304494": {"version": 5, "species": "homo_sapiens"},
    "ENST00000000000": {"version": 6, "species": "homo_sapiens"},
}

TARK_API_BASE_MAP = {
    "ENST00000304494": {"GRCH38_version": [10, 9, 8, 7, 6], "GRCH37_version": [5]},
    "ENST00000000000": {"GRCH38_version": [20, 19, 18], "GRCH37_version": [6, 5]},
}


@pytest.fixture(autouse=True)
def patch_retriever(monkeypatch):
    """retrieve all monkeypath"""
    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.fetch_gff3", _fetch_gff3)
    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl._get_reference_information",
        _get_reference_information,
    )
    monkeypatch.setattr(
        "mutalyzer_retriever.sources.ensembl._get_tark_versions", _get_tark_versions
    )
    monkeypatch.setattr("mutalyzer_retriever.sources.ensembl.fetch_json", _fetch_json)
    monkeypatch.setattr("mutalyzer_retriever.retriever.retrieve_raw", _retrieve_raw)


def _fetch_json(r_id, r_version, api_base, assembly, timeout):
    if api_base == TARK_API_BASE:
        return _get_content(f"data/{r_id}.{r_version}.tark_raw.model.json")


def _get_tark_versions(r_id, api_base, timeout=1):
    if api_base == TARK_API_BASE and r_id in TARK_API_BASE_MAP:
        return (
            TARK_API_BASE_MAP[r_id]["GRCH38_version"],
            TARK_API_BASE_MAP[r_id]["GRCH37_version"],
        )


def _fetch_gff3(feature_id, api_base, timeout=1):
    if api_base == API_BASE_GRCH37:
        return _get_content(
            f"data/{feature_id}.{API_BASE_GRCH37_MAP[feature_id]['version']}.gff3"
        )
    return _get_content(f"data/{feature_id}.gff3")


def _get_reference_information(r_id, api_base, timeout=1):
    if api_base == API_BASE and r_id in API_BASE_MAP:
        return API_BASE_MAP[r_id]
    if api_base == API_BASE_GRCH37 and r_id in API_BASE_GRCH37_MAP:
        return API_BASE_GRCH37_MAP[r_id]


def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    try:
        with open(str(data_file), "r") as file:
            content = file.read()
    except FileNotFoundError as exc:
        raise NoReferenceError({}, []) from exc
    return content


def _retrieve_raw(
    r_id,
    r_source=None,
    r_type=None,
    size_off=True,
    configuration_path=None,
    timeout=1,
):
    if r_type == "fasta":
        return _get_content("data/" + r_id + ".fasta"), "fasta", "ncbi"
    elif r_id.startswith("LRG_"):
        return _get_content("data/" + r_id + ".lrg"), "lrg", "lrg"
    elif r_type == "json":
        return (
            json.loads(_get_content("data/" + r_id + ".tark_raw.json")),
            "json",
            "ensembl_tark",
        )
    else:
        return _get_content("data/" + r_id + ".gff3"), "gff3", "ncbi"


references = {
    "ncbi": {
        "gff3": [
            "NM_078467.2",  # mRNA record, no mRNA feature: gene exons go to the mRNA
            "NM_152263.2",
            "NM_152263.3",
            "NM_000077.4",
            "NM_002001.2",
            "NG_012337.1",
            "NR_002196.2",  # gene exons (spliced) preferred over the lnc_RNA exon
            "L41870.1",
            "NG_007485.1",
            "NC_012920.1",
            "NG_009930.1",  # transl_except repeated on all the CDS slices
            "AA010203.1",
            "NP_060665.3",
            # Gene and parentless exons dropped; mRNAs with no gene; 2 split CDSs
            "D64137.1",
            # Gene and parentless exons dropped; mRNAs with no gene; 3 split CDSs
            "AB006684.1",
            "NM_004152.3",  # CDS slices with ribosomal slippage, kept apart
            "7",
            "M65131.1",  # mRNA record with gene-mRNA-exon/CDS nesting
            "XR_948219.2",  # gene-lnc_RNA-exon nesting
            "NR_023343.1",  # snRNA, ncRNA wrapper
            "NC_015338.1",
            "NR_003286.4",  # mol_type=rRNA, ncRNA wrapper
            "NR_114042.1",  # rRNA with no gene
            "NC_003468.2",  # mol_type=genomic RNA, not a transcript record
            # mRNA with no transcript_id; two LAT ncRNAs, one spanning the
            # origin of the circular genome; UL15 spliced CDS with no mRNA
            "NC_001806.2",
            "NC_011507.2",  # gene with no Name, id from its description (VP1)
            "LC832987.1",  # NS1 and UP1 genes split over two lines sharing an ID
            "PZ706538.1",  # HLA-A, no mRNA: gene exons (with UTRs) kept, split CDS
            "PZ290550.1",  # COX1, no mRNA: gene exons miss a CDS slice, so not used
            "NG_008047.2",  # transcript (misc_RNA) NR_126035.2 as ncRNA; RNU7-1 snRNA
            "PV394717.1",  # LINC03128: lncRNA (spelling), no gene
            "PZ268039.1",  # lncRNA (spelling) under the PTEN gene
            "NR_163594.1",  # SSPOP pseudogene: transcript with no transcript_id
            "NR_047663.1",  # ACTN3 transcript (misc_RNA) with no transcript_id
            "NR_055330.1",  # piRNA (C. elegans)
            "NR_002716.1",  # RNU2-1 gene with no RNA and no exons
            "NR_002761.3",  # RNU2-2P gene with two identical snRNAs
            "NR_029493.1",  # MIR21 precursor_RNA, mature miRNAs not kept
            "NG_008076.3",  # MIR1289-1 precursor_RNA NR_031620.1 and its mature miRNA
            "X54156.1",  # p53 pre-mRNA (primary_transcript) with no gene, not kept
        ]
    },
    "ensembl_rest": {
        "gff3": [
            "ENSG00000147889",
            "ENST00000383925",
            "ENST00000304494",
            "ENSG00000198899",
            "ENSG00000284190",
            # TP53: transcript (protein_coding_CDS_not_defined) as mRNA
            "ENSG00000141510",
        ]
    },
    "ensembl_tark": {
        "json": [
            "ENST00000383925.1",
            "ENST00000383925",
            "ENST00000304494",
            "ENST00000304494.10",
        ]
    },
    "lrg": {"lrg": ["LRG_11", "LRG_417", "LRG_857"]},
}
