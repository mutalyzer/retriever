from schema import Schema, And, Or, Optional

transcript_provider_schema = Schema({
    "name": And(str, lambda s: s.strip() != ""),
    Optional("transcript_id"): And(str, lambda s: s.strip() != ""),
    Optional("protein_id"): And(str, lambda s: s.strip() != ""),
    Optional("description"): And(str, lambda s: s.strip() != ""),
})

transcript_schema = Schema({
    Optional("tag"): And(str, lambda s: s.strip() != ""),
    "providers": [transcript_provider_schema],
})

gene_schema = Schema({
    Optional("hgnc_id"): And(str, lambda s: s.strip() != ""),
    Optional("ensembl_id"): And(str, lambda s: s.strip() != ""),
    Optional("ncbi_rsg_id"): And(str, lambda s: s.strip() != ""),
    "name": And(str, lambda s: s.strip() != ""),
    "transcripts": [transcript_schema],
    Optional("description"): And(str, lambda s: s.strip() != "")
})

assembly_schema = Schema({
    "name": And(str, lambda s: s.strip() != ""),
    "accession": And(str, lambda s: s.strip() != ""),
    Optional("description"): And(str, lambda s: s.strip() != "")
})

related_schema = Schema({
    "assemblies": [assembly_schema],
    "genes": [gene_schema]
})

