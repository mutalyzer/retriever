from schema import Schema, And, Or, Optional

transcript_provider_schema = Schema({
    "name": And(str, lambda s: s.strip() != ""),
    Optional("transcript_accession"): And(str, lambda s: s.strip() != ""),
    Optional("protein_accession"): And(str, lambda s: s.strip() != ""),
    Optional("description"): And(str, lambda s: s.strip() != ""),
})

transcript_schema = Schema({
    Optional("tag"): And(str, lambda s: s.strip() != ""),
    "providers": [transcript_provider_schema],
})

gene_provider_schema = Schema({
    "name": And(str, lambda s: s.strip() != ""),
    Optional("id"): And(str, lambda s: s.strip() != ""),
    Optional("accession"): And(str, lambda s: s.strip() != "")
})

gene_schema = Schema({
    Optional("hgnc_id"): And(str, lambda s: s.strip() != ""),
    "name": And(str, lambda s: s.strip() != ""),
    Optional("providers"): [gene_provider_schema],
    Optional("transcripts"): [transcript_schema],
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

