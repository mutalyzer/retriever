from schema import Schema, And, Or, Optional

gene_provider_schema = Schema({
    "name": And(str, lambda s: s.strip() != ""),
    Optional("accession"): Or(None, And(str, lambda s: s.strip() != ""))
})

transcript_provider_schema = Schema({
    "name": And(str, lambda s: s.strip() != ""),
    "transcript_id": Or(None, And(str, lambda s: s.strip() != "")),
    "protein_id": Or(None, And(str, lambda s: s.strip() != ""))
})

transcript_schema = Schema({
    Optional("tag"): Or(str, None),
    "providers": [transcript_provider_schema]
})

gene_schema = Schema({
    "hgnc_id": And(str, lambda s: s.strip() != ""),
    "name": And(str, lambda s: s.strip() != ""),
    "providers": [gene_provider_schema],
    "transcripts": [transcript_schema],
    Optional("comment"): Or(None, And(str, lambda s: s.strip() != "")),
})

assembly_schema = Schema({
    "assembly_name": And(str, lambda s: s.strip() != ""),
    "accession": And(str, lambda s: s.strip() != "")
})

related_schema = Schema({
    "assemblies": [assembly_schema],
    "genes": [gene_schema]
})

query_schema = Schema({
    "query": And(str, lambda s: s.strip() != ""),
    "organism": And(str, lambda s: s.strip() != ""),
    "moltype": And(str, lambda s: s.strip() != ""),
    "accession": And(str, lambda s: s.strip() != ""),
    "related": related_schema

})
