from schema import Schema, And, Or, Use

# Allow dict with string keys and non-empty string values
assembly_accession_schema = {str: And(str, len)}

gene_reference_standard_schema = {
    "Gene_Reference": Or(And(str, len), None)
}

product_feature_schema = {
    "transcript": Or(And(str, len), None),
    "ensembl_transcript": Or(And(str, len), None),
    "protein": Or(And(str, len), None),
    "ensembl_protein": Or(And(str, len), None),
    "tag": Or(And(str, len), None)
}

gene_feature_schema = {
    'Gene_Reference': Or([gene_reference_standard_schema], None),
    'Products': Or([product_feature_schema], None)
}

related_sequence_schema = Schema({
    'query': And(str, len),
    'organism': And(str, len),
    'moltype': And(str, len),
    "current_accession": And(str, len),
    'new_related': Or([{
        'Assemblies': assembly_accession_schema,
        "Genes": gene_feature_schema,
    }], None)
})












