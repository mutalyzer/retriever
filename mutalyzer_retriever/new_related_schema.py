from schema import Schema, And, Or, Optional

assembly_accession_schema = {str: And(str, len)}

transcript_product_feature_schema = {
    'NCBI_accession': Or(And(str, len), None),
    'ENSEMBL_accession': Or(And(str, len), None)
}

genes_product_feature_schema = {
    'NCBI_accession': Or(And(str, len), None),
    'ENSEMBL_accession': Or(And(str, len), None),
    'Products': Or([transcript_product_feature_schema], None),
     Optional('Tag') : Or(And(str, len), None)
}

genes_feature_schema = {
    'Gene_name': Or(And(str, len), None),
    'RefSeqGene': Or(And(str, len), None),
    'Products': Or([genes_product_feature_schema], None)
}

related_sequence_schema = Schema({
    'Query': And(str, len),
    'Organism': And(str, len),
    'Moltype': And(str, len),
    'Current_accession': And(str, len),
    'New_related': Or({
        'Assemblies': assembly_accession_schema,
        'Genes': Or([genes_feature_schema],None),
    }, None)
})

