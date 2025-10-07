"""
Module for NCBI Datsets response parsing.
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/
"""
from mutalyzer_retriever.util import DataSource, HUMAN_TAXON
from mutalyzer_retriever.reference import GRCH37

def _parse_assemblies(report):
    taxon_name = None
    assemblies = []

    gene = report.get("gene", {})
    taxon_name = taxon_name or gene.get("taxname")
    sequence_name = None

    for annotation in gene.get("annotations", []):
        for loc in annotation.get("genomic_locations", []):
            accession = loc.get("genomic_accession_version")
            assembly_name = annotation.get("assembly_name")
            sequence_name = loc.get("sequence_name")
            assembly_entry = {}
            if assembly_name:
                assembly_entry["name"] = assembly_name
            if accession:
                assembly_entry["accession"] = accession

            if assembly_entry and assembly_entry not in assemblies:
                assemblies.append(assembly_entry)

    # Add GRCh37 assembly if human
    if taxon_name.upper() == HUMAN_TAXON and sequence_name:
        grch37_acc = GRCH37.get(sequence_name)
        if grch37_acc:
            grch37_entry = {
                    "name": "GRCh37.p13",
                    "accession": grch37_acc,
            }
            assemblies.append(grch37_entry)
    return assemblies


def _parse_genes(dataset_report):
    if not dataset_report.get("gene"):
        return {}
    
    gene_entry = {}

    gene = dataset_report.get("gene", {})
    # get gene information from ncbi provider
    ncbi_gene = {}
    ncbi_id = gene.get("gene_id")
    if ncbi_id:
        ncbi_gene = {"name": DataSource.NCBI, "id": ncbi_id}

    ref_standards = gene.get("reference_standards", [])
    for ref in ref_standards:
        ref_accession = ref.get("gene_range", {}).get("accession_version")
        if ref_accession:
            ncbi_gene["accession"] = ref_accession

    # get gene information from ensembl provider
    ensembl_gene = {}
    ensembl_ids = gene.get("ensembl_gene_ids")
    ensembl_id = ensembl_ids[0] if ensembl_ids else None
    if ensembl_id:
        ensembl_gene = {
                "name": DataSource.ENSEMBL,
                "accession": ensembl_ids[0] if ensembl_ids else None,
            }
        
    # get gene information
    hgnc_id_raw = gene.get("nomenclature_authority", {}).get("identifier")
    if hgnc_id_raw and ":" in hgnc_id_raw:
        hgnc_id = hgnc_id_raw.split(":")[1] if "HGNC" in hgnc_id_raw else None
        if hgnc_id:
            gene_entry["hgnc_id"] = hgnc_id 
    symbol = gene.get("symbol")
    if symbol:
        gene_entry["name"] = symbol
    description = gene.get("description")
    if description:
        gene_entry["description"] = description
    providers = [p for p in (ncbi_gene, ensembl_gene) if len(p) > 1]
    if providers:
        gene_entry["providers"] = providers

    return gene_entry


def _parse_transcripts(product_report):
    product = product_report.get("product", {})
    gene_products = []

    for transcript in product.get("transcripts", []):
        
        # build ncbi transcripts and protein units
        ncbi_transcript_acc = transcript.get("accession_version")
        ncbi_desc = transcript.get("name")
        ncbi_protein = transcript.get("protein", {})
        ncbi_protein_acc = ncbi_protein.get("accession_version")
        ncbi_transcript = {}
        if ncbi_transcript_acc:
            ncbi_transcript = {
                "name": DataSource.NCBI,
                "transcript_accession": ncbi_transcript_acc
            }
        if ncbi_desc:
            ncbi_transcript["description"] = ncbi_desc
        if ncbi_protein_acc:
            ncbi_transcript["protein_accession"] = ncbi_protein_acc

        
        # build ensembl transcripts and protein units
        ensembl_transcript = {}
        ensembl_acc = transcript.get("ensembl_transcript")
        ensembl_protein_acc = ncbi_protein.get("ensembl_protein")
        if ensembl_acc:
            ensembl_transcript = {
                    "name": DataSource.ENSEMBL,
                    "transcript_accession": ensembl_acc
                }
        if ensembl_protein_acc:
            ensembl_transcript["protein_accession"] = ensembl_protein_acc

        providers = [
            p for p in (ncbi_transcript, ensembl_transcript) if len(p) > 1
        ]
        if not providers:
            continue

        product_entry = {"providers": providers}
        if  transcript.get("select_category"):
            product_entry["tag"] = transcript.get("select_category")
        gene_products.append(product_entry)

    return gene_products


def empty_reports(reports):
    if reports.get("reports"):
        return True

def _parse_dataset_report(dataset_report):
    """
    Parses a gene dataset report from the NCBI Datasets API.

    Args:
        json_report (dict): JSON object returned from the NCBI Datasets API.

    Returns:
        dictionary:
            - taxon_name (str): .
            - assemblies (list):
            - genes (list):
    """
    if not empty_reports(dataset_report):
        return {}
    
    output = {}
    taxname = dataset_report.get("reports")[0].get("taxname")
    if taxname:
        output["taxname"] = taxname    

    genes = []
    for report in dataset_report.get("reports", []):
        # Extract genomic assemblies and genes
        assemblies = _parse_assemblies(report)
        if assemblies:
            output["assemblies"] = assemblies
        gene = _parse_genes(report)
        if gene:
            genes.append(gene)
    if genes:
        output["genes"] = genes

    return output


def _parse_product_report(product_report):
    """
    Parse NCBI product report JSON into a dictionary.
    Args:
        data (dict): JSON data from ncbi product API.
    Returns:
        dict: Mapping gene_symbol -> list of transcript info dicts
    """
    if not empty_reports(product_report):
        return {}

    output = {}
    taxname = product_report.get("reports")[0].get("taxname")
    if taxname:
        output["taxname"] = taxname
    
    genes = []
    for report in product_report.get("reports", []):
        product = {}
        transcripts = _parse_transcripts(report)
        if transcripts:
            product["transcripts"] = transcripts

        gene_symbol = report.get("product", {}).get("symbol")
        if gene_symbol:
            product["name"] = gene_symbol
        genes.append(product)
    output["genes"] = genes

    return output


def _merge_datasets(genomic_related, product_related):
    """
    Merges genomic and product-related from Datasets.
    """
    if not (product_related and genomic_related) or not genomic_related:
        return {}

    related = {}

    merged_genes = []
    if genomic_related.get("assemblies"):
        related["assemblies"] = genomic_related.get("assemblies")

    for genomic_gene in genomic_related.get("genes", []):
        symbol = genomic_gene.get("name")
        for product_gene in product_related.get("genes", []):
            if symbol == product_gene.get("name"):
                transcripts = product_gene.get("transcripts")
                if transcripts:
                    gene_products = {"transcripts": transcripts}
                    merged_gene =  genomic_gene | gene_products
                    merged_genes.append(merged_gene)
            if merged_genes:
                related["genes"] = merged_genes

    # Sort genes alphabetically by name
    related["genes"].sort(key=lambda g: g.get("name", ""))

    return related
