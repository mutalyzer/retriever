"""
Module for NCBI Datsets response parsing.
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/
"""
from mutalyzer_retriever.reference import GRCH37
from mutalyzer_retriever.util import HUMAN_TAXON, DataSource


def parse_assemblies(dataset_report):
    """
    Args:
        report (dict): One gene report from the NCBI Datasets dataset api endpoint.
    Returns:
        list[dict]: A list of assemblies, where each entry may contain:
            - "name" (str): Assembly name (e.g., "GRCh38.p14").
            - "accession" (str): Assembly RefSeq accession.
    """
    if not isinstance(dataset_report, dict):
        return []

    assemblies = []

    gene = dataset_report.get("gene", {})
    taxon_name = gene.get("taxname", "")
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

    # Add GRCh37 assembly for human since NCBI Datasets omits this version
    if taxon_name.upper() == HUMAN_TAXON and sequence_name:
        grch37_acc = GRCH37.get(sequence_name)
        if grch37_acc:
            assemblies.append({
                    "name": "GRCh37.p13",
                    "accession": grch37_acc,
            })

    return assemblies


def parse_gene(dataset_report):
    """
    Args:
        dataset_report (dict): One gene report from the NCBI Datasets dataset API.

    Returns:
        dict: Parsed gene data, which may include:
            - "hgnc_id" (str): HGNC identifier.
            - "name" (str): Gene symbol.
            - "description" (str): Gene description.
            - "providers" (list[dict]): List of providers (Ensembl or/and NCBI) for the gene
              with their corresponding accessions and ID.
    """
    if not isinstance(dataset_report, dict):
        return {}
    gene = dataset_report.get("gene", {})
    if not gene:
        return {}

    gene_entry = {}

    # get gene information from ncbi provider
    ncbi_gene = {}
    ncbi_id = gene.get("gene_id")
    if ncbi_id:
        ncbi_gene = {"name": DataSource.NCBI, "id": ncbi_id}

        for ref in gene.get("reference_standards", []):
            ref_accession = ref.get("gene_range", {}).get("accession_version")
            if ref_accession:
                ncbi_gene["accession"] = ref_accession
                break

    # get gene information from ensembl provider
    ensembl_gene = {}
    ensembl_ids = gene.get("ensembl_gene_ids") or []
    if ensembl_ids:
        ensembl_gene = {
                "name": DataSource.ENSEMBL,
                "accession": ensembl_ids[0]
            }

    # get gene information
    hgnc_id_raw = gene.get("nomenclature_authority", {}).get("identifier")
    if hgnc_id_raw and "HGNC:" in hgnc_id_raw:
        gene_entry["hgnc_id"] = hgnc_id_raw.split(":")[1]

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


def parse_transcripts(product_report):
    """
    Parse transcript and protein information from an NCBI Datasets product report.
    Args:
        product_report (dict): One gene report from the NCBI Datasets API.

    Returns:
        list[dict]: A list of transcript entries, where each entry may contain:
            - "providers" (list[dict]): List of providers (Ensembl or/and NCBI) for the transcripts and proteins
              with their corresponding accessions and ID.
            - "tag" (str, optional): Transcript tag, e.g., "MANE Select".
    """
    if not isinstance(product_report, dict):
        return []

    product = product_report.get("product", {})
    gene_products = []

    for transcript in product.get("transcripts", []):
        # Build ncbi transcripts and protein units
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


        # Build ensembl transcripts and protein units
        ensembl_transcript = {}
        ensembl_transcript_acc = transcript.get("ensembl_transcript")
        ensembl_protein_acc = ncbi_protein.get("ensembl_protein")
        if ensembl_transcript_acc:
            ensembl_transcript = {
                    "name": DataSource.ENSEMBL,
                    "transcript_accession": ensembl_transcript_acc
                }
        if ensembl_protein_acc:
            ensembl_transcript["protein_accession"] = ensembl_protein_acc

        # Combine providers
        providers = [
            p for p in (ncbi_transcript, ensembl_transcript) if len(p) > 1
        ]
        if not providers:
            continue

        product_entry = {"providers": providers}

        transcript_tag = transcript.get("select_category")
        if transcript_tag:
            product_entry["tag"] = transcript_tag
        gene_products.append(product_entry)

    return gene_products


def parse_dataset_report(dataset_report):
    """
    Parse transcript and protein information from an NCBI Datasets dataset report.
    Args:
        dataset_report (dict): One gene report from the NCBI Datasets API.
    Returns:
        dict:
            - taxon_name (str): Species name.
            - assemblies (list): Parsed genomic assemblies.
            - genes (list): Parsed gene entries with their products.
    """
    if not isinstance(dataset_report, dict):
        return {}
    reports = dataset_report.get("reports", [])
    if not reports:
        return {}

    output = {}

    # Extract taxon name and assemblies from the first report.
    # The input acccession for related module is from a single species
    # on the same chromosome, so parsing only the first report is sufficient.
    first_report = reports[0]
    taxname = first_report.get("taxname")
    if taxname:
        output["taxname"] = taxname
    assemblies = parse_assemblies(first_report)
    if assemblies:
        output["assemblies"] = assemblies

    # Parse genes and their products
    genes = []
    for report in reports:
        gene = parse_gene(report)
        if gene:
            genes.append(gene)
    if genes:
        output["genes"] = genes

    return output


def parse_product_report(product_report):
    """
    Parse a gene product report from the NCBI Datasets product report API.
    Args:
        product_report (dict): One report for one gene's products from the NCBI Datasets product API endpoint.
    Returns:
        dict:
            - "taxname" (str): Species or taxon name.
            - "genes" (list[dict]): List of gene entries, each containing:
                - "name" (str): Gene symbol.
                - "transcripts" (list[dict]): Parsed transcript and protein mappings.
    """
    if not isinstance(product_report, dict):
        return {}

    reports = product_report.get("reports", [])
    if not reports:
        return {}

    output = {}

    # Extract taxon name from the first report
    first_report = reports[0]
    taxname = first_report.get("taxname")
    if taxname:
        output["taxname"] = taxname

    genes = []
    for report in reports:
        product = {}

        gene_symbol = report.get("product", {}).get("symbol")
        if gene_symbol:
            product["name"] = gene_symbol

        transcripts = parse_transcripts(report)
        if transcripts:
            product["transcripts"] = transcripts

        if product:
            genes.append(product)
    output["genes"] = genes

    return output


def merge_datasets(genomic_related, product_related):
    """
    Merges parsed genomic and product-related from datasets.
    """
    if not product_related or not genomic_related:
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
                else:
                    merged_gene = genomic_gene
                merged_genes.append(merged_gene)
            if merged_genes:
                related["genes"] = merged_genes

    # Sort genes alphabetically by name
    if related.get("genes"):
        related["genes"].sort(key=lambda g: g.get("name", ""))

    return related
