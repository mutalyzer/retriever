import json
import re
from copy import deepcopy

import requests

from mutalyzer_retriever.client import EnsemblClient, NCBIClient, NCBIEutilsClient
from mutalyzer_retriever.configuration import cache_url, settings
from mutalyzer_retriever.parsers import datasets, ensembl_gene_lookup
from mutalyzer_retriever.request import request
from mutalyzer_retriever.util import (
    DEFAULT_TIMEOUT,
    HUMAN_TAXON,
    DataSource,
    MoleculeType,
)


def get_cds_to_mrna(cds_id, timeout=DEFAULT_TIMEOUT):
    def _get_from_api_cache():
        api_url = cache_url()
        if api_url:
            url = api_url + "/cds_to_mrna/" + cds_id
            try:
                annotations = json.loads(requests.get(url, timeout=timeout).text)
            except Exception:
                return
            if annotations.get("mrna_id"):
                return annotations["mrna_id"]

    mrna_id = _get_from_api_cache()
    if mrna_id:
        return mrna_id

    ncbi = _fetch_ncbi_datasets_gene_accession(cds_id, timeout)
    if (
        ncbi.get("genes")
        and len(ncbi["genes"]) == 1
        and ncbi["genes"][0].get("gene")
        and ncbi["genes"][0]["gene"].get("transcripts")
    ):
        transcripts = ncbi["genes"][0]["gene"]["transcripts"]
        mrna_ids = set()
        for transcript in transcripts:
            if (
                transcript.get("accession_version")
                and transcript.get("protein")
                and transcript["protein"].get("accession_version") == cds_id
            ):
                mrna_ids.add(transcript["accession_version"])
        return sorted(list(mrna_ids))


def _fetch_ncbi_datasets_gene_accession(accession_id, timeout=TimeoutError):
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/gene/accession/{accession_id}/product_report"
    return json.loads(request(url=url, timeout=timeout))


def filter_report_from_other_genes(gene_symbol: str, reports: dict):
    "Return a report matching the given gene symbol from the NCBI reports data."
    # NCBI datasets would return related genes when query by gene name
    # e.g, query with CYP2D6 would get info of CYP2D7
    # https://api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol/CYP2D6D%2CCYP2D6/taxon/9606/product_report
    for report in reports.get("reports", {}):
        for key, value in report.items():
            if isinstance(value, dict) and value.get("symbol") in [gene_symbol, gene_symbol.upper()]:
                return {"reports": [{key: value}]}
    return {}


def _convert_locations(accession: str, locations):
    """
    Check if input locations are in the format of 'start_end' or single points,
    and return a normalized string like: accession:start_end;accession:start_end
    """
    if locations is None:
        raise ValueError(f"Unkown location on {accession}.")

    pattern_range = re.compile(r"^\d+_\d+$")
    pattern_point = re.compile(r"^\d+$")

    valid_locations = []

    for item in locations.split(";"):
        item = item.strip()
        if pattern_range.match(item):
            raw_start, raw_end = map(int, item.split("_"))
            start = min(raw_start, raw_end)
            end = max(raw_start, raw_end)
        elif pattern_point.match(item):
            start = int(item)
            end = start + 1
        else:
            raise NameError(
                f"Invalid location format: '{item}'. Expected format: point or range 'point;start_end'."
            )

        valid_locations.append(f"{accession}:{start}-{end}")

    return valid_locations


def _merge_provider(ensembl_provider, ncbi_provider):
    """
    Merge Ensembl provider information into a list of NCBI providers.

    Args:
        ensembl_provider (dict): The provider entry from Ensembl to merge.
        ncbi_provider (list): A list of provider dictionaries from NCBI.

    Returns:
        list or dict: A new list of providers with the Ensembl entry merged in,
                      or the Ensembl provider dict if NCBI providers were not present.
    """

    if not ncbi_provider:
        return ensembl_provider
    providers = deepcopy(ncbi_provider)
    if len(providers) > 1:
        for i, p in enumerate(ncbi_provider):
            if p.get("name") == DataSource.ENSEMBL:
                providers[i] = ensembl_provider
    else:
        providers.append(ensembl_provider)
    return providers


def ncbi_match(ncbi_data, ensembl_id):
    """
    Check if an NCBI data entry contains a provider from Ensembl
    that matches the given Ensembl ID.

    Args:
        ncbi_data (dict): A gene or transcript entry from NCBI data.
        ensembl_id (str): An Ensembl accession to match.

    Returns:
        list: The list of providers if a match is found, otherwise an empty list.
    """
    providers = ncbi_data.get("providers", [])
    for p in providers:
        if (
            p.get("name") == DataSource.ENSEMBL and
            (
                p.get("accession") == ensembl_id or
                p.get("transcript_accession") == ensembl_id
            )
        ):
            return providers
    return []



def _merge_assemblies(ensembl_related, ncbi_related):
    "Merge two lists of assemblies gathered from ensembl and ncbi"
    return (
        ensembl_related.get("assemblies", []) +
        ncbi_related.get("assemblies", [])
    )


def _merge_transcripts(ensembl_related, ncbi_gene):
    """
    Merge transcript data from Ensembl and NCBI sources for a given gene.

    Args:
        ensembl_related (dict): Ensembl related for a single gene.
        ncbi_gene (dict): NCBI related for the same gene.

    Returns:
        list: A merged list of transcript entries with combined provider data.
    """
    ensembl_transcripts = ensembl_related.get("transcripts") or []
    ncbi_transcripts = ncbi_gene.get("transcripts") or []

    if not ensembl_transcripts:
        return ncbi_transcripts

    merged = []

    for ensembl_t in ensembl_transcripts:
        ensembl_accession = ensembl_t.get("transcript_accession")
        matched = False
        # shape ensembl gene data and merge gene info from two sources
        ensembl_entry = {"name": DataSource.ENSEMBL, **ensembl_t}

        for ncbi_t in ncbi_transcripts:
            transcript = deepcopy(ncbi_t)
            if (transcript
                and len(transcript.get("providers",[])) > 1
                and ncbi_match(ncbi_t, ensembl_accession)
            ):
                transcript["providers"] = _merge_provider(
                    ensembl_entry, ncbi_t.get("providers")
                )
                merged.append(transcript)
                matched = True
                break
            if (
                transcript
                and transcript.get("providers")
                and len(transcript["providers"]) == 1
                and ncbi_t not in merged
            ):
                merged.append(ncbi_t)

        if not matched:
            # No match found in NCBI, add Ensembl-only entry
            merged.append({"providers": [ensembl_entry]})
    return merged


def _merge_gene(ensembl_related, ncbi_related):
    "Merge two lists of related genes gathered from ensembl and ncbi"
    ensembl_gene_name = ensembl_related.get("name")
    ensembl_gene_accession = ensembl_related.get("accession")
    if not (ensembl_gene_name and ensembl_gene_accession):
        return ncbi_related.get("genes", [])
    for ncbi_gene in ncbi_related.get("genes", []):
        if ncbi_gene.get("name") == ensembl_gene_name:
            # shape ensembl gene data and merge gene info from two sources
            ensembl_entry = {
                "name": DataSource.ENSEMBL,
                "accession": ensembl_gene_accession,
            }
            gene = {}
            for key, value in ncbi_gene.items():
                if key == "providers":
                    gene["providers"] = _merge_provider(
                        ensembl_entry, ncbi_gene["providers"]
                    )
                elif key == "transcripts":
                    gene["transcripts"] = _merge_transcripts(
                        ensembl_related, ncbi_gene
                    )
                else:
                    gene[key] = value
            return [gene]
    return []


def _merge(ensembl_related, ncbi_related):
    """Merge related gathered from two sources"""
    merged = {}
    merged_assemblies = _merge_assemblies(ensembl_related, ncbi_related)
    if merged_assemblies:
        merged["assemblies"] = merged_assemblies
    merged_genes = _merge_gene(ensembl_related, ncbi_related)
    if merged_genes:
        merged["genes"] = merged_genes
    return merged


def _get_related_by_gene_symbol_from_ncbi(gene_symbol, taxon_name=HUMAN_TAXON):
    """
    Given a gene symbol, return a set of related sequence accessions (genomic and/or products).
    Returns related_dict, or {} if nothing found.
    """
    if not gene_symbol:
        return {}
    related = {}

    client = NCBIClient(timeout=DEFAULT_TIMEOUT)
    dataset_response = client.get_gene_symbol_dataset_report(
        gene_symbol, taxon_name
    )
    filtered_dataset_response = filter_report_from_other_genes(gene_symbol, dataset_response)
    parsed_dataset = datasets.parse_dataset_report(filtered_dataset_response)
    product_response = client.get_gene_symbol_product_report(
        gene_symbol, taxon_name
    )
    filtered_product_response = filter_report_from_other_genes(gene_symbol, product_response)
    parsed_product = datasets.parse_product_report(filtered_product_response)
    related = datasets.merge_datasets(
        parsed_dataset, parsed_product
    )
    return related


def _get_related_by_gene_symbol_from_ensembl(gene_symbol):
    """Fetch related gene data from Ensembl using a gene symbol"""
    if not gene_symbol:
        return {}
    client = EnsemblClient(timeout=DEFAULT_TIMEOUT)
    gene_lookup_response = client.lookup_symbol(gene_symbol)
    return ensembl_gene_lookup.parse_ensembl_gene_lookup_json(
        gene_lookup_response
    )


def _get_related_by_gene_symbol(gene_symbol):
    """
    Fetch related genomic and product sequences from NCBI endpoint using a gene symbol.
    This function contacts the NCBI Datasets and Ensembl to gathered related from both sources.
    Merge and filter the related from two sources.

    Args:
        accession (str): A RefSeq accession.
    Returns:
        related (dict or None): A dictionary of related sequences; otherwise {}.
    Raises:
        RuntimeError: If the NCBI Datasets API is unavailable or returns an invalid response.
    """
    ncbi_related = _get_related_by_gene_symbol_from_ncbi(gene_symbol, taxon_name=HUMAN_TAXON)
    ensembl_related = _get_related_by_gene_symbol_from_ensembl(gene_symbol)
    related = _merge(ensembl_related, ncbi_related)
    related = filter_related(gene_symbol, related)
    return related


def _parse_ensembl_transcript_lookup_json(ensembl_transcript_json):
    """Parse an Ensembl transcript JSON object into related gene data."""
    ensembl_gene_id = ensembl_transcript_json.get("Parent")
    client = EnsemblClient(timeout=DEFAULT_TIMEOUT)
    ensembl_gene_json = client.lookup_id(ensembl_gene_id, expand=1)
    return _parse_ensembl(ensembl_gene_json)


def _parse_ensembl_protein_lookup_json(ensembl_protein_json):
    """Parse an Ensembl protein JSON object into related gene data."""
    ensembl_transcript_id = ensembl_protein_json.get("Parent")
    client = EnsemblClient(timeout=DEFAULT_TIMEOUT)
    ensembl_transcript_json = client.lookup_id(
        ensembl_transcript_id, expand=1
    )
    return _parse_ensembl(ensembl_transcript_json)


def _parse_ensembl(ensembl_json):
    """Dispatch Ensembl JSON parsing based on object type (Gene, Transcript, or Translation)."""
    if not ensembl_json:
        return {}
    obj_type = ensembl_json.get("object_type")
    if obj_type == "Gene" and ensembl_json.get("Transcript"):
        return ensembl_gene_lookup.parse_ensembl_gene_lookup_json(
            ensembl_json
        )
    if obj_type == "Transcript" and ensembl_json.get("Parent"):
        return _parse_ensembl_transcript_lookup_json(ensembl_json)
    if obj_type == "Translation" and ensembl_json.get("Parent"):
        return _parse_ensembl_protein_lookup_json(ensembl_json)

    raise ValueError(f"Unsupported Ensembl object: {obj_type}")


def _get_ensembl_gene_info(accession_base, moltype):
    """
    Get gene-related information from Ensembl for a given accession.

    Args:
        accession_base (str): Ensembl accession (gene, transcript, or protein).
        moltype (str): Molecule type, determines the level of detail to fetch.
                       Should be one of MoleculeType.DNA, RNA, or PROTEIN.
    Returns:
        dict: Parsed gene-related information from Ensembl, including transcripts
              and/or proteins depending on the accession type.
    Raises:
        ValueError: If the molecule type is unsupported or the Ensembl lookup fails.
    """

    if moltype not in {
        MoleculeType.DNA,
        MoleculeType.RNA,
        MoleculeType.PROTEIN,
    }:
        raise ValueError(f"Unsupported molecule type: {moltype}")

    client = EnsemblClient(timeout=DEFAULT_TIMEOUT)
    expand_flag = 0 if moltype == MoleculeType.PROTEIN else 1

    ensembl_lookup_json = client.lookup_id(accession_base, expand=expand_flag)

    return _parse_ensembl(ensembl_lookup_json)


def filter_assemblies(related_assemblies):
    """
    Filter out non-chromosomal genomic accessions from a list of assemblies.
    Only retains assemblies with accessions that start with 'NC_'.
    Args:
        related_assemblies (list[dict]): A list of assembly dictionaries, each containing key "accession".
    Returns:
        list[dict]: A filtered list containing only assemblies with 'NC_' accessions.
    """
    return [
        assembly
        for assembly in related_assemblies
        if assembly.get("accession", "").startswith("NC_")
    ]


def filter_genes(accession_base, related_genes):
    """
    Filter a list of gene entries by updating their transcript lists.

    Args:
        accession_base (str): The base accession to match (no version).
        related_genes (list): List of gene dictionaries with products info.

    Returns:
        list: A list of genes, each optionally containing filtered transcripts.
    """
    if not accession_base:
        return []

    filtered_genes = []
    for gene in related_genes:
        # Make a copy of this gene, only make changes in the 'transcripts'.
        filtered_gene = {k: v for k, v in gene.items() if k != "transcripts"}
        filtered_transcripts = filter_transcripts(accession_base, gene)
        if filtered_transcripts:
            filtered_gene["transcripts"] = filtered_transcripts

        filtered_genes.append(filtered_gene)

    return filtered_genes


def filter_transcripts(accession_base, gene):
    """
    This function filters products of one gene and keeps only transcripts/proteins that either:
        - Match the base accession (ignoring version), or
        - Have a MANE Select or other tag.

    Args:
        accession_base (str): The base accession to match (no version).
        gene (dict): A dictionary contains a gene's products and proteins from different sources.

    Returns:
        list: A list of filtered products
    """
    filtered_transcripts = []

    for transcript in gene.get("transcripts", []):
        if "tag" in transcript:
            filtered_transcripts.append(transcript)
        else:
            for provider in transcript.get("providers", []):
                products = [
                    (provider.get("transcript_accession") or "").split(".")[0],
                    (provider.get("protein_accession") or "").split(".")[0],
                ]
                if accession_base in products:
                    filtered_transcripts.append(transcript)
                    break

    return filtered_transcripts


def filter_related(accession_base, related):
    """
    Filter related data for a given accession by removing non-chromosomal assemblies
    and less primary gene products.

    Args:
        accession_base (str): The base accession to match against related gene products.
        related (dict): A dictionary that may contain:
            - "assemblies" (list): List of assembly records.
            - "genes" (list): List of gene records, each with possible transcripts/proteins.

    Returns:
        dict: A filtered dictionary containing:
            - "assemblies": Only chromosomal assemblies.
            - "genes": Genes with at least MANE Select transcripts or products.
    """


    filtered = {}
    filtered_assemblies = filter_assemblies(related.get("assemblies", []))
    if filtered_assemblies:
        filtered["assemblies"] = filtered_assemblies
    filtered_gene_products = filter_genes(
        accession_base, related.get("genes", [])
    )
    if filtered_gene_products:
        filtered["genes"] = filtered_gene_products
    return filtered


def _get_gene_related(gene_ids):
    client = NCBIClient(timeout=DEFAULT_TIMEOUT)
    product_response = client.get_gene_id_product_report(gene_ids)
    parsed_product = datasets.parse_product_report(product_response)
    dataset_response = client.get_gene_id_dataset_report(gene_ids)
    parsed_dataset = datasets.parse_dataset_report(dataset_response)
    taxname = product_response.get("taxname")
    return taxname, datasets.merge_datasets(parsed_dataset, parsed_product)


def _parse_genome_annotation_report(genome_annotation_report):
    """
    Extract gene IDs and taxon name from a genome annotation report,
    and retrieve related gene and assembly data.

    Args:
        genome_annotation_report (dict): The genome annotation response from NCBI,
            expected to contain a list of reports with gene annotations.

    Returns:
        list:
            - gene_ids: NCBI gene Ids.
    """

    gene_ids = []
    for report in genome_annotation_report.get("reports", []):
        annotation = report.get("annotation", {})
        gene_id = annotation.get("gene_id")
        if gene_id is not None:
            gene_ids.append(gene_id)
    return gene_ids


def _get_related_by_chr_location(accession, locations):
    """
    Retrieve related gene and assembly information for a chromosomal NCBI accession.

    This function uses the accession and genomic locations to fetch the corresponding
    genome annotation report from NCBI datasets, and parses it to extract the taxon name and
    related gene/assembly data.

    Args:
        accession (str): A chromosomal NCBI accession.
        locations (str): Genomic locations associated with the accession.

    Returns:
        tuple: A tuple containing:
            - taxon_name (str): The organism name.
            - related (dict): Related gene and assembly data parsed from the annotation.

    Raises:
        ValueError: If the assembly accession cannot be determined for the given input.
    """

    assembly_accession = _get_assembly_accession(accession)
    if not assembly_accession:
        raise ValueError(
            f"Assembly accession could not be determined for {accession}"
        )

    client = NCBIClient(timeout=DEFAULT_TIMEOUT)
    annotation_response = client.get_genome_annotation_report(
        assembly_accession, _convert_locations(accession, locations)
    )
    gene_ids = _parse_genome_annotation_report(annotation_response)
    if gene_ids:
        taxon_name, related = _get_gene_related(gene_ids)
        return taxon_name, related
    return None, {}


def _gene_ids_from_ng(accession):
    """
    Return a list of NCBI Gene IDs associated with an NG_ (RefSeqGene) accession.
    """
    client = NCBIEutilsClient(timeout=DEFAULT_TIMEOUT)
    data = client.elink(accession, dbfrom="nuccore", db="gene")

    gene_ids = []
    for linkset in data.get("linksets", []):
        for linkdb in linkset.get("linksetdbs", []):
            if linkdb.get("dbto") == "gene":
                gene_ids.extend(linkdb.get("links", []))

    return list(dict.fromkeys(gene_ids))


def _get_related_by_ng_accession(accession):
    accession_base = accession.split(".")[0]
    gene_ids = _gene_ids_from_ng(accession_base)
    if gene_ids:
        return _get_gene_related(gene_ids)
    return None, {}


def _get_assembly_accession(accession):
    """
    Retrieve the assembly accession corresponding to a given NCBI chromosome accession.
    """
    client = NCBIClient(timeout=DEFAULT_TIMEOUT)
    return client.get_assembly_accession(accession)


def _get_related_by_accession_from_ncbi(accession):
    """
    Fetch related genomic and product sequences from NCBI endpoint using a RefSeq accession.
    This function contacts the NCBI Datasets API endpoints to gather both dataset report and product report
    for a given transcript or protein accession.

    Args:
        accession (str): A RefSeq accession.
    Returns:
        tuple:
            - taxon_name (str or None): The organism name associated with input accession.
            - related (dict or None): A dictionary of related sequences; otherwise {}.
    Raises:
        RuntimeError: If the NCBI Datasets API is unavailable or returns an invalid response.
    """

    related_from_ncbi = {}

    client = NCBIClient(timeout=DEFAULT_TIMEOUT)
    product_report = client.get_accession_product_report(accession)
    dataset_report = client.get_accession_dataset_report(accession)

    parsed_products = datasets.parse_product_report(product_report)
    parsed_dataset = datasets.parse_dataset_report(dataset_report)

    related_from_ncbi = datasets.merge_datasets(parsed_dataset, parsed_products)
    taxname = product_report.get("taxname")
    return taxname, related_from_ncbi


def _get_related_by_ensembl_id(accession, moltype):
    """
    Retrieve and filter related gene and assembly data for a given Ensembl accession
    by querying both Ensembl and NCBI sources.

    Args:
        accession (str): An Ensembl accession.
        moltype (str): The molecule type, one of MoleculeType.DNA, RNA, or PROTEIN.

    Returns:
        dict: A dictionary containing related data, with optional keys:
            - "assemblies": List of chromosomal assemblies (RefSeq NC_ accessions).
            - "genes": List of related genes with filtered transcripts or proteins.
    """

    accession_base = accession.split(".")[0]

    # Get taxon_name and related from ensembl
    ensembl_related = _get_ensembl_gene_info(accession_base, moltype)

    ncbi_related = {}
    # Get related from ncbi using gene symbol and taxname
    gene_symbol = ensembl_related.get("name")
    taxname = ensembl_related.get("taxon_name")
    if gene_symbol and taxname:
        ncbi_related = _get_related_by_gene_symbol_from_ncbi(gene_symbol, taxname)

    related = _merge(ensembl_related, ncbi_related)
    if taxname and taxname.upper() == HUMAN_TAXON:
        return filter_related(accession_base, related)

    return related


def _get_related_by_ncbi_id(accession, moltype, locations):
    """
    Retrieve and filter related gene and assembly data for a given NCBI accession
    by querying NCBI and Ensembl sources.

    Args:
        accession (str): An NCBI accession.
        moltype (str): The molecule type, one of MoleculeType.DNA, RNA, or PROTEIN.
        locations (str): Genomic location(s) if applicable (used for DNA/chromosomal).

    Returns:
        dict: A dictionary containing filtered related data, with optional keys:
            - "assemblies": List of chromosomal assemblies (RefSeq NC_ accessions).
            - "genes": List of related genes with filtered transcripts or proteins.

    Raises:
        NameError: If the accession cannot be retrieved or resolved via NCBI.
    """
    accession_base = accession.split(".")[0]

    # Get taxon_name and related from ncbi datasets
    if moltype in [MoleculeType.RNA, MoleculeType.PROTEIN]:
        taxon_name, ncbi_related = _get_related_by_accession_from_ncbi(accession)
    elif moltype == MoleculeType.DNA and "NC_" in accession:
        taxon_name, ncbi_related = _get_related_by_chr_location(accession, locations)
    elif moltype == MoleculeType.DNA and "NG_" in accession:
        taxon_name, ncbi_related = _get_related_by_ng_accession(accession)
    else:
        raise NameError(f"Could not retrieve {accession} from NCBI.")

    # Get related from ensembl using ensembl gene id.
    related = {}
    if ncbi_related:
        # All ensembl gene IDs
        ensembl_genes_id = [
            provider["accession"]
            for gene in ncbi_related.get("genes", [])
            for provider in gene.get("providers", [])
            if provider.get("name") == DataSource.ENSEMBL
        ]

        if not ensembl_genes_id:
            related = ncbi_related
            if taxon_name and taxon_name.upper() == HUMAN_TAXON:
                related = filter_related(accession_base, related)
            return related

        all_genes = []
        all_assemblies = []
        for ensembl_gene in ensembl_genes_id:
            ensembl_gene_related = _get_ensembl_gene_info(ensembl_gene, moltype=MoleculeType.DNA)
            if not taxon_name:
                taxon_name = ensembl_gene_related.get("taxon_name")
            if ensembl_gene_related:
                related_entry = _merge(ensembl_gene_related, ncbi_related)
                if not all_assemblies and related_entry.get("assemblies"):
                    all_assemblies = related_entry.get("assemblies")
                    related = {"assemblies": all_assemblies}
                if related_entry.get("genes"):
                    all_genes.extend(related_entry.get("genes"))
        if all_genes:
            related["genes"] = all_genes

        if taxon_name and taxon_name.upper() == HUMAN_TAXON:
            related = filter_related(accession_base, related)

    return related


def parse_ensembl_id(seq_id):
    """
    Parse a sequence id to determine its source and molecule type.
    Args:
        seq_id (str): The sequence ID string to parse.
    Returns:
        tuple:
            - DataSource.ENSEMBL (enum): if the sequence id is in the format of
            a recognised ENSEMBL accession
            - moltype (str): One of 'dna', 'rna', 'protein' and others.
        Returns (DataSource.ENSEMBL, MoleculeType.UNKNOWN) if format is invalid or not supported.
        Returns (None, None) is it is not from Ensembl.
    """
    # Ensembl declares its identifiers should be in the form of
    # ENS[species prefix][feature type prefix][a unique eleven digit number]
    # See at https://www.ensembl.org/info/genome/stable_ids/index.html

    if not seq_id.startswith("ENS"):
        return None, None
    ensembl_feature_map = {
        "E": "exon",
        "FM": "protein family",
        "G": "dna",
        "GT": "gene tree",
        "P": "protein",
        "R": "regulatory feature",
        "T": "rna",
    }
    ensembl_pattern = re.compile(
        r"^ENS[A-Z]*?(FM|GT|G|T|P|R|E)\d{11}(?:\.\d+)?$"
    )
    match = ensembl_pattern.match(seq_id)
    if match:
        prefix = match.group(1)
        moltype = ensembl_feature_map.get(prefix, MoleculeType.UNKNOWN)
        return DataSource.ENSEMBL, moltype
    return DataSource.ENSEMBL, MoleculeType.UNKNOWN


def parse_ncbi_id(seq_id):
    """
    Parse a sequence id to determine its source and molecule type.
    Args:
        seq_id (str): The sequence ID string to parse.
    Returns:
        tuple:
            - DataSource.NCBI (enum): if the sequence id is in the format of
            a recognised NCBI RefSeq accession
            - moltype (str): One of 'dna', 'rna', 'protein' or unknown.
        Returns (None, None) if the sequence ID does not match RefSeq format.
    """

    # NCBI RefSeq prefix history:
    # https://www.ncbi.nlm.nih.gov/books/NBK21091/table/
    # ch18.T.refseq_accession_numbers_and_mole/?report=objectonly
    refseq_moltype_map = {
        # DNA
        "AC": "dna",
        "NC": "dna",
        "NG": "dna",
        "NT": "dna",
        "NW": "dna",
        "NZ": "dna",
        # RNA
        "NM": "rna",
        "XM": "rna",
        "NR": "rna",
        "XR": "rna",
        # Protein
        "NP": "protein",
        "YP": "protein",
        "XP": "protein",
        "WP": "protein",
        "AP": "protein",
    }
    refseq_pattern = re.compile(r"^([A-Z]{2})_\d+(?:\.\d+)?$")
    match = refseq_pattern.match(seq_id)
    if match:
        prefix = match.group(1)
        moltype = refseq_moltype_map.get(prefix, "unknown")
        return DataSource.NCBI, moltype
    return None, None


def detect_sequence_source(seq_id):
    """
    Detects the source and molecular type of a sequence ID.
    Args:
        seq_id (str): The sequence ID string to evaluate.
    Returns: (source: str, moltype: str):
                source: DataSource.ENSEMBL, DataSource.NCBI, or DataSource.OTHER.
                moltype: MoleculeType.DNA, MoleculeType.RNA, MoleculeType.PROTEIN, MoleculeType.UNKNOWN.
    """
    source, moltype = parse_ensembl_id(seq_id)
    if source:
        return source, moltype

    source, moltype = parse_ncbi_id(seq_id)
    if source:
        return source, moltype

    return DataSource.OTHER, MoleculeType.UNKNOWN


def get_related(accession, locations=None):
    """
    Retrieve related assembly/gene/transcript/protein information based
    on accession or gene symbol
    Args:
        accession (str): A sequence accession (e.g., RefSeq, Ensembl ID) or a human gene symbol.
        locations, optional (str): A point or a range on chromosome, in the format of
            '10000;120000_130000'. Defaults to None.

    Returns:
        related (dict): A dictionary containing related information retrieved from Ensembl, NCBI,

    Raises:
        NameError: If the given accession is not from NCBI RefSeq or ENSEMBL.
    """
    accession = accession.upper().strip()

    source, moltype = detect_sequence_source(accession)
    if source == DataSource.ENSEMBL and moltype != MoleculeType.UNKNOWN:
        return _get_related_by_ensembl_id(accession, moltype)
    if source == DataSource.NCBI and moltype != MoleculeType.UNKNOWN:
        return _get_related_by_ncbi_id(accession, moltype, locations)
    if source == DataSource.OTHER:
        return _get_related_by_gene_symbol(accession)
    raise NameError(f"Could not retrieve related for {accession}.")
