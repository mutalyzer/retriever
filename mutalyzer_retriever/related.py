import json
import re
import requests
from copy import deepcopy
from mutalyzer_retriever.util import (
    DataSource,
    MoleculeType,
    HUMAN_TAXON,
    DEFAULT_TIMEOUT,
)
from mutalyzer_retriever.configuration import cache_url, settings
from mutalyzer_retriever.request import Http400, request
from mutalyzer_retriever.reference import GRCH37
from mutalyzer_retriever.client import NCBIClient, EnsemblClient
from mutalyzer_retriever.parsers import datasets, ensembl_gene_lookup


def get_cds_to_mrna(cds_id, timeout=DEFAULT_TIMEOUT):
    def _get_from_api_cache():
        api_url = cache_url()
        if api_url:
            url = api_url + "/cds_to_mrna/" + cds_id
            try:
                annotations = json.loads(requests.get(url).text)
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


def _validate_locations(accession: str, locations):
    """
    Check if input locations are in the format of 'start_end' or single points,
    and return a normalized string like: accession:start_end;accession:start_end
    """
    pattern_range = re.compile(r"^\d+_\d+$")
    pattern_point = re.compile(r"^\d+$")

    valid_locations = []

    if locations == "0":
        raise ValueError(f"Unkown location on chromosomes.")

    for item in locations.split(";"):
        item = item.strip()
        if pattern_range.match(item):
            try:
                raw_start, raw_end = map(int, item.split("_"))
                start, end = sorted([raw_start, raw_end])
            except (ValueError, TypeError):
                raise NameError(
                    f"Invalid range format: '{item}'. Expected numeric values like '100_200'."
                )
        elif pattern_point.match(item):
            start = int(item)
            end = start + 1
        else:
            raise NameError(
                f"Invalid location format: '{item}'. Expected format: porint or range 'start_end'."
            )

        valid_locations.append(f"{accession}:{start}-{end}")

    return valid_locations


def _merge_provider(ensembl_provider, ncbi_provider):
    if not ncbi_provider:
        return ensembl_provider
    providers = deepcopy(ncbi_provider)
    for i, p in enumerate(ncbi_provider):
        if p.get("name") == DataSource.ENSEMBL:
            providers[i] = ensembl_provider
    return providers


def ncbi_match(ncbi_data, ensembl_id):
    providers = ncbi_data.get("providers", [])
    for p in providers:
        if p.get("name") == DataSource.ENSEMBL and (
            p.get("accession") or p.get("transcript_accession") == ensembl_id
        ):
            return providers


def _merge_assemblies(ensembl_related, ncbi_related):
    "Merge two lists of assemblies gathered from ensembl and ncbi"
    return ensembl_related.get("assemblies", []) + ncbi_related.get(
        "assemblies", []
    )


def _merge_transcripts(ensembl_related, ncbi_gene):
    "Merge two lists of transcripts gathered from ensembl and ncbi"
    ensembl_transcripts = ensembl_related.get("transcripts", [])
    ncbi_transcripts = ncbi_gene.get("transcripts", [])

    if ensembl_transcripts is None:
        return ncbi_transcripts

    merged = []

    ensembl_transcripts = ensembl_related.get("transcripts", [])
    for ensembl_t in ensembl_transcripts:
        ensembl_accession = ensembl_t.get("transcript_accession")
        matched = False
        # shape ensembl gene data and merge gene info from two sources
        ensembl_entry = {**ensembl_t, "name": DataSource.ENSEMBL}

        for ncbi_t in ncbi_transcripts:
            transcript = deepcopy(ncbi_t)
            if len(transcript.get("providers")) > 1 and ncbi_match(
                ncbi_t, ensembl_accession
            ):
                transcript["providers"] = _merge_provider(
                    ensembl_entry, ncbi_t.get("providers")
                )
                merged.append(transcript)
                matched = True
                break
            if (
                transcript.get("providers")
                and len(transcript["providers"]) == 1
                and {"providers": transcript.get("providers")} not in merged
            ):
                merged.append({"providers": transcript.get("providers")})

        if not matched:
            merged.append({"providers": [ensembl_entry]})
    return merged


def _merge_gene(ensembl_related, ncbi_related):
    # Add checks for empty
    "Merge two lists of related genes gathered from ensembl and ncbi"
    ensembl_gene_name = ensembl_related.get("name")
    ensembl_gene_accession = ensembl_related.get("accession")
    if not (ensembl_gene_name and ensembl_gene_accession):
        return ncbi_related.get("genes", [])

    for ncbi_gene in ncbi_related.get("genes", []):
        if ncbi_gene.get("name") == ensembl_gene_name:
            # shape ensembl gene data and merge gene info from two sources
            ensembl_entry = {
                "accession": ensembl_gene_accession,
                "name": DataSource.ENSEMBL,
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
    merged = {}
    if _merge_assemblies(ensembl_related, ncbi_related):
        merged["assemblies"] = _merge_assemblies(
            ensembl_related, ncbi_related
        )
    if _merge_gene(ensembl_related, ncbi_related):
        merged["genes"] = _merge_gene(ensembl_related, ncbi_related)
    return merged


def _get_related_by_gene_symbol_from_ncbi(
    gene_symbol, taxon_name=HUMAN_TAXON
):
    """
    Given a gene symbol, return a set of related sequence accessions (genomic and/or products).
    Returns related_dict, or (None, None) if nothing found.
    """
    if not gene_symbol:
        return None, {}
    related = {}

    client = NCBIClient(timeout=DEFAULT_TIMEOUT)
    dataset_response = client.get_gene_symbol_dataset_report(
        gene_symbol, taxon_name
    )
    parsed_dataset = datasets._parse_dataset_report(dataset_response)
    product_response = client.get_gene_symbol_product_report(
        gene_symbol, taxon_name
    )
    parsed_product = datasets._parse_product_report(product_response)
    related = datasets._merge_datasets(
        parsed_dataset, parsed_product
    )
    return related


def _get_related_by_gene_symbol_from_ensembl(gene_symbol):
    if not gene_symbol:
        return {}
    client = EnsemblClient(timeout=DEFAULT_TIMEOUT)
    gene_lookup_response = client.lookup_symbol(gene_symbol)
    return ensembl_gene_lookup._parse_ensembl_gene_lookup_json(
        gene_lookup_response
    )


def _get_related_by_gene_symbol(gene_symbol):
    ncbi_related = _get_related_by_gene_symbol_from_ncbi(
        gene_symbol, taxon_name=HUMAN_TAXON
    )
    ensembl_related = _get_related_by_gene_symbol_from_ensembl(gene_symbol)
    related = _merge(ensembl_related, ncbi_related)
    related = filter_related(gene_symbol, related)
    return related


def _parse_ensembl_transcript_lookup_json(ensembl_transcript_json):
    ensembl_gene_id = ensembl_transcript_json.get("Parent")
    client = EnsemblClient(timeout=DEFAULT_TIMEOUT)
    ensembl_gene_json = client.lookup_id(ensembl_gene_id, expand=1)
    return _parse_ensembl(ensembl_gene_json)


def _parse_ensembl_protein_lookup_json(ensembl_protein_json):
    ensembl_transcript_id = ensembl_protein_json.get("Parent")
    client = EnsemblClient(timeout=DEFAULT_TIMEOUT)
    ensembl_transcript_json = client.lookup_id(
        ensembl_transcript_id, expand=1
    )
    return _parse_ensembl(ensembl_transcript_json)


def _parse_ensembl(ensembl_json):
    """
    Dispatch parsing based on Ensembl object type.
    """
    obj_type = ensembl_json.get("object_type")
    if not obj_type:
        return {}
    if obj_type == "Gene" and ensembl_json.get("Transcript"):
        return ensembl_gene_lookup._parse_ensembl_gene_lookup_json(
            ensembl_json
        )
    if obj_type == "Transcript" and ensembl_json.get("Parent"):
        return _parse_ensembl_transcript_lookup_json(ensembl_json)
    if obj_type == "Translation" and ensembl_json.get("Parent"):
        return _parse_ensembl_protein_lookup_json(ensembl_json)

    raise ValueError(f"Unsupported Ensembl object: {obj_type}")


def fetch_ensembl_gene_info(accession_base, moltype):
    """
    Fetch gene-related information for a given Ensembl accession.

    Args:
        accession_base (str): Ensembl accession, support for gene, transcript, protein.
        moltype (str): Molecule type — features parsed from ensembl accession prefix.

    Returns:
        dict:

    Raises:
        ValueError: If the Ensembl lookup fails or returns nothing.
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
    Filter out non-chromosomal accessions under assemblies.
    Args: A list of assemblies
    Returns: A list of filtered assemblies
    """
    filtered_assemblies = []
    for assembly in related_assemblies:
        if assembly.get("accession", "unknown").startswith("NC_"):
            filtered_assemblies.append(assembly)
    return filtered_assemblies


def filter_gene_products(accession_base, related_genes):
    """
    Filter out accessions from every gene's products (transcripts and proteins).
    Keep the products with MANE Select tag or the same as the queried accession.
    Args:
        accession_base (str):
        related_genes (list): related genes
    Returns: A list of related genes of filtered products
    """
    filtered_genes = []

    for gene in related_genes:
        filtered_transcripts = []
        for transcript in gene.get("transcripts", []):
            matching_providers = [
                provider
                for provider in transcript.get("providers", [])
                if accession_base
                in {
                    (provider.get("transcript_accession") or "").split(".")[
                        0
                    ],
                    (provider.get("protein_accession") or "").split(".")[0],
                }
            ]
            if "tag" in transcript or matching_providers:
                filtered_transcripts.append(transcript)

        gene["transcripts"] = filtered_transcripts
        filtered_genes.append(gene)
    return filtered_genes


def filter_related(accession_base, related):

    filtered = {}
    filtered_assemblies = filter_assemblies(related.get("assemblies", []))
    if filtered_assemblies:
        filtered["assemblies"] = filtered_assemblies
    filterd_gene_products = filter_gene_products(
        accession_base, related.get("genes", [])
    )
    if filterd_gene_products:
        filtered["genes"] = filterd_gene_products
    return filtered


def _get_gene_related(gene_ids):
    client = NCBIClient(timeout=DEFAULT_TIMEOUT)
    product_response = client.get_gene_id_product_report(gene_ids)
    parsed_product = datasets._parse_product_report(product_response)
    dataset_response = client.get_gene_id_dataset_report(gene_ids)
    parsed_dataset = datasets._parse_dataset_report(dataset_response)
    taxname = product_response.get("taxname")
    return taxname, datasets._merge_datasets(parsed_dataset, parsed_product)


def _parse_genome_annotation_report(genome_annotation_report):
    gene_ids = []
    for report in genome_annotation_report.get("reports", []):
        annotation = report.get("annotation", {})
        gene_id = annotation.get("gene_id")
        taxon_name = annotation.get("taxname")
        if gene_id is not None:
            gene_ids.append(gene_id)
    if gene_ids:
        taxon_name, related = _get_gene_related(gene_ids)
        return taxon_name, related

    return None, {}


def _get_related_by_chr_location(accession, locations):
    taxon_name = None
    related = {}
    assembly_accession = _get_assembly_accession(accession)
    if not assembly_accession:
        raise ValueError(
            f"Assembly accession could not be determined for {accession}"
        )

    client = NCBIClient(timeout=DEFAULT_TIMEOUT)
    annotation_response = client.get_genome_annotation_report(
        assembly_accession, _validate_locations(accession, locations)
    )
    taxon_name, related = _parse_genome_annotation_report(annotation_response)

    return taxon_name, related


def _get_assembly_accession(accession):
    client = NCBIClient(timeout=DEFAULT_TIMEOUT)
    return client.get_assembly_accession(accession)


def _get_related_by_accession_from_ncbi(accession):
    """
    Fetch related genomic and product sequences from NCBI endpoint using a RefSeq accession.
    This function contacts the NCBI Datasets API to gather both genomic and product-level
    annotations for a given transcript or protein accession.

    Args:
        accession (str): A RefSeq accession.
        timeout (int, optional): Timeout, defaults to DEFAULT_TIMEOUT.
    Returns:
        tuple:
            - taxon_name (str or None): The organism name associated with the accession.
            - related (dict or None): A dictionary of related sequences; otherwise None.
    Raises:
        RuntimeError: If the NCBI Datasets API is unavailable or returns an invalid response.
    """

    related_from_ncbi = {}

    client = NCBIClient(timeout=DEFAULT_TIMEOUT)
    product_report = client.get_accession_product_report(accession)
    dataset_report = client.get_accession_dataset_report(accession)

    parsed_products = datasets._parse_product_report(product_report)
    parsed_dataset = datasets._parse_dataset_report(dataset_report)

    related_from_ncbi = datasets._merge_datasets(
        parsed_dataset, parsed_products
    )
    taxname = product_report.get("taxname")
    return taxname, related_from_ncbi


def _get_related_by_ensembl_id(accession, moltype):
    """
    Given an Ensembl accession, return filtered related from Ensembl and NCBI
    Args:
        accession (str): An Ensembl accession.

    Returns:
        related (dict): A dictionary containing related sequences from Ensembl, NCBI
    """
    accession_base = accession.split(".")[0]

    # Get taxon_name and related from ensembl
    ensembl_related = fetch_ensembl_gene_info(accession_base, moltype)

    ncbi_related = {}
    # Get related from ncbi using gene symbol and taxname
    gene_symbol = ensembl_related.get("name")
    taxname = ensembl_related.get("taxon_name")
    if gene_symbol and taxname:
        ncbi_related = _get_related_by_gene_symbol_from_ncbi(
            gene_symbol, taxname
        )

    related = _merge(ensembl_related, ncbi_related)
    if taxname and taxname.upper() == HUMAN_TAXON:
        return filter_related(accession_base, related)

    return related


def _get_related_by_ncbi_id(accession, moltype, locations):
    related = {}
    accession_base = accession.split(".")[0]

    # Get taxon_name and related from ncbi datasets
    if moltype in [MoleculeType.RNA, MoleculeType.PROTEIN]:
        taxon_name, ncbi_related = _get_related_by_accession_from_ncbi(
            accession
        )
    elif moltype == MoleculeType.DNA and "NC_" in accession:
        taxon_name, ncbi_related = _get_related_by_chr_location(
            accession, locations
        )
    else:
        raise NameError(f"Could not retrieve {accession} from NCBI.")

    # Get related from ensembl using ensembl gene id.
    if ncbi_related:
        ensembl_genes_id = [
            provider["accession"]
            for gene in ncbi_related.get("genes", [])
            for provider in gene.get("providers", [])
            if provider.get("name") == DataSource.ENSEMBL
        ]
        all_genes = []
        all_assemblies = []
        for ensembl_gene in ensembl_genes_id:
            ensembl_gene_related = fetch_ensembl_gene_info(
                ensembl_gene, moltype="dna"
            )
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
    # add doc string

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
                source: DataSource.ENSEMBL, DataSource.NCBI, or "other"
                moltype: "dna", "rna", "protein", "unknown" and others
    """
    seq_id = seq_id.strip()

    source, moltype = parse_ensembl_id(seq_id)
    if source:
        return source, moltype

    source, moltype = parse_ncbi_id(seq_id)
    if source:
        return source, moltype

    return DataSource.OTHER, MoleculeType.UNKNOWN


def get_related(accession, locations="0"):
    """
    Retrieve related assembly/gene/transcript/protein information based
    on accession or gene symbol
    Args:
        accession (str): A sequence accession (e.g., RefSeq, Ensembl ID) or a human gene symbol.
        locations, optional (str): A point or a range on chromosome, in the format of
            '10000;120000_130000'. Defaults to '0'.

    Returns:
        related (dict): A dictionary containing related information retrieved from Ensembl, NCBI,

    Raises:
        NameError: If the given accession is not from NCBI RefSeq or ENSEMBL.
    """
    accession = accession.upper()

    source, moltype = detect_sequence_source(accession)
    if source == DataSource.ENSEMBL and moltype != MoleculeType.UNKNOWN:
        return _get_related_by_ensembl_id(accession, moltype)
    if source == DataSource.NCBI and moltype != MoleculeType.UNKNOWN:
        return _get_related_by_ncbi_id(accession, moltype, locations)
    if source == DataSource.OTHER:
        return _get_related_by_gene_symbol(accession)
    raise NameError(f"Could not retrieve related for {accession}.")
