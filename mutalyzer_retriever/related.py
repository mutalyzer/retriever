import json
import re
import requests
from copy import deepcopy
from mutalyzer_retriever.util import DataSource, MoleculeType, HUMAN_TAXON, DEFAULT_TIMEOUT
from mutalyzer_retriever.configuration import cache_url, settings
from mutalyzer_retriever.request import Http400, request
from mutalyzer_retriever.reference import GRCH37
from mutalyzer_retriever.client import NCBIClient, EnsemblClient
from mutalyzer_retriever.parsers.ensembl_lookup import _parse_ensembl_gene_lookup_json


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


def _validate_locations(accession:str, locations):
    """
    Check if input locations are in the format of 'start_end' or single points,
    and return a normalized string like: accession:start_end;accession:start_end
    """
    pattern_range = re.compile(r'^\d+_\d+$')
    pattern_point = re.compile(r'^\d+$')

    valid_locations = []

    if locations == "0":
        raise ValueError(f"Unkown location on chromosomes.")

    for item in locations.split(';'):
        item = item.strip()
        if pattern_range.match(item):
            try:
                raw_start, raw_end = map(int, item.split('_'))
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

def merge_provider(ensembl_provider, ncbi_provider):
    if not ncbi_provider:
        return ensembl_provider
    else:
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
    return ensembl_related.get("assemblies", []) + ncbi_related.get("assemblies", [])


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
            if ncbi_match(ncbi_t, ensembl_accession):
                transcript["providers"] = merge_provider(ensembl_entry, ncbi_t.get("providers"))
                merged.append(transcript)
                matched = True
                break
        if not matched:
            merged.append({"providers":[ensembl_entry]})
    return merged


def _merge_gene(ensembl_related, ncbi_related):
    "Merge two lists of related genes gathered from ensembl and ncbi"
    ensembl_gene_name = ensembl_related.get("name")
    ensembl_gene_accession = ensembl_related.get("accession")
    if not (ensembl_gene_name and ensembl_gene_accession):
        return ncbi_related.get("genes", [])

    merged = []
    for ncbi_gene in ncbi_related.get("genes", []):
        if ncbi_gene.get("name") == ensembl_gene_name:
            # shape ensembl gene data and merge gene info from two sources
            gene = deepcopy(ncbi_gene)
            ensembl_entry = {"accession": ensembl_gene_accession,
                              "name": DataSource.ENSEMBL}
            gene["providers"] = merge_provider(ensembl_entry, ncbi_gene.get("providers"))
            # merge transcripts
            gene["transcripts"] = _merge_transcripts(ensembl_related, ncbi_gene)
            merged.append(gene)
    return merged


def _merge(ensembl_related, ncbi_related):
    return {
        "assemblies": _merge_assemblies(ensembl_related, ncbi_related),
        "genes": _merge_gene(ensembl_related, ncbi_related),
    }


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
        return
    
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
        ncbi_transcript = transcript.get("accession_version")
        ncbi_desc = transcript.get("name")
        ncbi_protein = transcript.get("protein", {})
        ncbi_protein_acc = ncbi_protein.get("accession_version")
        ensembl_transcript = transcript.get("ensembl_transcript")
        ensembl_protein = ncbi_protein.get("ensembl_protein")

        ncbi_transcript = {
                "name": DataSource.NCBI,
                "transcript_accession": ncbi_transcript,
                "protein_accession": ncbi_protein_acc,
                "description": ncbi_desc,
            }

        ensembl_transcript = {
                "name": DataSource.ENSEMBL,
                "transcript_accession": ensembl_transcript,
                "protein_accession": ensembl_protein,
            }

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
        return
    
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
        return

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
    related = {}
    taxnme = None
    if genomic_related is None and product_related is None:
        return None, None

    if genomic_related:
        merged_genes = []
        taxnme = genomic_related.get("taxname")
        if genomic_related.get("assemblies"):
            related["assemblies"] = genomic_related.get("assemblies")


        for genomic_gene in genomic_related.get("genes", []):
            symbol = genomic_gene.get("name")
            for product_gene in product_related.get("genes", []):
                if symbol == product_gene.get("name"):
                    transcripts = product_gene.get("transcripts", [])
                    gene_products = {"transcripts": transcripts}

                    merged_gene =  genomic_gene | gene_products
                    merged_genes.append(merged_gene)
                if merged_genes:
                    related["genes"] = merged_genes
        # Sort genes alphabetically by name
        related["genes"].sort(key=lambda g: g.get("name", ""))

        return taxnme, related
    return None, None



def _get_related_by_gene_symbol_from_ncbi(gene_symbol, taxon_name=HUMAN_TAXON):
    """
    Given a gene symbol, return a set of related sequence accessions (genomic and/or products).
    Returns (taxon_name, related_dict), or (None, None) if nothing found.
    """
    if not gene_symbol:
        return None, None
    related = {}

    client = NCBIClient(timeout=DEFAULT_TIMEOUT)

    dataset_response = client.get_gene_symbol_dataset_report(gene_symbol, taxon_name)
    parsed_dataset = _parse_dataset_report(
        dataset_response
    )
    product_response = client.get_gene_symbol_product_report(gene_symbol, taxon_name)
    parsed_product = _parse_product_report(
        product_response
    )
    related = _merge_datasets(parsed_dataset, parsed_product)
    return related


def _get_related_by_gene_symbol_from_ensembl(
    gene_symbol
):
    if not gene_symbol:
        return None
    client = EnsemblClient(timeout=DEFAULT_TIMEOUT)
    gene_lookup_response = client.lookup_symbol(gene_symbol)
    return _parse_ensembl_gene_lookup_json(gene_lookup_response)
 

def _get_related_by_gene_symbol(gene_symbol):
    _, ncbi_related = _get_related_by_gene_symbol_from_ncbi(
        gene_symbol, taxon_name=HUMAN_TAXON
    )
    ensembl_related = _get_related_by_gene_symbol_from_ensembl(
        gene_symbol
    )
    related = _merge(ensembl_related, ncbi_related)
    
    related = filter_related(gene_symbol, related)
    return related


def _parse_ensembl_gene_lookup_json(response):
    transcripts = []
    for transcript in response.get("Transcript", []):
        transcript_id = transcript.get("id")
        transcript_version = transcript.get("version")

        if not transcript_id or not transcript_version:
            continue
        t = {
                "transcript_accession": f"{transcript_id}.{transcript_version}",
                "description": transcript.get("display_name"),
            }

        translation = transcript.get("Translation", {})
        protein_id = translation.get("id")
        protein_version = translation.get("version")
        if protein_id and protein_version:
            protein_accession = f"{protein_id}.{protein_version}"
            t["protein_accession"] = protein_accession

        transcripts.append(t)

    gene_symbol = response.get("display_name")
    taxon_name = response.get("species").replace("_", " ").upper()    
    output = {
        "taxon_name": taxon_name,
        "name":gene_symbol
    }
    if transcripts:
        output["transcripts"] = transcripts
    gene_id = response.get("id")
    if gene_id:
        output["accession"] = gene_id
    return output


def _parse_ensembl_transcript_lookup_json(ensembl_transcript_json):
    ensembl_gene_id = ensembl_transcript_json.get("Parent")
    client = EnsemblClient(timeout=DEFAULT_TIMEOUT)
    ensembl_gene_json = client.lookup_id(ensembl_gene_id)
    return _parse_ensembl_lookup_json(ensembl_gene_json)


def _parse_ensembl_protein_lookup_json(ensembl_protein_json):
    ensembl_transcript_id = ensembl_protein_json.get("Parent")
    client = EnsemblClient(timeout=DEFAULT_TIMEOUT)
    ensembl_transcript_json = client.lookup_id(ensembl_transcript_id)
    return _parse_ensembl_lookup_json(ensembl_transcript_json)


def _parse_ensembl_lookup_json(ensembl_json):
    """
    Dispatch parsing based on Ensembl object type.
    """
    obj_type = ensembl_json.get("object_type")
    if obj_type == "Gene" and ensembl_json.get("Transcript"):
        return _parse_ensembl_gene_lookup_json(ensembl_json)
    if obj_type == "Transcript" and ensembl_json.get("Parent"):
        return _parse_ensembl_transcript_lookup_json(ensembl_json)
    if obj_type == "Translation" and ensembl_json.get("Parent"):
        return _parse_ensembl_protein_lookup_json(ensembl_json)

    raise ValueError(f"Unsupported or malformed Ensembl object: {obj_type}")



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
    if moltype not in {MoleculeType.DNA, MoleculeType.RNA, MoleculeType.PROTEIN}:
        raise ValueError(f"Unsupported molecule type: {moltype}")

    client = EnsemblClient(timeout=DEFAULT_TIMEOUT)
    expand_flag = 0 if moltype == MoleculeType.PROTEIN else 1

    ensembl_lookup_json = client.lookup_id(
        accession_base,
        expand=expand_flag
    )

    if not ensembl_lookup_json:
        raise ValueError(f"No data returned for {accession_base} from ENSEMBL.")

    return _parse_ensembl_lookup_json(ensembl_lookup_json)


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
                if accession_base in {
                    (provider.get("transcript_accession") or "").split(".")[0],
                    (provider.get("protein_accession") or "").split(".")[0]
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
    filterd_gene_products = filter_gene_products(accession_base, related.get("genes", []))
    if filterd_gene_products:
        filtered["genes"] = filterd_gene_products
    return filtered


def _get_gene_related(gene_ids):
    client = NCBIClient(timeout=DEFAULT_TIMEOUT)
    product_response = client.get_gene_id_product_report(gene_ids)
    parsed_product = _parse_product_report(product_response)
    dataset_response = client.get_gene_id_dataset_report(gene_ids)
    parsed_dataset = _parse_dataset_report(dataset_response)
    return _merge_datasets(parsed_dataset, parsed_product)


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
    annotation_response = client.get_genome_annotation_report(assembly_accession, _validate_locations(accession, locations))
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

    related = {}
    products = None
    assemblies = None
    taxon_name = None

    client = NCBIClient(timeout=DEFAULT_TIMEOUT)
    product_report = client.get_accession_product_report(accession)
    dataset_report = client.get_accession_dataset_report(accession)

    _, products = _parse_product_report(product_report)
    taxon_name, assemblies = _parse_dataset_report(dataset_report)

    related = _merge_datasets(assemblies, products)
    return taxon_name, related


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
    taxname = ensembl_related.get("taxon_name")
    ncbi_related = {}
    # Get related from ncbi using gene symbol and taxname
    genes = ensembl_related.get("genes")
    if genes and isinstance(genes, list) and genes[0].get("name"):
        gene_symbol = genes[0]["name"]
        _, ncbi_related = _get_related_by_gene_symbol_from_ncbi(gene_symbol, taxname)

    related = _merge(ensembl_related, ncbi_related)
    if taxname and taxname.upper() == HUMAN_TAXON:
        return filter_related(accession_base, related)
    
    return related



def _get_related_by_ncbi_id(accession, moltype, locations):
    related = {}
    accession_base = accession.split(".")[0]

    # Try to get taxon_name and related from ncbi datasets
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

        for ensembl_gene in ensembl_genes_id:
            ensembl_gene_related = fetch_ensembl_gene_info(ensembl_gene, moltype="dna")
            if ensembl_gene_related:
                related = _merge(ensembl_gene_related, ncbi_related)
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
    ensembl_pattern = re.compile(r"^ENS[A-Z]*?(FM|GT|G|T|P|R|E)\d{11}(?:\.\d+)?$")
    match = ensembl_pattern.match(seq_id)
    if match:
        prefix = match.group(1)
        moltype = ensembl_feature_map.get(prefix, MoleculeType.UNKNOWN)
        return DataSource.ENSEMBL, moltype
    return DataSource.ENSEMBL, MoleculeType.UNKNOWN


def parse_ncbi_id(seq_id):
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


# get rid of timeout
def get_related(accession, locations=None):
    """
    Retrieve related assembly/gene/transcript/protein information based
    on accession or gene symbol
    Args:
        accession (str): A sequence accession (e.g., RefSeq, Ensembl ID) or a human gene symbol.
        locations (list of list[int, int], optional): Genomic location ranges
            in the form [[start, end], ...]. Defaults to '0'.
    Returns:
        related (dict): A dictionary containing related information retrieved from Ensembl, NCBI,

    Raises:
        NameError: If the given accession is not from NCBI RefSeq or ENSEMBL.
    """

    if locations is None:
        locations = "0"
    accession = accession.upper()
 
    source, moltype = detect_sequence_source(accession)
    if source == DataSource.ENSEMBL and moltype != MoleculeType.UNKNOWN:
        return _get_related_by_ensembl_id(accession, moltype)
    if source == DataSource.NCBI and moltype != MoleculeType.UNKNOWN:
        return _get_related_by_ncbi_id(
            accession, moltype, locations=locations
        )
    if source == DataSource.OTHER:
        return _get_related_by_gene_symbol(accession)
    raise NameError(f"Could not retrieve related for {accession}.")
