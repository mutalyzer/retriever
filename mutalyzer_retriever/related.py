import json
import re
from urllib.parse import quote
import requests
from mutalyzer_retriever.util import DataSource, MoleculeType, HUMAN_TAXON, DEFAULT_TIMEOUT
from mutalyzer_retriever.configuration import cache_url, settings
from mutalyzer_retriever.request import Http400, request
from mutalyzer_retriever.reference import GRCH37

class BaseAPIClient:
    """Base class for API clients"""
    
    def __init__(self, base_url: str, timeout: int = DEFAULT_TIMEOUT):
        self.base_url = base_url
        self.timeout = timeout
    
    def _make_request(self, url: str, params: dict | None = None):
        """Make HTTP request"""
        response = request(url=url, params=params, timeout=self.timeout)
        return json.loads(response)

class NCBIClient(BaseAPIClient):
    """Client for NCBI API operations"""
    
    def __init__(self, timeout: int = DEFAULT_TIMEOUT):
        base_url = settings.get("NCBI_DATASETS_API")
        super().__init__(base_url, timeout)
    
    def get_accession_dataset_report(self, accession: str):
        """Fetch dataset report for given accession"""
        url = f"{self.base_url}/gene/accession/{accession}/dataset_report"
        return self._make_request(url)
    
    def get_accession_product_report(self, accession: str):
        """Fetch product report for given accession"""
        url = f"{self.base_url}/gene/accession/{accession}/product_report"
        return self._make_request(url)
    
    def get_gene_id_dataset_report(self, gene_ids: list[str]):
        """Fetch dataset report for gene IDs"""
        gene_id_str = quote(",".join(map(str, gene_ids)))
        url = f"{self.base_url}/gene/id/{gene_id_str}/dataset_report"
        return self._make_request(url)
    
    def get_gene_id_product_report(self, gene_ids: list[str]):
        """Fetch product report for gene IDs"""
        gene_id_str = quote(",".join(map(str, gene_ids)))
        url = f"{self.base_url}/gene/id/{gene_id_str}/product_report"
        return self._make_request(url)
    
    def get_gene_symbol_dataset_report(self, gene_symbol: str, taxname: str = HUMAN_TAXON):
        """Fetch dataset report for gene symbol"""
        taxname_url_str = quote(taxname, safe="")
        url = f"{self.base_url}/gene/symbol/{gene_symbol}/taxon/{taxname_url_str}/dataset_report"
        return self._make_request(url)
    
    def get_gene_symbol_product_report(self, gene_symbol: str, taxname: str = HUMAN_TAXON):
        """Fetch product report for gene symbol"""
        taxname_url_str = quote(taxname, safe="")
        url = f"{self.base_url}/gene/symbol/{gene_symbol}/taxon/{taxname_url_str}/product_report"
        return self._make_request(url)
    
    def get_assembly_accession(self, accession: str):
        """Get assembly accession for sequence accession"""
        url = f"{self.base_url}/genome/sequence_accession/{accession}/sequence_assemblies"
        response = self._make_request(url)
        accessions = response.get("accessions")
        if isinstance(accessions, list) and accessions:
            return accessions[0]
        return None
    
    def get_genome_annotation_report(self, assembly_accession: str, locations: list[str]):
        """Get genome annotation report for assembly and locations"""
        url = f"{self.base_url}/genome/accession/{assembly_accession}/annotation_report"
        params = [("locations", loc) for loc in locations]
        return self._make_request(url, params)


class EnsemblClient(BaseAPIClient):
    """Client for Ensembl API operations"""
    
    def __init__(self, timeout: int = DEFAULT_TIMEOUT):
        super().__init__(settings.get("ENSEMBL_API"), timeout)
    
    def lookup_symbol(self, gene_symbol: str, species: str = "homo_sapiens"):
        """Lookup gene by symbol"""
        url = f"{self.base_url}/lookup/symbol/{species}/{gene_symbol}?content-type=application/json;expand=1"
        return self._make_request(url)
    
    def lookup_id(self, accession_base: str, expand: int = 1):
        """Lookup by Ensembl ID"""
        url = f"{self.base_url}/lookup/id/{accession_base}?content-type=application/json;expand={expand}"
        return self._make_request(url)


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


def clean_dict(d):
    """Remove keys with None, empty string, or empty list values."""
    return {k: v for k, v in d.items() if v not in (None, "", [])}


def _validate_locations(accession:str, locations: str):
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



def _merge_assemblies(ensembl_assemblies: list, ncbi_assemblies: list):
    "Merge two lists of assemblies gathered from ensembl and ncbi"
    return ensembl_assemblies + ncbi_assemblies


def _merge_transcripts(ensembl_transcripts: list, ncbi_transcripts: list):
    "Merge two lists of transcripts gathered from ensembl and ncbi"
    merged_transcripts = []
    for transcript in ncbi_transcripts:
        merged_transcripts.append(transcript.copy())


    for ensembl_t in ensembl_transcripts:
        ensembl_acc = ensembl_t.get("transcript_accession")
        matched = False

        for transcript in merged_transcripts:
            providers = transcript.get("providers", [])
            for i, provider in enumerate(providers):
                if provider.get("name") == "ENSEMBL" and provider.get("transcript_accession") == ensembl_acc:
                    providers[i] = clean_dict(ensembl_t)
                    matched = True
                    break  
            if matched:
                break  

        if not matched:
            merged_transcripts.append({"providers": [clean_dict(ensembl_t)]})

    return merged_transcripts


def _merge_gene(ensembl_gene, ncbi_genes):
    "Merge two lists of related genes gathered from ensembl and ncbi"
    genes = []
    for gene in ncbi_genes:
        gene_name = gene.get("name", "")
        if gene_name == next(g for g in ensembl_gene if g != "providers"):
            ensembl_gene_provider = ensembl_gene.get("providers")
            if ensembl_gene_provider and all(
                ensembl_gene_provider != p for p in gene.get("providers", [])
            ):
                gene["providers"].append(ensembl_gene_provider)
        gene["transcripts"] = _merge_transcripts(
            ensembl_gene.get(gene_name, []), gene.get("transcripts", [])
        )
        genes.append(gene)
    return genes


def _merge(ensembl_related, ncbi_related):
    # Merge gene sets related from ENSEMBL and NCBI
    related = {
        "assemblies": _merge_assemblies(
            ensembl_related.get("assemblies", []),
            ncbi_related.get("assemblies", []),
        ),
        "genes": _merge_gene(ensembl_related, ncbi_related.get("genes", [])),
    }
    return related


def _parse_assemblies(dataset_report):
    taxname = None
    assemblies = []

    gene = dataset_report.get("gene", {})
    taxname = taxname or gene.get("taxname")
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
    if taxname == "Homo sapiens" and sequence_name:
        grch37_acc = GRCH37.get(sequence_name)
        if grch37_acc:
            grch37_entry = {
                    "name": "GRCh37.p13",
                    "accession": grch37_acc,
            }
            assemblies.append(grch37_entry)
    return assemblies


def _parse_genes(dataset_report):

    gene = dataset_report.get("gene", {})
    symbol = gene.get("symbol")
    description = gene.get("description")
    ensembl_ids = gene.get("ensembl_gene_ids")
    hgnc_id = None

    # Reference and provider info
    ncbi_id = gene.get("gene_id")
    hgnc_id_raw = gene.get("nomenclature_authority", {}).get("identifier")
    if hgnc_id_raw and ":" in hgnc_id_raw:
        hgnc_id = hgnc_id_raw.split(":")[1] if "HGNC" in hgnc_id_raw else None

    ref_standards = gene.get("reference_standards", [])
    ref_accession = None
    for ref in ref_standards:
        ref_accession = ref.get("gene_range", {}).get("accession_version")

    ncbi_gene = clean_dict(
        {
            "name": "NCBI",
            "id": ncbi_id,
            "accession": ref_accession,
        }
    )

    ensembl_gene = clean_dict(
        {
            "name": "ENSEMBL",
            "accession": ensembl_ids[0] if ensembl_ids else None,
        }
    )

    providers = [p for p in (ncbi_gene, ensembl_gene) if len(p) > 1]
    gene_entry = clean_dict(
        {
            "name": symbol,
            "hgnc_id": hgnc_id,
            "description": description,
            "providers": providers,
        }
    )
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

        ncbi_transcript = clean_dict(
            {
                "name": "NCBI",
                "transcript_accession": ncbi_transcript,
                "protein_accession": ncbi_protein_acc,
                "description": ncbi_desc,
            }
        )

        ensembl_transcript = clean_dict(
            {
                "name": "ENSEMBL",
                "transcript_accession": ensembl_transcript,
                "protein_accession": ensembl_protein,
            }
        )

        providers = [
            p for p in (ncbi_transcript, ensembl_transcript) if len(p) > 1
        ]
        if not providers:
            continue

        product_entry = clean_dict(
            {
                "providers": providers,
                "tag": transcript.get("select_category"),
            }
        )
        gene_products.append(product_entry)
    return gene_products


def _parse_dataset_report(dataset_report):
    """
    Parses a gene dataset report from the NCBI Datasets API.

    Args:
        json_report (dict): JSON object returned from the NCBI Datasets API.

    Returns:
        tuple:
            - taxname (str): The scientific name of the organism (e.g., "Homo sapiens").
            - data (dict): A dictionary containing:
                - "assemblies" (list): A list of assembly dictionaries.
                - "genes" (list): A list of gene dictionaries.
    """
    assemblies = []
    genes = []
    taxname = None

    for report in dataset_report.get("reports", []):
        taxname = taxname or report.get("gene", {}).get("taxname", "unknown")
        # Extract genomic assemblies and genes
        assemblies = _parse_assemblies(report)
        gene = _parse_genes(report)
        if gene:
            genes.append(gene)
    return taxname, {"assemblies": assemblies, "genes": genes}


def _parse_product_report(product_report):
    """
    Parse NCBI product report JSON into a dictionary.
    Args:
        data (dict): JSON data from ncbi product API.
    Returns:
        dict: Mapping gene_symbol -> list of transcript info dicts
    """
    genes_dict = {}
    taxname = None

    for report in product_report.get("reports", []):
        taxname = report.get("taxname")
        product = report.get("product", {})
        symbol = product.get("symbol")
        gene_products = _parse_transcripts(report)
        if gene_products:
            genes_dict[symbol] = gene_products

    return taxname, genes_dict


def _merge_datasets(genomic_related, product_related):
    """
    Merges genomic and product-related from Datasets.
    """
    if genomic_related is None and product_related is None:
        return
    if genomic_related:
        related = {
            "assemblies": [
                clean_dict(a)
                for a in genomic_related.get("assemblies", [])
                if clean_dict(a)
            ],
            "genes": [],
        }

        for gene in genomic_related.get("genes", []):
            symbol = gene.get("name")
            transcripts = product_related.get(symbol, [])

            gene_info = clean_dict(
                {
                    "name": symbol,
                    "hgnc_id": gene.get("hgnc_id"),
                    "description": gene.get("description"),
                    "transcripts": transcripts,
                    "providers": gene.get("providers"),
                }
            )

            if gene_info:
                related["genes"].append(gene_info)

        # Sort genes alphabetically by name
        related["genes"].sort(key=lambda g: g.get("name", ""))

        return related
    return None


def filter_report_from_other_genes(gene_symbol: str, reports: dict):
    # NCBI datasets would return related genes when query by gene name
    # e.g, query with CYP2D6 would get info of CYP2D7
    # https://api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol/CYP2D6D%2CCYP2D6/taxon/9606/product_report
    if not reports:
        return None
    for report in reports.get("reports"):
        for key, value in report.items():
            if value.get("symbol") == gene_symbol.upper():
                return {"reports": [{key: value}]}
    return None


def _get_related_by_gene_symbol_from_ncbi(gene_symbol, taxname="Homo Sapiens"):
    """
    Given a gene symbol, return a set of related sequence accessions (genomic and/or products).
    Returns (taxname, related_dict), or (None, None) if nothing found.
    """
    if not gene_symbol:
        return None, None
    related = {}
    products = None
    assemblies = None

    client = NCBIClient(timeout=DEFAULT_TIMEOUT)

    dataset_response = client.get_gene_symbol_dataset_report(gene_symbol, taxname)
    taxname, assemblies = _parse_dataset_report(
        dataset_response
    )
    product_response = client.get_gene_symbol_product_report(gene_symbol, taxname)
    taxname, products = _parse_product_report(
        product_response
    )

    related = _merge_datasets(assemblies, products)
    return taxname, related


def _get_related_by_gene_symbol_from_ensembl(
    gene_symbol, taxname=HUMAN_TAXON
):
    if not gene_symbol:
        return None
    client = EnsemblClient(timeout=DEFAULT_TIMEOUT)
    gene_lookup_response = client.lookup_symbol(gene_symbol)
    taxname, ensembl_gene_id, gene_name, ensemble_transcripts = _parse_ensembl_lookup_json(gene_lookup_response)
    return taxname, {
        gene_name: ensemble_transcripts,
        "providers": {"name": "ENSEMBL", "accession": ensembl_gene_id},
    }


def _get_related_by_gene_symbol(gene_symbol):
    _, ncbi_related = _get_related_by_gene_symbol_from_ncbi(
        gene_symbol, taxname=HUMAN_TAXON
    )
    _, ensembl_related = _get_related_by_gene_symbol_from_ensembl(
        gene_symbol, taxname=HUMAN_TAXON
    )
    related = _merge(ensembl_related, ncbi_related)
    related = filter_related(gene_symbol, related)
    return related


def _parse_ensembl_gene_lookup_json(ensembl_gene_json):
    ensembl_related_transcripts = []
    gene_symbol = ensembl_gene_json.get("display_name")
    gene_id = ensembl_gene_json.get("id")
    taxname = ensembl_gene_json.get("species").replace("_", " ")

    for transcript in ensembl_gene_json.get("Transcript", []):
        translation = transcript.get("Translation", {})

        transcript_id = transcript.get("id")
        transcript_version = transcript.get("version")
        protein_id = translation.get("id")
        protein_version = translation.get("version")

        if not transcript_id or not transcript_version:
            continue
        if protein_id and protein_version:
            protein_accession = f"{protein_id}.{protein_version}"
        else:
            protein_accession = None

        ensembl_related_transcripts.append(
            {
                "name": "ENSEMBL",
                "transcript_accession": f"{transcript_id}.{transcript_version}",
                "protein_accession": protein_accession,
                "description": transcript.get("display_name"),
            }
        )

    return taxname, gene_id, gene_symbol, ensembl_related_transcripts


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
    has_transcripts = bool(ensembl_json.get("Transcript"))
    has_parent = bool(ensembl_json.get("Parent"))

    if obj_type == "Gene" and has_transcripts:
        return _parse_ensembl_gene_lookup_json(ensembl_json)
    elif obj_type == "Transcript" and has_parent:
        return _parse_ensembl_transcript_lookup_json(ensembl_json)
    elif obj_type == "Translation" and has_parent:
        return _parse_ensembl_protein_lookup_json(ensembl_json)

    raise ValueError(f"Unsupported or malformed Ensembl object: {obj_type}")



def fetch_ensembl_gene_info(accession_base: str, moltype: str) -> tuple[str, dict]:
    """
    Fetch gene-related information for a given Ensembl accession.

    Args:
        accession_base (str): Ensembl accession, support for gene, transcript, protein.
        moltype (str): Molecule type — features parsed from ensembl accession prefix.

    Returns:
        tuple: (taxname, related_info_dict)
            - taxname (str): Taxonomic name (e.g., "Homo sapiens")
            - related_info_dict (dict): Dictionary including gene and its products.

    Raises:
        ValueError: If the Ensembl lookup fails or returns nothing.
    """
    if moltype not in {MoleculeType.DNA, MoleculeType.RNA, MoleculeType.PROTEIN}:
        raise ValueError(f"Unsupported molecule type: {moltype}")

    client = EnsemblClient(timeout=DEFAULT_TIMEOUT)
    ensembl_lookup_json = client.lookup_id(
        accession_base,
        expand=0 if moltype == "protein" else 1
    )

    if not ensembl_lookup_json:
        raise ValueError(f"No data returned for {accession_base} from ENSEMBL.")

    parsed = _parse_ensembl_lookup_json(ensembl_lookup_json)
    if not parsed:
        raise ValueError(f"Failed to parse response for {accession_base}.")

    taxname, ensembl_id, gene, related_transcripts = parsed

    return taxname, {
        gene: related_transcripts,
        "providers": {
            "name": "ENSEMBL",
            "accession": ensembl_id
        }
    }




def filter_assemblies(related_assemblies: list):
    """
    Filter out non-chromosomal accessions under assemblies.
    Args: A list of assemblies
    Returns: A list of filtered assemblies
    """
    filtered_assemblies = []
    for assembly in related_assemblies:
        if assembly.get("accession", "unknown").startswith("NC_"):
            filtered_assemblies.append(clean_dict(assembly))
    return filtered_assemblies


def filter_gene_products(accession_base, related_genes: list):
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
        filtered_genes.append(clean_dict(gene))
    return filtered_genes


def filter_related(accession_base, related):
    return {
        "assemblies": filter_assemblies(related.get("assemblies", [])),
        "genes": filter_gene_products(
            accession_base, related.get("genes", [])
        ),
    }


def _get_gene_related(gene_ids):
    client = NCBIClient(timeout=DEFAULT_TIMEOUT)
    product_response = client.get_gene_id_product_report(gene_ids)
    taxname, products = _parse_product_report(product_response)
    dataset_response = client.get_gene_id_dataset_report(gene_ids)
    taxname, assemblies = _parse_dataset_report(dataset_response)
    return taxname, _merge_datasets(assemblies, products)


def _parse_genome_annotation_report(genome_annotation_report):
    gene_ids = []
    for report in genome_annotation_report.get("reports", []):
        annotation = report.get("annotation", {})
        gene_id = annotation.get("gene_id")
        taxname = annotation.get("taxname")
        if gene_id is not None:
            gene_ids.append(gene_id)
    if gene_ids:
        taxname, related = _get_gene_related(gene_ids)
        return taxname, related
    return None, {}


def _get_related_by_chr_location(accession, locations):
    taxname = None
    related = {}
    assembly_accession = _get_assembly_accession(accession)
    if not assembly_accession:
        raise ValueError(
            f"Assembly accession could not be determined for {accession}"
        )

    client = NCBIClient(timeout=DEFAULT_TIMEOUT)
    annotation_response = client.get_genome_annotation_report(assembly_accession, _validate_locations(accession, locations))
    taxname, related = _parse_genome_annotation_report(annotation_response)

    return taxname, related


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
            - taxname (str or None): The organism name associated with the accession.
            - related (dict or None): A dictionary of related sequences; otherwise None.
    Raises:
        RuntimeError: If the NCBI Datasets API is unavailable or returns an invalid response.
    """

    related = {}
    products = None
    assemblies = None
    taxname = None

    client = NCBIClient(timeout=DEFAULT_TIMEOUT)
    product_report = client.get_accession_product_report(accession)
    dataset_report = client.get_accession_dataset_report(accession)

    products = _parse_product_report(product_report)
    taxname, assemblies = _parse_dataset_report(dataset_report)

    related = _merge_datasets(assemblies, products)
    return taxname, related


def _get_related_by_ensembl_id(accession_base: str, moltype:str):
    """
    Given an ensembl accession, return filtered related from Ensembl and NCBI
    Args:
        accession_base (str): An ensembl accession base.

    Returns:
        related (dict): A dictionary containing related sequences from Ensembl, NCBI
    """
    
    # Get taxname and related from ensembl
    taxname, ensembl_related = fetch_ensembl_gene_info(accession_base, moltype)
    # Get related from ncbi using gene symbol and taxname
    gene_symbol = next(iter(ensembl_related))
    _, ncbi_related = _get_related_by_gene_symbol_from_ncbi(
        gene_symbol, taxname
    )
    # Merge and filter related from two sources
    if ensembl_related or ncbi_related:
        related = _merge(ensembl_related, ncbi_related)
        if taxname.upper() == HUMAN_TAXON:
            return filter_related(accession_base, related)
        return related


def _get_related_by_ncbi_id(accession, moltype, locations):
    related = {}
    accession_base = accession.split(".")[0]

    # Try to get taxname and related from ncbi datasets
    if moltype in [MoleculeType.RNA, MoleculeType.PROTEIN]:
        taxname, ncbi_related = _get_related_by_accession_from_ncbi(
            accession
        )
    elif moltype == MoleculeType.DNA and "NC_" in accession:
        taxname, ncbi_related = _get_related_by_chr_location(
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
            if provider.get("name") == "ENSEMBL"
        ]

        for ensembl_gene in ensembl_genes_id:
            taxname, ensembl_gene_related = fetch_ensembl_gene_info(ensembl_gene, moltype="dna")
            if ensembl_gene_related:
                related = _merge(ensembl_gene_related, ncbi_related)
        if taxname and taxname.upper() == "HOMO SAPIENS":
            related = filter_related(accession_base, related)

    return related


def detect_sequence_source(seq_id: str) -> tuple[str, str]:
    """
    Detects the source and molecular type of a sequence ID.
    Args:
        seq_id (str): The sequence ID string to evaluate.
    Returns: (source: str, moltype: str):
                source: "ensembl", "ncbi", or "other"
                moltype: "dna", "rna", "protein", "unknown" and others
    """
    seq_id = seq_id.strip()

    # Ensembl declares its identifiers should be in the form of
    # ENS[species prefix][feature type prefix][a unique eleven digit number]
    # See at https://www.ensembl.org/info/genome/stable_ids/index.html
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
        moltype = ensembl_feature_map.get(prefix, "unknown")
        return "ensembl", moltype

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
        return "ncbi", moltype

    return "other", "unknown"


def get_related(accession, locations=None, timeout=DEFAULT_TIMEOUT):
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
    accession_base = accession.split(".")[0]

    source, moltype = detect_sequence_source(accession)
    if source == "ensembl" and moltype != "unkown":
        return _get_related_by_ensembl_id(accession_base, moltype)
    elif source == "ncbi" and moltype != "unknown":
        return _get_related_by_ncbi_id(
            accession, moltype, locations=locations
        )
    elif source == "other":
        return _get_related_by_gene_symbol(accession)
    else:
        raise NameError(f"Could not retrieve related for {accession}.")
