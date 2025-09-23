import json
import re
from urllib.parse import quote
import requests
from mutalyzer_retriever.configuration import cache_url, settings
from mutalyzer_retriever.request import Http400, request
from mutalyzer_retriever.reference import GRCH37


DEFAULT_TIMEOUT = 10


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


def _valid_locations(accession:str, locations: str):
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



def _merge_assemblies(ebi_assemblies: list, ncbi_assemblies: list):
    "Merge two lists of assemblies gathered from ensembl and ncbi"
    return ebi_assemblies + ncbi_assemblies


def _merge_transcripts(ebi_transcripts: list, ncbi_transcripts: list):
    "Merge two lists of transcripts gathered from ensembl and ncbi"
    for ebi_t in ebi_transcripts:
        matched = False
        for transcript in ncbi_transcripts:
            for i, provider in enumerate(transcript.get("providers", [])):
                if provider.get("name") == "ENSEMBL" and provider.get(
                    "transcript_accession"
                ) == ebi_t.get("transcript_accession"):
                    transcript["providers"][i] = ebi_t
                    matched = True
        if not matched:
            ncbi_transcripts.append({"providers": [clean_dict(ebi_t)]})
    return ncbi_transcripts


def _merge_gene(ebi_gene, ncbi_genes):
    "Merge two lists of related genes gathered from ensembl and ncbi"
    genes = []
    for gene in ncbi_genes:
        gene_name = gene.get("name", "")
        if gene_name == next(g for g in ebi_gene if g != "providers"):
            ebi_gene_provider = ebi_gene.get("providers")
            if ebi_gene_provider and all(
                ebi_gene_provider != p for p in gene.get("providers", [])
            ):
                gene["providers"].append(ebi_gene_provider)
        gene["transcripts"] = _merge_transcripts(
            ebi_gene.get(gene_name, []), gene.get("transcripts", [])
        )
        genes.append(gene)
    return genes


def _merge(ebi_related, ncbi_related):
    # Merge gene sets related from ENSEMBL and NCBI
    related = {
        "assemblies": _merge_assemblies(
            ebi_related.get("assemblies", []),
            ncbi_related.get("assemblies", []),
        ),
        "genes": _merge_gene(ebi_related, ncbi_related.get("genes", [])),
    }
    return related


def _get_grch37_chr_accession(chrid):
    # Use info from reference
    return GRCH37.get(chrid)


def _parse_assemblies(dataset_report):
    taxname = None
    assemblies = []

    gene = dataset_report.get("gene", {})
    taxname = taxname or gene.get("taxname", "Unknown")
    sequence_name = None

    for annotation in gene.get("annotations", []):
        for loc in annotation.get("genomic_locations", []):
            accession = loc.get("genomic_accession_version")
            assembly_name = annotation.get("assembly_name")
            sequence_name = loc.get("sequence_name")

            assembly_entry = clean_dict(
                {
                    "name": assembly_name,
                    "accession": accession,
                }
            )

            if assembly_entry and assembly_entry not in assemblies:
                assemblies.append(assembly_entry)

    # Add GRCh37 assembly if human
    if taxname == "Homo sapiens" and sequence_name:
        grch37_acc = _get_grch37_chr_accession(
            sequence_name
        )
        grch37_entry = clean_dict(
            {
                "name": "GRCh37.p13",
                "accession": grch37_acc,
            }
        )
        if grch37_entry and grch37_entry not in assemblies:
            assemblies.append(grch37_entry)
    return assemblies


def _parse_genes(dataset_report):

    gene = dataset_report.get("gene", {})
    symbol = gene.get("symbol")
    description = gene.get("description")
    ensembl_ids = gene.get("ensembl_gene_ids")

    # Reference and provider info
    ncbi_id = gene.get("gene_id")
    hgnc_id_raw = gene.get("nomenclature_authority", {}).get("identifier")
    hgnc_id = hgnc_id_raw.split(":")[1] if hgnc_id_raw else None

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
    genomic_related = None
    product_related = None

    base_api = settings.get("NCBI_DATASETS_API")
    taxname_url_str = quote(taxname, safe="")
    dataset_report_url = (
        f"{base_api}/gene/symbol/{gene_symbol}/taxon/"
        f"{taxname_url_str}/dataset_report"
    )

    product_report_url = (
        f"{base_api}/gene/symbol/{gene_symbol}/taxon/"
        f"{taxname_url_str}/product_report"
    )


    taxname, genomic_related = _get_genomic_related_from_datasets(
        dataset_report_url
    )
    taxname, product_related = _get_product_related_from_datasets(
        product_report_url
    )

    if product_related or genomic_related:
        related = _merge_datasets(genomic_related, product_related)
        return taxname, related
    return None, None


def _get_related_by_gene_symbol_from_ensembl(
    gene_symbol, taxname="homo Sapiens"
):
    if not gene_symbol:
        return None
    url = (
        f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}"
        f"?content-type=application/json;expand=1"
    )
    try:
        ensembl_gene_lookup_json = json.loads(
            request(url=url, timeout=DEFAULT_TIMEOUT)
        )
        parsed = _parse_ebi_lookup_json(ensembl_gene_lookup_json)
        if parsed:
            taxname, ensembl_id, gene, ebi_related_transcripts = parsed
            return taxname, {
                gene: ebi_related_transcripts,
                "providers": {"name": "ENSEMBL", "accession": ensembl_id},
            }
    except Http400 as e:
        if e.status_code == 400:
            raise ValueError(
                f"Cannot fetch: invalid gene symbol '{gene_symbol}'"
            ) from e


def _get_related_by_gene_symbol(gene_symbol):
    _, ncbi_related = _get_related_by_gene_symbol_from_ncbi(
        gene_symbol, taxname="Homo Sapiens"
    )
    _, ebi_related = _get_related_by_gene_symbol_from_ensembl(
        gene_symbol, taxname="Homo Sapiens"
    )
    related = _merge(ebi_related, ncbi_related)
    related = filter_related(gene_symbol, related)
    return related


def _fetch_ebi_lookup_grch38(accession_base: str, expand=1):
    if not accession_base:
        return
    url = (
        f"https://rest.ensembl.org/lookup/id/{accession_base}"
        f"?content-type=application/json;expand={expand}"
    )

    return json.loads(request(url=url, timeout=DEFAULT_TIMEOUT))


def _parse_ebi_gene_lookup_json(ebi_gene_json):
    ebi_related_transcripts = []
    gene_symbol = ebi_gene_json.get("display_name")
    gene_id = ebi_gene_json.get("id")
    taxname = ebi_gene_json.get("species").replace("_", " ")

    for transcript in ebi_gene_json.get("Transcript", []):
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

        ebi_related_transcripts.append(
            {
                "name": "ENSEMBL",
                "transcript_accession": f"{transcript_id}.{transcript_version}",
                "protein_accession": protein_accession,
                "description": transcript.get("display_name"),
            }
        )

    return taxname, gene_id, gene_symbol, ebi_related_transcripts


def _parse_ebi_transcript_lookup_json(ebi_transcript_json):
    ebi_gene_json = _fetch_ebi_lookup_grch38(ebi_transcript_json.get("Parent"))
    (
        taxname,
        gene_id,
        gene_symbol,
        ebi_related_transcripts,
    ) = _parse_ebi_lookup_json(ebi_gene_json)
    return taxname, gene_id, gene_symbol, ebi_related_transcripts


def _parse_ebi_protein_lookup_json(ebi_protein_json):
    ebi_transcript_json = _fetch_ebi_lookup_grch38(
        ebi_protein_json.get("Parent")
    )
    (
        taxname,
        gene_id,
        gene_symbol,
        ebi_related_transcripts,
    ) = _parse_ebi_transcript_lookup_json(ebi_transcript_json)
    return taxname, gene_id, gene_symbol, ebi_related_transcripts


def _parse_ebi_lookup_json(ebi_json):
    if ebi_json.get("object_type") == "Gene" and ebi_json.get("Transcript"):
        return _parse_ebi_gene_lookup_json(ebi_json)
    if ebi_json.get("object_type") == "Transcript" and ebi_json.get("Parent"):
        return _parse_ebi_transcript_lookup_json(ebi_json)
    if ebi_json.get("object_type") == "Translation" and ebi_json.get("Parent"):
        return _parse_ebi_protein_lookup_json(ebi_json)


def fetch_ensembl_gene_info(accession_base: str, moltype: str) -> tuple[str, dict]:
    """
    Fetch gene-related information for a given Ensembl accession.

    Args:
        accession_base (str): Ensembl accession, support for gene, transcript, protein.
        moltype (str): Molecule type — one of "dna", "rna", "protein" or other

    Returns:
        tuple: (taxname, related_info_dict)
            - taxname (str): Taxonomic name (e.g., "Homo sapiens")
            - related_info_dict (dict): Dictionary including gene and its products.

    Raises:
        ValueError: If the Ensembl lookup fails or returns nothing.
    """
    if moltype not in {"dna", "rna", "protein"}:
        raise ValueError(f"Unsupported molecule type: {moltype}")

    try:
        # Fetch data from EBI for this accession
        ebi_lookup_json = _fetch_ebi_lookup_grch38(
            accession_base,
            expand=0 if moltype == "protein" else 1
        )

        if not ebi_lookup_json:
            raise ValueError(f"No data returned for {accession_base} from ENSEMBL.")

        parsed = _parse_ebi_lookup_json(ebi_lookup_json)
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

    except Exception as e:
        raise RuntimeError(f"Error fetching data for {accession_base} from ENSEMBL: {e}")



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


def _fetch_related_from_ncbi_dataset_report(gene_ids):
    base_api = settings.get("NCBI_DATASETS_API")
    gene_id_str = quote(",".join(map(str, gene_ids)))
    url = f"{base_api}/gene/id/{gene_id_str}/dataset_report"
    genes_dataset_report_json = json.loads(
        request(url=url, timeout=DEFAULT_TIMEOUT)
    )
    taxname, dataset_related = _parse_dataset_report(genes_dataset_report_json)
    return taxname, dataset_related


def _fetch_related_from_ncbi_product_report(gene_ids):
    base_api = settings.get("NCBI_DATASETS_API")
    gene_id_str = quote(",".join(map(str, gene_ids)))
    url = f"{base_api}/gene/id/{gene_id_str}/product_report"
    genes_product_report_json = json.loads(
        request(url=url, timeout=DEFAULT_TIMEOUT)
    )
    if genes_product_report_json:
        return _parse_product_report(genes_product_report_json)


def _get_gene_related(gene_ids):
    _, genomic_related = _fetch_related_from_ncbi_dataset_report(
        gene_ids
    )
    _, product_related = _fetch_related_from_ncbi_product_report(
        gene_ids
    )
    return genomic_related, product_related


def _parse_genome_annotation_report(genome_annotation_report):
    gene_ids = []
    for report in genome_annotation_report.get("reports", []):
        annotation = report.get("annotation", {})
        gene_id = annotation.get("gene_id")
        taxname = annotation.get("taxname")
        if gene_id is not None:
            gene_ids.append(gene_id)
    if gene_ids:
        genomic_related, product_related = _get_gene_related(gene_ids)
        related = _merge_datasets(genomic_related, product_related)
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

    api_base = settings.get("NCBI_DATASETS_API")
    endpoint = f"genome/accession/{assembly_accession}/annotation_report"
    params = [("locations", loc) for loc in _valid_locations(accession, locations)]
    genome_response = json.loads(
        request(
            url=f"{api_base}/{endpoint}",
            params=params,
            timeout=DEFAULT_TIMEOUT,
        )
    )
    taxname, related = _parse_genome_annotation_report(genome_response)

    return taxname, related


def _get_assembly_accession(accession):
    api_base = settings.get("NCBI_DATASETS_API")
    endpoint = f"genome/sequence_accession/{accession}/sequence_assemblies"
    assembly_json = json.loads(
        request(url=f"{api_base}/{endpoint}", timeout=DEFAULT_TIMEOUT)
    )
    accessions = assembly_json.get("accessions")
    if isinstance(accessions, list) and accessions:
        return accessions[0]
    return None


def _get_genomic_related_from_datasets(dataset_report_url):
    try:
        dataset_json = json.loads(
            request(url=dataset_report_url, timeout=DEFAULT_TIMEOUT)
        )
        if dataset_json:
            taxname, genomic_related = _parse_dataset_report(dataset_json)
            return taxname, genomic_related
    except Exception as e:
        raise RuntimeError(f"Error fetching related accessions: {e}") from e
    return None, None


def _get_product_related_from_datasets(product_report_url):
    try:
        product_json = json.loads(request(url=product_report_url, timeout=DEFAULT_TIMEOUT))
        if product_json:
            taxname, product_related = _parse_product_report(product_json)
            return taxname, product_related
    except Exception as e:
        raise RuntimeError(f"Error fetching related accessions: {e}") from e
    return None, None


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
    genomic_related = None
    product_related = None
    taxname = None

    base_url = settings.get("NCBI_DATASETS_API")

    accession_without_versions = accession.split(".")[0]
    dataset_report_url = (
        f"{base_url}/gene/accession/{accession_without_versions}/dataset_report"
    )
    product_related_url = (
        f"{base_url}/gene/accession/{accession_without_versions}/product_report"
    )
    taxname, genomic_related = _get_genomic_related_from_datasets(
        dataset_report_url
    )
    taxname, product_related = _get_product_related_from_datasets(
        product_related_url
    )
    if genomic_related or product_related:
        related = _merge_datasets(genomic_related, product_related)
        return taxname, related
    return None, None


def _get_related_by_ensembl_id(accession_base: str, moltype:str):
    """
    Given an ensembl accession, return filtered related from Ensembl and NCBI
    Args:
        accession_base (str): An ensembl accession base.

    Returns:
        related (dict): A dictionary containing related sequences from Ensembl, NCBI
    """
    related = {}
    # Get taxname and related from ensembl
    taxname, ebi_related = fetch_ensembl_gene_info(accession_base, moltype)
    # Get related from ncbi using gene symbol and taxname
    gene_symbol = next(iter(ebi_related))
    _, ncbi_related = _get_related_by_gene_symbol_from_ncbi(
        gene_symbol, taxname
    )
    # Merge and filter related from two sources
    if ebi_related or ncbi_related:
        related = _merge(ebi_related, ncbi_related)
        if taxname.upper() == "HOMO SAPIENS":
            related = filter_related(accession_base, related)
    return related


def _get_related_by_ncbi_id(accession, moltype, locations):
    related = {}
    accession_base = accession.split(".")[0]

    # Try to get taxname and related from ncbi datasets
    try:
        if moltype.upper() in ["RNA", "PROTEIN"]:
            taxname, ncbi_related = _get_related_by_accession_from_ncbi(
                accession
            )
        elif moltype.upper() == "DNA" and "NC_" in accession:
            taxname, ncbi_related = _get_related_by_chr_location(
                accession, locations
            )
        else:
            raise NameError(f"Could not retrieve related for {accession}.")
    except Exception as e:
        raise RuntimeError(f"Error fetching related accessions: {e}") from e
    # Get related from ensembl using ensembl gene id.
    if ncbi_related:
        ensembl_genes_id = [
            provider["accession"]
            for gene in ncbi_related.get("genes", [])
            for provider in gene.get("providers", [])
            if provider.get("name") == "ENSEMBL"
        ]

        for ensembl_gene in ensembl_genes_id:
            taxname, ebi_gene_related = fetch_ensembl_gene_info(ensembl_gene, moltype="dna")
            if ebi_gene_related:
                related = _merge(ebi_gene_related, ncbi_related)
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
                moltype: "dna", "rna", "protein" or "unknown"
    """
    seq_id = seq_id.strip()

    # Ensembl declares its identifiers should be in the form of
    # ENS[species prefix][feature type prefix][a unique eleven digit number]
    # See at https://www.ensembl.org/info/genome/stable_ids/index.html
    ensembl_feature_map = {
        "E": "exon",
        "FM": "Ensembl protein family",
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
    related = {}

    source, moltype = detect_sequence_source(accession)
    if source == "ensembl" and moltype != "unkown":
        related = _get_related_by_ensembl_id(accession_base, moltype)
    elif source == "ncbi" and moltype != "unknown":
        related = _get_related_by_ncbi_id(
            accession, moltype, locations=locations
        )
    elif source == "other":
        related = _get_related_by_gene_symbol(accession)
    else:
        raise NameError(f"Could not retrieve related for {accession}.")

    return related
