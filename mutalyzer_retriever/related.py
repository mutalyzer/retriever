# A script to get related references from ncbi and ebi
import argparse
import ast
import json
import re
from urllib.parse import quote
from mutalyzer_retriever.request import Http400, request


def ncbi_urls(endpoint: str):
    urls = {
        "ELink": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi",
        "ESummary": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        "EFetch": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        "Datasets_gene": "https://api.ncbi.nlm.nih.gov/datasets/v2/gene",
        "Datasets_genome": "https://api.ncbi.nlm.nih.gov/datasets/v2/genome",
    }
    return urls.get(endpoint)


def clean_dict(d):
    """Remove keys with None, empty string, or empty list values."""
    return {k: v for k, v in d.items() if v not in (None, "", [])}


def _merge_assemblies(ebi_assemblies, ncbi_assemblies):
    assemblies = []
    assemblies.extend(ebi_assemblies)
    assemblies.extend(ncbi_assemblies)
    return list({frozenset(a.items()): a for a in assemblies}.values())


def _merge_transcripts(ebi_transcripts, ncbi_transcripts):
    for ebi_t in ebi_transcripts:
        matched = False
        for transcript in ncbi_transcripts:
            for i, provider in enumerate(transcript.get("providers", [])):
                if (
                    provider.get("name") == "ENSEMBL" and
                    provider.get("transcript_accession") == ebi_t.get("transcript_accession")
                ):
                    transcript["providers"][i] = ebi_t
                    matched = True
        if not matched:
            ncbi_transcripts.append({"providers": [clean_dict(ebi_t)]})
    return ncbi_transcripts


def _merge_gene(ebi_gene, ncbi_genes):
    genes = []
    for gene in ncbi_genes:
        gene_name = gene.get("name", "unknown")
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
    # Merge gene sets related from ENSEMBL into NCBI
    related = {
        "assemblies": _merge_assemblies(
            ebi_related.get("assemblies", []),
            ncbi_related.get("assemblies", []),
        ),
        "genes": _merge_gene(
                    ebi_related, ncbi_related.get("genes", [])
                    )
    }
    return related


def _get_grch37_chr_accession(chrid, timeout=10):
    url = (
        f"{ncbi_urls('Datasets_genome')}/accession/GCF_000001405.25/"
        f"sequence_reports?chromosomes={chrid}&"
        f"role_filters=assembled-molecule"
    )
    response = json.loads(request(url=url, timeout=timeout))
    return response.get("reports", [{}])[0].get("refseq_accession")


def _parse_dataset_report(json_report):
    """
    Parse NCBI gene dataset report.
    Args:
        json_report (dict): JSON data from NCBI dataset API.
    Returns:
        tuple: 
            - taxname (str)
            - dict: {
                "assemblies": list of assembly dicts,
                "genes": list of gene dicts
              }
    """
    assemblies = []
    genes = []
    taxname = None

    for report in json_report.get("reports", []):
        gene = report.get("gene", {})
        symbol = gene.get("symbol")
        description = gene.get("description")
        ensembl_ids = gene.get("ensembl_gene_ids")
        taxname = taxname or gene.get("taxname", "Unknown")
        sequence_name = None

        # Extract genomic assemblies
        for annotation in gene.get("annotations", []):
            for loc in annotation.get("genomic_locations", []):
                accession = loc.get("genomic_accession_version")
                sequence_name = loc.get("sequence_name")
                assembly_name = annotation.get("assembly_name")

                assembly_entry = clean_dict({
                    "assembly_name": assembly_name,
                    "accession": accession,
                })

                if assembly_entry and assembly_entry not in assemblies:
                    assemblies.append(assembly_entry)

        # Add GRCh37 assembly if human
        if taxname == "Homo sapiens" and sequence_name:
            grch37_acc = _get_grch37_chr_accession(sequence_name)
            grch37_entry = clean_dict({
                "assembly_name": "GRCh37.p13",
                "accession": grch37_acc,
            })

            if grch37_entry and grch37_entry not in assemblies:
                assemblies.append(grch37_entry)

        # Reference and provider info
        ncbi_id = gene.get("gene_id")
        hgnc_id_raw = gene.get("nomenclature_authority", {}).get("identifier")
        hgnc_id = hgnc_id_raw.split(":")[1] if hgnc_id_raw else None

        ref_standards = gene.get("reference_standards", [])
        ref_accession = None
        for ref in ref_standards:
            ref_accession = ref.get("gene_range", {}).get("accession_version")

        ncbi = clean_dict({
            "name": "NCBI",
            "id": ncbi_id,
            "accession": ref_accession,
        })

        ensembl = clean_dict({
            "name": "ENSEMBL",
            "accession": ensembl_ids[0] if ensembl_ids else None,
        })

        providers = [p for p in (ncbi, ensembl) if len(p) > 1]
        if not providers:
            continue

        gene_entry = clean_dict({
            "name": symbol,
            "hgnc_id": hgnc_id,
            "description": description,
            "providers": providers,
        })

        if gene_entry:
            genes.append(gene_entry)

    return taxname, {"assemblies": assemblies, "genes": genes}



def _parse_product_report(data):
    """
    Parse NCBI product report JSON into a dictionary.
    Args:
        data (dict): JSON data from ncbi product API.
    Returns:
        dict: Mapping gene_symbol -> list of transcript info dicts
    """
    genes_dict = {}

    for report in data.get("reports", []):
        product = report.get("product", {})
        symbol = product.get("symbol", "UNKNOWN")
        gene_products = []

        for transcript in product.get("transcripts", []):
            ncbi_transcript = transcript.get("accession_version")
            ncbi_desc = transcript.get("name")
            ncbi_protein = transcript.get("protein", {})
            ncbi_protein_acc = ncbi_protein.get("accession_version")
            ensembl_transcript = transcript.get("ensembl_transcript")
            ensembl_protein = ncbi_protein.get("ensembl_protein")

            ncbi = clean_dict(
                {
                    "name": "NCBI",
                    "transcript_accession": ncbi_transcript,
                    "protein_accession": ncbi_protein_acc,
                    "description": ncbi_desc,
                }
            )

            ensembl = clean_dict(
                {
                    "name": "ENSEMBL",
                    "transcript_accession": ensembl_transcript,
                    "protein_accession": ensembl_protein,
                }
            )

            providers = [p for p in (ncbi, ensembl) if len(p) > 1]
            if not providers:
                continue

            product_entry = clean_dict(
                {
                    "providers": providers,
                    "tag": transcript.get("select_category"),
                }
            )
            gene_products.append(product_entry)

        if gene_products:
            genes_dict[symbol] = gene_products

    return genes_dict


def _merge_genome_products(genomic_related, product_related):
    """
    Merges genomic and product-related gene data.
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
                    "refseqgene": gene.get("refseqgene"),
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
        return
    for report in reports.get("reports"):
        for key, value in report.items():
            if value.get("symbol") == gene_symbol.upper():
                return {"reports": [{key: value}]}


def _get_related_by_gene_symbol_from_ncbi(
    gene_symbol, taxname="Homo Sapiens", timeout=10
):
    """
    Given a gene symbol, return a set of related sequence accessions (genomic and/or products).
    Returns (taxname, related_dict), or (None, None) if nothing found.
    """
    if not gene_symbol:
        return None, None
    related = {}
    genomic_related = None
    product_related = None

    base_url = ncbi_urls("Datasets_gene")
    taxname_url_str = quote(taxname, safe="")
    dataset_report_url = f"{base_url}/symbol/{gene_symbol}/taxon/{taxname_url_str}/dataset_report"

    # Fetch and parse genomic related data
    dataset_json = json.loads(
        request(url=dataset_report_url, timeout=timeout)
    )
    dataset_json = filter_report_from_other_genes(gene_symbol, dataset_json)
    _, genomic_related = (
        _parse_dataset_report(dataset_json)
        if dataset_json
        else (None, None)
    )

    # Fetch and parse product related data
    product_url = (
        f"{base_url}/symbol/{gene_symbol}/taxon/9606/product_report"
    )
    product_json = json.loads(request(url=product_url, timeout=timeout))
    product_json = filter_report_from_other_genes(gene_symbol, product_json)
    product_related = (
        _parse_product_report(product_json)
        if product_json
        else (None, None)
    )

    if product_related or genomic_related:
        related = _merge_genome_products(genomic_related, product_related)
        return "Homo sapiens", related
    return None, None


def _get_related_by_gene_symbol_ensembl(
    gene_symbol, taxname="homo Sapiens", timeout=10
):
    if not gene_symbol:
        return
    url = (
        f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}"
        f"?content-type=application/json;expand=1"
    )
    try:
        ensembl_gene_lookup_json = json.loads(request(url=url, timeout=timeout))
        parsed = _parse_ebi_lookup_json(ensembl_gene_lookup_json)
        if parsed:
            taxname, ensembl_id, gene, ebi_related_transcripts = parsed
            return taxname, {
                gene: ebi_related_transcripts,
                "providers": {"name": "ENSEMBL", "accession": ensembl_id},
            }        
    except Http400 as e:
        if e.status_code == 400:
            raise ValueError(f"Cannot fetch: invalid gene symbol '{gene_symbol}'") from e 



def _get_related_by_gene_symbol(gene_symbol, timeout=10):
    _, ncbi_related = _get_related_by_gene_symbol_from_ncbi(
        gene_symbol, taxname="Homo Sapiens", timeout=10
    )
    _, ebi_related = _get_related_by_gene_symbol_ensembl(
        gene_symbol, taxname="Homo Sapiens", timeout=10
    )
    related = _merge(ebi_related, ncbi_related)
    related = filter_related(gene_symbol, related)
    return related


def _fetch_ebi_lookup_grch38(accession_base: str, expand=1, timeout=10):
    if not accession_base:
        return
    url = (
        f"https://rest.ensembl.org/lookup/id/{accession_base}"
        f"?content-type=application/json;expand={expand}"
    )

    return json.loads(request(url=url, timeout=timeout))


def _parse_ebi_gene_lookup_json(ebi_gene_json):
    ebi_related_transcripts = []
    gene_symbol = ebi_gene_json.get("display_name")
    gene_id = ebi_gene_json.get('id')
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
    ebi_gene_json = _fetch_ebi_lookup_grch38(
        ebi_transcript_json.get("Parent")
    )
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
    if ebi_json.get("object_type") == "Transcript" and ebi_json.get(
        "Parent"
    ):
        return _parse_ebi_transcript_lookup_json(ebi_json)
    if ebi_json.get("object_type") == "Translation" and ebi_json.get(
        "Parent"
    ):
        return _parse_ebi_protein_lookup_json(ebi_json)


def fetch_ensembl_gene_info(accession_base: str):
    gene = None
    ebi_related_transcripts = []
    taxname = None

    try:
        ebi_lookup_json = _fetch_ebi_lookup_grch38(accession_base)
    except Http400 as http_e:
        try:
            error_json = http_e.response.json()
            if (
                error_json.get("error")
                == "Expand option only available for Genes and Transcripts"
            ):
                ebi_lookup_json = _fetch_ebi_lookup_grch38(
                    accession_base, expand=0
                )
                parsed = _parse_ebi_lookup_json(ebi_lookup_json)
                if parsed:
                    (
                        taxname,
                        ensembl_id,
                        gene,
                        ebi_related_transcripts,
                    ) = parsed
                    return taxname, {gene: ebi_related_transcripts}
            else:
                raise
        except Exception as parse_e:
            raise RuntimeError(
                f"Failed to parse error : {parse_e}"
            ) from parse_e

    if ebi_lookup_json:
        parsed = _parse_ebi_lookup_json(ebi_lookup_json)
        if parsed:
            taxname, ensembl_id, gene, ebi_related_transcripts = parsed
        else:
            raise ValueError(
                f"Failed to retrieve related data from ENSEMBL for {accession_base}"
            )
    return taxname, {gene: ebi_related_transcripts}


def filter_assemblies(related_assemblies:list):
    filtered_assemblies = []

    for assembly in related_assemblies:
        if assembly.get("accession", "Unknown").startswith("NC_"):
            filtered_assemblies.append(clean_dict(assembly))
    return filtered_assemblies


def filter_gene_products(ID_base, related_genes:list):
    filtered_genes = []

    for gene in related_genes:
        filtered_transcripts = []
        for transcript in gene.get("transcripts", []):
            matching_providers = [
                provider
                for provider in transcript.get("providers", [])
                if ID_base
                == (provider.get("transcript_accession") or "").split(".")[
                    0
                ]
                or ID_base
                == (provider.get("protein_accession") or "").split(".")[0]
            ]
            filtered_transcript = {}
            if "tag" in transcript or matching_providers:
                filtered_transcript["providers"] = transcript.get(
                    "providers", []
                )
                filtered_transcripts.append(transcript)

        gene["transcripts"] = filtered_transcripts
        filtered_genes.append(clean_dict(gene))
    return filtered_genes


def filter_related(ID_base, related):
    return {
        "assemblies": filter_assemblies(related.get("assemblies", [])),
        "genes": filter_gene_products(ID_base, related.get("genes", []))
    }


def _fetch_related_from_ncbi_dataset_report(gene_ids, timeout):
    base_url = ncbi_urls("Datasets_gene")
    gene_id_str = quote(",".join(map(str, gene_ids)))
    url = f"{base_url}/id/{gene_id_str}/dataset_report"
    genes_dataset_report_json = json.loads(
        request(url=url, timeout=timeout)
    )
    taxname, dataset_related = _parse_dataset_report(
        genes_dataset_report_json
    )
    return taxname, dataset_related


def _fetch_related_from_ncbi_product_report(gene_ids, taxname, timeout):
    base_url = ncbi_urls("Datasets_gene")
    gene_id_str = quote(",".join(map(str, gene_ids)))
    url = f"{base_url}/id/{gene_id_str}/product_report"
    genes_product_report_json = json.loads(
        request(url=url, timeout=timeout)
    )
    if genes_product_report_json:
        return _parse_product_report(genes_product_report_json)


def _get_gene_related(gene_ids, timeout):
    taxname, genomic_related = _fetch_related_from_ncbi_dataset_report(
        gene_ids, timeout
    )

    product_related = _fetch_related_from_ncbi_product_report(
        gene_ids, taxname, timeout
    )
    return genomic_related, product_related


def _get_related_by_chr_location(accession, locations, timeout):
    gene_ids = []
    related = {}
    assembly_accession = _get_assembly_accession(accession, timeout=timeout)
    if not assembly_accession:
        raise ValueError(
            f"Assembly accession could not be determined for {accession}"
        )
    location_params = [
        f"{accession}:{start}-{end}" for start, end in locations
    ]
    encoded_locations = [quote(loc) for loc in location_params]
    location_string = "&".join(
        [f"locations={loc}" for loc in encoded_locations]
    )
    url = (
        f"{ncbi_urls('Datasets_genome')}/accession/{assembly_accession}/"
        f"annotation_report?{location_string}"
    )

    taxname = None
    genome_response = json.loads(request(url=url, timeout=timeout))
    for report in genome_response.get("reports", []):
        annotation = report.get("annotation", {})
        gene_id = annotation.get("gene_id")
        taxname = annotation.get("taxname")
        if gene_id is not None:
            gene_ids.append(gene_id)
    if gene_ids:
        genomic_related, product_related = _get_gene_related(
            gene_ids, timeout=timeout
        )
        related = _merge_genome_products(genomic_related, product_related)
    return taxname, related


def _get_assembly_accession(accession, timeout):
    url = f"{ncbi_urls('Datasets_genome')}/sequence_accession/{accession}/sequence_assemblies"
    assembly_json = json.loads(request(url=url, timeout=timeout))
    accessions = assembly_json.get("accessions")
    if isinstance(accessions, list) and accessions:
        return accessions[0]
    return None


def _get_related_by_accession_from_ncbi(accessions, timeout):
    """
    Input: A sequence (transcripts or proteins) accession.
    Output: A set of related sequence accessions.
    """
    related = {}
    genomic_related = None
    product_related = None
    taxname = None

    base_url = ncbi_urls("Datasets_gene")
    accessions_without_versions = accessions.split(".")[0]
    dataset_report_url = (
        f"{base_url}/accession/{accessions_without_versions}/dataset_report"
    )

    dataset_json = json.loads(
        request(url=dataset_report_url, timeout=timeout)
    )
    if dataset_json:
        taxname, genomic_related = _parse_dataset_report(dataset_json)

    product_related_url = (
        f"{base_url}/accession/{accessions_without_versions}/product_report"
    )
    product_json = json.loads(
        request(url=product_related_url, timeout=timeout)
    )
    if product_json:
        product_related = _parse_product_report(product_json)
    if genomic_related or product_related:
        related = _merge_genome_products(genomic_related, product_related)
        return taxname, related
    return None, None


def _get_related_by_ensembl_id(accession_base):
    taxname, ebi_related = fetch_ensembl_gene_info(accession_base)
    gene_symbol = next(iter(ebi_related))
    _, ncbi_related = _get_related_by_gene_symbol_from_ncbi(
        gene_symbol, taxname
    )
    if ebi_related or ncbi_related:
        related = _merge(ebi_related, ncbi_related)
        if taxname.upper() == "HOMO SAPIENS":
            related = filter_related(accession_base, related)
        return related


def _get_related_by_ncbi_id(
    accession, moltype, locations=[0, 0], timeout=30
):
    accession_base = accession.split(".")[0]
    try:
        if moltype.upper() in ["RNA", "PROTEIN"]:
            taxname, ncbi_related = _get_related_by_accession_from_ncbi(
                accession, timeout=timeout
            )
        elif moltype.upper() == "DNA" and "NC_" in accession:
            taxname, ncbi_related = _get_related_by_chr_location(
                accession, locations, timeout
            )
        else:
            raise NameError(f"Could not retrieve related for {accession}.")
    except Exception as e:
        raise RuntimeError(f"Error fetching related accessions: {e}") from e
    if ncbi_related:
        ebi_related = {}
        ensembl_genes_id = [
            provider["accession"]
            for gene in ncbi_related.get("genes", [])
            for provider in gene.get("providers", [])
            if provider.get("name") == "ENSEMBL"
        ]

        for ensembl_gene in ensembl_genes_id:
            taxname, ebi_gene_related = fetch_ensembl_gene_info(
                ensembl_gene
            )
            if ebi_gene_related:
                related = _merge(ebi_gene_related, ncbi_related)
        if taxname and taxname.upper() == "HOMO SAPIENS":
            related = filter_related(accession_base, related)

        return related


def detect_sequence_source(seq_id: str):
    """
    Detects if the input string is an Ensembl ID or a RefSeq accession or .
    Input: seq_id (str): The sequence ID string to evaluate.
    Returns: tuple: (source: str, moltype: str)
            source: "ensembl", "ncbi", or "other"
            moltype: "dna", "rna", "protein" or "unknown"
    """
    seq_id = seq_id.strip()

    # Ensembl declares its identifiers should be in the form of
    # ENS[species prefix][feature type prefix][a unique eleven digit number]
    # See at https://www.ensembl.org/info/genome/stable_ids/index.html
    if re.match(r"^ENS[A-Z0-9]+(?:\.\d+)?$", seq_id):
        return "ensembl", "unknown"

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

    match = re.match(r"^([A-Z]{2})_\d+(\.\d+)?$", seq_id)
    if match:
        prefix = match.group(1)
        moltype = refseq_moltype_map.get(prefix, "unknown")
        return "ncbi", moltype

    return "other", "unknown"


def get_new_related(accession, locations=[[0, 0]], timeout=10):
    """
    Retrieve related assembly/gene/transcript/protein information based
    on accession or gene symbol
    Args:
        accession (str): A sequence accession (e.g., RefSeq, Ensembl ID) or a human gene symbol.
        locations (list of list[int], optional): Genomic location ranges in the form [[start, end], ...].
            Defaults to [[0, 0]].
        timeout (int, optional): Timeout in seconds for remote API calls. Defaults to 10.

    Returns:
        dict: A dictionary containing related information retrieved from Ensembl, NCBI,

    Raises:
        NameError: If the given accession is not from NCBI RefSeq or ENSEMBL.
    """
    accession = accession.upper()
    accession_base = accession.split(".")[0]
    related = {}

    source, moltype = detect_sequence_source(accession)

    if source == "ensembl":
        related = _get_related_by_ensembl_id(accession_base)
    elif source == "ncbi" and moltype != "unknown":
        related = _get_related_by_ncbi_id(
            accession, moltype, locations=locations
        )
    elif source == "other":
        related = _get_related_by_gene_symbol(accession)
    else:
        raise NameError(f"Could not retrieve related for {accession}.")

    return related


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get related sequences")
    parser.add_argument(
        "accession",
        help="A sequence accession, support gene name, ncbi or ebi accession.",
    )
    parser.add_argument(
        "locations",
        nargs="?",
        type=ast.literal_eval,
        help="A list of locations, e.g. [[112088000, 112088000]]",
    )

    args = parser.parse_args()

    print(
        json.dumps(
            get_new_related(args.accession, args.locations), indent=2
        )
    )
