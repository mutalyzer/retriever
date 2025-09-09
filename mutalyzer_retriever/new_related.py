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


def merge_ncbi_ebi(ebi_related, ncbi_related):
    related = {
        "assemblies": [],
        "genes": [],
    }

    if ncbi_related:
        related = {
            "assemblies": ncbi_related["assemblies"],
            "genes": ncbi_related["genes"],
        }

        if not related["genes"]:
            return related
        # Build a mapping from Ensembl transcript_accession to transcript entry
        for gene in ncbi_related["genes"]:
            gene_symbol = gene.get("name")
            ncbi_transcripts_from_one_gene = gene.get("transcripts")
            transcript_index = {}
            for transcript in ncbi_transcripts_from_one_gene:
                for provider in transcript.get("providers", []):
                    tid = provider.get("transcript_accession")
                    provider_name = provider.get("name")
                    if tid and provider_name == "ENSEMBL":
                        transcript_index[tid] = transcript
            if gene_symbol in ebi_related:
                for ebi_transcripts_from_one_gene in ebi_related[gene_symbol]:
                    ebi_tid = ebi_transcripts_from_one_gene.get("transcript_accession")
                    if not ebi_tid:
                        continue
                    provider_entry = clean_dict(
                        {
                            "name": "ENSEMBL",
                            "transcript_accession": ebi_tid,
                            "protein_accession": ebi_transcripts_from_one_gene.get("protein_accession"),
                            "description": ebi_transcripts_from_one_gene.get("description"),
                        }
                    )

                    if ebi_tid in transcript_index:
                        existing = transcript_index[ebi_tid]
                        providers = existing.get("providers", [])

                        replaced = False
                        for i, p in enumerate(providers):
                            if p.get("name") == "ENSEMBL":
                                providers[i] = provider_entry
                                replaced = True
                                break
                        if not replaced:
                            providers.append(provider_entry)
                    else:
                        new_transcript = clean_dict(
                            {
                                "providers": [provider_entry],
                                "tag": ebi_transcripts_from_one_gene.get("tag"),
                            }
                        )
                        ncbi_transcripts_from_one_gene.append(new_transcript)
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
    Parses ncbi gene dataset report and returns taxname and cleaned data for assemblies and genes.
    Args:
        data (dict): JSON data from ncbi dataset API.

    Returns:
        taxname,
        dict: {"assemblies": list, "genes": list}
    """
    assemblies = []
    genes = []
    taxname = None

    for report in json_report.get("reports", []):
        gene = report.get("gene", {})
        symbol = gene.get("symbol", [])
        ref_accession = None
        gene_description = gene.get("description", [])
        ensembl_gene = gene.get("ensembl_gene_ids")

        # Extract taxname (only once)
        if not taxname:
            taxname = gene.get("taxname", "Unknown")

        sequence_name = None

        for annotation in gene.get("annotations", []):
            for loc in annotation.get("genomic_locations", []):
                genomic_acc = loc.get("genomic_accession_version")
                sequence_name = loc.get("sequence_name")
                assembly_name = annotation.get("assembly_name")

                clean_assembly = clean_dict(
                    {
                        "assembly_name": assembly_name,
                        "accession": genomic_acc,
                    }
                )
                if clean_assembly and clean_assembly not in assemblies:
                    assemblies.append(clean_assembly)

        # Add GRCh37 chromosome accession
        if taxname == "Homo sapiens":
            grch37_acc = _get_grch37_chr_accession(sequence_name)
            grch37_entry = clean_dict(
                {"accession": grch37_acc, "assembly_name": "GRCh37.p13"}
            )

            if grch37_entry and grch37_entry not in assemblies:
                assemblies.append(grch37_entry)

        # Handle reference standard accessions
        ref_acc = gene.get("reference_standards", [])
        hgnc_id_str = gene.get("nomenclature_authority", {}).get("identifier")
        hgnc_id = hgnc_id_str.split(":")[1] if hgnc_id_str else None

        ensembl_gene_id = ensembl_gene[0] if ensembl_gene else None
        ncbi_gene_id = gene.get("gene_id", [])
        for ref in ref_acc:
            ref_accession = ref.get("gene_range", {}).get("accession_version")

        # Build provider dicts, removing empty fields
        ncbi = clean_dict(
            {
                "name": "ncbi",
                "id": ncbi_gene_id,
                "accession": ref_accession,
            }
        )

        ensembl = clean_dict(
            {
                "name": "ENSEMBL",
                "accession": ensembl_gene_id,
            }
        )

        providers = [p for p in (ncbi, ensembl) if len(p) > 1]
        if not providers:
            continue  # skip if both are empty

        gene_entry = clean_dict(
            {
                "name": symbol,
                "hgnc_id": hgnc_id,
                "description": gene_description,
                "providers": providers,
            }
        )
        if gene_entry:
            genes.append(gene_entry)
    return taxname, {"assemblies": assemblies, "genes": genes}


def _parse_product_report(data):
    """
    Parse ncbi product report json into a dict:
    Args:
        data (dict): JSON data from ncbi product API.

    Returns:
        dict: Mapping gene_symbol -> list of transcript info dicts
    """
    genes_dict = {}

    for report in data.get("reports", []):
        product = report.get("product", {})
        transcripts = product.get("transcripts", [])
        symbol = product.get("symbol", "UNKNOWN")
        gene_products = []

        transcripts = product.get("transcripts", [])

        for transcript in transcripts:
            ncbi_transcript_acc = transcript.get("accession_version")
            ncbi_trancript_description = transcript.get("name")
            ncbi_protein = transcript.get("protein", {})
            ncbi_protein_acc = ncbi_protein.get("accession_version")
            ensembl_transcript_acc = transcript.get("ensembl_transcript")
            ensembl_protein_acc = ncbi_protein.get("ensembl_protein")

            # Build provider dicts, removing empty fields
            ncbi = clean_dict(
                {
                    "name": "ncbi",
                    "transcript_accession": ncbi_transcript_acc,
                    "protein_accession": ncbi_protein_acc,
                    "description": ncbi_trancript_description,
                }
            )

            ensembl = clean_dict(
                {
                    "name": "ENSEMBL",
                    "transcript_accession": ensembl_transcript_acc,
                    "protein_accession": ensembl_protein_acc,
                }
            )

            providers = [p for p in (ncbi, ensembl) if len(p) > 1]
            if not providers:
                continue  # skip if both are empty

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


def _merge_related(genomic_related, product_related):
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
            if not transcripts:
                continue

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


def _get_related_by_gene_symbol(gene_symbol, taxname="Homo Sapiens", timeout=10):
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
    dataset_report_url = (
        f"{base_url}/symbol/{gene_symbol}/taxon/{taxname_url_str}/dataset_report"
    )

    # Fetch and parse genomic related data
    dataset_json = json.loads(request(url=dataset_report_url, timeout=timeout))
    _, genomic_related = (
        _parse_dataset_report(dataset_json) if dataset_json else (None, None)
    )

    # Fetch and parse product related data
    product_url = (
        f"{base_url}/symbol/{gene_symbol}/taxon/{taxname_url_str}/product_report"
    )
    product_json = json.loads(request(url=product_url, timeout=timeout))
    product_related = (
        _parse_product_report(product_json) if product_json else (None, None)
    )

    if product_related or genomic_related:
        related = _merge_related(genomic_related, product_related)
        return "Homo sapiens", related
    return None, None


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

    return taxname, gene_symbol, ebi_related_transcripts


def _parse_ebi_transcript_lookup_json(ebi_transcript_json):
    ebi_gene_json = _fetch_ebi_lookup_grch38(ebi_transcript_json.get("Parent"))
    taxname, gene_symbol, ebi_related_transcripts = _parse_ebi_lookup_json(
        ebi_gene_json
    )
    return taxname, gene_symbol, ebi_related_transcripts


def _parse_ebi_protein_lookup_json(ebi_protein_json):
    ebi_transcript_json = _fetch_ebi_lookup_grch38(ebi_protein_json.get("Parent"))
    taxname, gene_symbol, ebi_related_transcripts = _parse_ebi_transcript_lookup_json(
        ebi_transcript_json
    )
    return taxname, gene_symbol, ebi_related_transcripts


def _parse_ebi_lookup_json(ebi_json):
    if ebi_json.get("object_type") == "Gene" and ebi_json.get("Transcript"):
        return _parse_ebi_gene_lookup_json(ebi_json)
    if ebi_json.get("object_type") == "Transcript" and ebi_json.get("Parent"):
        return _parse_ebi_transcript_lookup_json(ebi_json)
    if ebi_json.get("object_type") == "Translation" and ebi_json.get("Parent"):
        return _parse_ebi_protein_lookup_json(ebi_json)


def _get_related_by_accession_from_ebi(accession_base: str):
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
                ebi_lookup_json = _fetch_ebi_lookup_grch38(accession_base, expand=0)
                parsed = _parse_ebi_lookup_json(ebi_lookup_json)
                if parsed:
                    taxname, gene, ebi_related_transcripts = parsed
                    return taxname, {gene: ebi_related_transcripts}
            else:
                raise
        except Exception as parse_e:
            raise RuntimeError(f"Failed to parse error : {parse_e}") from parse_e

    if ebi_lookup_json:
        parsed = _parse_ebi_lookup_json(ebi_lookup_json)
        if parsed:
            taxname, gene, ebi_related_transcripts = parsed
        else:
            raise ValueError(
                f"Failed to retrieve related data from ENSEMBL for {accession_base}"
            )
    return taxname, {gene: ebi_related_transcripts}


def filter_related(ID_base, related):
    filtered_genes = []

    for gene in related.get("genes", []):
        filtered_transcripts = []
        for transcript in gene.get("transcripts", []):
            has_tag = "tag" in transcript
            matching_providers = [
                provider
                for provider in transcript.get("providers", [])
                if ID_base == (provider.get("transcript_accession") or "").split(".")[0]
                or ID_base == (provider.get("protein_accession") or "").split(".")[0]
            ]

            if has_tag or matching_providers:
                filtered_transcript = dict(transcript)  # shallow copy
                if has_tag:
                    filtered_transcript["providers"] = transcript.get("providers", [])
                else:
                    filtered_transcript["providers"] = matching_providers
                filtered_transcripts.append(filtered_transcript)
        if filtered_transcripts:
            filtered_gene = dict(gene)
            filtered_gene["transcripts"] = filtered_transcripts
            filtered_genes.append(filtered_gene)

    return {
        "assemblies": related.get("assemblies"),
        "genes": filtered_genes,
    }


def _fetch_related_from_ncbi_dataset_report(gene_ids, timeout):
    base_url = ncbi_urls("Datasets_gene")
    gene_id_str = quote(",".join(map(str, gene_ids)))
    url = f"{base_url}/id/{gene_id_str}/dataset_report"
    genes_dataset_report_json = json.loads(request(url=url, timeout=timeout))
    taxname, dataset_related = _parse_dataset_report(genes_dataset_report_json)
    return taxname, dataset_related


def _fetch_related_from_ncbi_product_report(gene_ids, taxname, timeout):
    base_url = ncbi_urls("Datasets_gene")
    gene_id_str = quote(",".join(map(str, gene_ids)))
    url = f"{base_url}/id/{gene_id_str}/product_report"
    genes_product_report_json = json.loads(request(url=url, timeout=timeout))
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
        raise ValueError(f"Assembly accession could not be determined for {accession}")
    location_params = [f"{accession}:{start}-{end}" for start, end in locations]
    encoded_locations = [quote(loc) for loc in location_params]
    location_string = "&".join([f"locations={loc}" for loc in encoded_locations])
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
        genomic_related, product_related = _get_gene_related(gene_ids, timeout=timeout)
        related = _merge_related(genomic_related, product_related)
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

    dataset_json = json.loads(request(url=dataset_report_url, timeout=timeout))
    if dataset_json:
        taxname, genomic_related = _parse_dataset_report(dataset_json)

    product_related_url = (
        f"{base_url}/accession/{accessions_without_versions}/product_report"
    )
    product_json = json.loads(request(url=product_related_url, timeout=timeout))
    if product_json:
        product_related = _parse_product_report(product_json)
    if genomic_related or product_related:
        related = _merge_related(genomic_related, product_related)
        return taxname, related
    return None, None


def _get_related_ensembl(accession_base):
    taxname, ebi_related = _get_related_by_accession_from_ebi(accession_base)
    gene_symbol = next(iter(ebi_related))
    _, ncbi_related = _get_related_by_gene_symbol(gene_symbol, taxname)
    if ebi_related:
        related = merge_ncbi_ebi(ebi_related, ncbi_related)
        if taxname.upper() == "HOMO SAPIENS":
            related = filter_related(accession_base, related)
        return related


def _get_related_ncbi(accession, moltype, locations=[0, 0], timeout=30):
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
            _, ebi_gene_related = _get_related_by_accession_from_ebi(ensembl_gene)
            ebi_related.update(ebi_gene_related)
        if ebi_related:
            related = merge_ncbi_ebi(ebi_related, ncbi_related)
        if taxname and taxname.upper() == "HOMO SAPIENS":
            related = filter_related(accession_base, related)

    return related


def detect_sequence_source(seq_id: str):
    """
    Detects whether the input string is an Ensembl ID or a RefSeq accession.
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

    # ncbi archived historical references for refseq prefix at
    # https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly

    # Maybe here only gives refseq supported type?
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


def get_new_related(accession, locations=[0, 0]):
    """
    Input: A refseq identifier or human gene name.
    Output: A dictionary of related.
    """
    accession_base = accession.split(".")[0]
    related = {}

    source, moltype = detect_sequence_source(accession)

    if source == "ensembl":
        related = _get_related_ensembl(accession_base)
    elif source == "ncbi" and moltype != "unknown":
        related = _get_related_ncbi(accession, moltype, locations=locations)
    elif source == "other":
        _, related = _get_related_by_gene_symbol(accession)
    else:
        raise NameError(f"Could not retrieve related for {accession}.")

    return related


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get related sequences")
    parser.add_argument(
        "ID",
        help="A sequence ID, support gene name, ncbi or ebi accession.",
    )
    parser.add_argument(
        "locations",
        nargs="?",
        type=ast.literal_eval,  # ← safe parsing of Python literals
        help="A list of locations, e.g. [[112088000, 112088000]]",
    )

    args = parser.parse_args()

    print(json.dumps(get_new_related(args.ID, args.locations), indent=2))
