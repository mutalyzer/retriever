# A script to get related references from NCBI and EBI

import json
import requests
import jsonschema
import argparse
from urllib.parse import quote
from jsonschema import validate
from mutalyzer_retriever.related import _get_summary_result_one
from mutalyzer_retriever.request import Http400, RequestErrors, request


def NCBI_URLs(endpoint: str):
    URLs = {
        "ELink": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi",
        "ESummary": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        "EFetch": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        "Datasets_gene": "https://api.ncbi.nlm.nih.gov/datasets/v2/gene",
        "Datasets_genome": "https://api.ncbi.nlm.nih.gov/datasets/v2/genome"
    }
    return URLs.get(endpoint)

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
        ncbi_transcripts = ncbi_related["genes"][0]["transcripts"]
        transcript_index = {}
        for transcript in ncbi_transcripts:
            for provider in transcript.get("providers", []):
                tid = provider.get("transcript_accession")
                provider_name = provider.get("name")
                if tid and provider_name == "ENSEMBL":
                    transcript_index[tid] = transcript
        
        for ebi_transcript in ebi_related:
            ebi_tid = ebi_transcript.get('transcript_accession')
            if not ebi_tid:
                continue
            provider_entry = clean_dict({
                "name": "ENSEMBL",
                "transcript_accession": ebi_tid,
                "protein_accession": ebi_transcript.get("protein_accession"),
                "description": ebi_transcript.get("description")
            })
            
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
                new_transcript = clean_dict({
                    "providers": [provider_entry],
                    "tag": ebi_transcript.get("tag")
                })
                ncbi_transcripts.append(new_transcript)
    return related




def _get_grch37_chr_accession(chr, timeout=10):
    url = f"{NCBI_URLs('Datasets_genome')}/accession/GCF_000001405.25/sequence_reports?chromosomes={chr}&role_filters=assembled-molecule" 
    response = json.loads(request(url=url, timeout=timeout))
    return response.get("reports", [{}])[0].get("refseq_accession")


def _parse_dataset_report(json_report):
    """
    Parses NCBI gene dataset report and returns taxname and cleaned data for assemblies and genes.
    Args:
        data (dict): JSON data from NCBI dataset API.

    Returns:
        taxname,
        dict: {"assemblies": list, "genes": list}    
    """
    assemblies = []
    genes = []
    ensembl_genes = []
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

                clean_assembly = clean_dict({
                    "assembly_name": assembly_name,
                    "accession": genomic_acc
                })
                if clean_assembly:
                    assemblies.append(clean_assembly)

        # Add GRCh37 chromosome accession
        if taxname == "Homo sapiens":
            grch37_acc = _get_grch37_chr_accession(sequence_name)
            grch37_entry = clean_dict({
                "accession": grch37_acc,
                "accession_name": "GRCh37.p13"
            })

            if grch37_entry:
                assemblies.append(grch37_entry)

        # Handle reference standard accessions
        ref_acc = gene.get("reference_standards", [])
        hgnc_id_str = gene.get("nomenclature_authority", {}).get("identifier")
        hgnc_id = hgnc_id_str.split(':')[1] if hgnc_id_str else None

        ensembl_gene_id = ensembl_gene[0] if ensembl_gene else None
        ncbi_gene_id = gene.get("gene_id", [])

        # Build provider dicts, removing empty fields
        ncbi = clean_dict({
            "name": "NCBI",
            "accession": ncbi_gene_id,
        })

        ensembl = clean_dict({
            "name": "ENSEMBL",
            "accession": ensembl_gene_id,
        })

        providers = [p for p in (ncbi, ensembl) if len(p) > 1]
        if not providers:
            continue  # skip if both are empty        

        for ref in ref_acc:
            ref_accession = ref.get("gene_range", {}).get("accession_version")

        gene_entry = clean_dict({
            "name": symbol,
            "refseqgene": ref_accession,
            "hgnc_id": hgnc_id,
            "description": gene_description,
            "providers": providers
        })
        if gene_entry:
            genes.append(gene_entry)
    return taxname, {"assemblies": assemblies, "genes": genes}


def _parse_product_report(data):
    """
    Parse NCBI product report json into a dict:
    Args:
        data (dict): JSON data from NCBI product API.

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
            ncbi = clean_dict({
                "name": "NCBI",
                "transcript_accession": ncbi_transcript_acc,
                "protein_accession": ncbi_protein_acc,
                "description": ncbi_trancript_description
            })

            ensembl = clean_dict({
                "name": "ENSEMBL",
                "transcript_accession": ensembl_transcript_acc,
                "protein_accession": ensembl_protein_acc
            })

            providers = [p for p in (ncbi, ensembl) if len(p) > 1]
            if not providers:
                continue  # skip if both are empty

            product_entry = clean_dict({
                "providers": providers,
                "tag": transcript.get("select_category")
            })

            gene_products.append(product_entry)

        if gene_products:
            genes_dict[symbol] = gene_products

    return genes_dict


def _merge_related(genomic_related, product_related):
    """
    Merges genomic and product-related gene data.
    """
    related = {
        "assemblies": [clean_dict(a) for a in genomic_related.get("assemblies", []) if clean_dict(a)],
        "genes": []
    }
 
    for gene in genomic_related.get("genes", []):
        symbol = gene.get("name")
  

        transcripts = product_related.get(symbol, [])
 
        if not transcripts:
            continue  

        gene_info = clean_dict({
            "name": symbol,
            "hgnc_id": gene.get("hgnc_id"),
            "refseqgene": gene.get("refseqgene"),
            "description": gene.get("description"),
            "transcripts": transcripts
        })

        if gene_info:
            related["genes"].append(gene_info)

    # Sort genes alphabetically by name
    related["genes"].sort(key=lambda g: g.get("name", ""))

    return related


def _get_related_by_gene_symbol(gene_symbol,taxname="Homo Sapiens", timeout=10):
    """
    Given a gene symbol, return a set of related sequence accessions (genomic and/or products).
    Returns (taxname, related_dict), or (None, None) if nothing found.
    """
    if not gene_symbol:
        return None, None
    related = dict()
    genomic_related = None
    product_related = None

    base_url = NCBI_URLs("Datasets_gene")
    taxname_url_str = quote(taxname, safe="")
    dataset_report_url = f"{base_url}/symbol/{gene_symbol}/taxon/{taxname_url_str}/dataset_report"

    # Fetch and parse genomic related data
    dataset_json = json.loads(request(url=dataset_report_url, timeout=timeout))
    _, genomic_related = _parse_dataset_report(dataset_json) if dataset_json else (None, None)


    # Fetch and parse product related data
    product_url = f"{base_url}/symbol/{gene_symbol}/taxon/{taxname_url_str}/product_report"
    product_json = json.loads(request(url=product_url, timeout=timeout))
    product_related = _parse_product_report(product_json) if product_json else (None, None)

    if product_related or genomic_related:     
        related = _merge_related(genomic_related, product_related)
        return "Homo sapiens", related   
    return None, None



def get_uids(linkset:dict):
    uids_from_linkset = []
    for linksetdb in linkset["linksetdbs"]:
        uids_from_linkset.extend(linksetdb["links"])
    return uids_from_linkset


def _fetch_ncbi_esummary(db, query_id, timeout=10):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {"db": db, "id": query_id, "retmode": "json"}
    return json.loads(request(url=url, params=params, timeout=timeout))


def _get_summary_result_one(summary):
    if (
        summary.get("error")
        or not summary["result"].get("uids")
        or len(summary["result"]["uids"]) != 1
    ):
        return {}
    return summary["result"][summary["result"]["uids"][0]]


def _fetch_ebi_lookup_grch38(accession_base: str, expand= 1, timeout=10):
    if not accession_base:
        return
    url = f"https://rest.ensembl.org/lookup/id/{accession_base}?content-type=application/json;expand={expand}"
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

        ebi_related_transcripts.append({
            "name": "ENSEMBL",
            "transcript_accession": f"{transcript_id}.{transcript_version}",
            "protein_accession": protein_accession,
            "description": transcript.get("display_name")
        })

    return taxname, gene_symbol, ebi_related_transcripts    


def _parse_ebi_transcript_lookup_json(ebi_transcript_json):
    ebi_gene_json = _fetch_ebi_lookup_grch38(ebi_transcript_json.get("Parent"))
    taxname, gene_symbol, ebi_related_transcripts = _parse_ebi_lookup_json(ebi_gene_json)
    return taxname, gene_symbol, ebi_related_transcripts

def _parse_ebi_protein_lookup_json(ebi_protein_json):
    ebi_transcript_json = _fetch_ebi_lookup_grch38(ebi_protein_json.get("Parent"))
    taxname, gene_symbol, ebi_related_transcripts = _parse_ebi_transcript_lookup_json(ebi_transcript_json)
    return taxname, gene_symbol, ebi_related_transcripts
                                                   

def _parse_ebi_lookup_json(ebi_json):
    if ebi_json.get("object_type") == "Gene" and ebi_json.get("Transcript"):
        return _parse_ebi_gene_lookup_json(ebi_json)
    elif ebi_json.get("object_type") == "Transcript" and ebi_json.get("Parent"):
        return _parse_ebi_transcript_lookup_json(ebi_json)
    elif ebi_json.get("object_type") == "Translation" and ebi_json.get("Parent"):
        return _parse_ebi_protein_lookup_json(ebi_json)


def _get_related_by_accession_from_EBI(accession_base: str):
    gene = None
    ebi_related_transcripts = []
    taxname = None

    try:
        ebi_lookup_json = _fetch_ebi_lookup_grch38(accession_base)
    except Http400 as e:
        try:
            error_json = e.response.json()
            if error_json.get("error") == "Expand option only available for Genes and Transcripts":
                ebi_lookup_json = _fetch_ebi_lookup_grch38(accession_base, expand=0)
                parsed = _parse_ebi_lookup_json(ebi_lookup_json)
                if parsed:
                    taxname, gene, ebi_related_transcripts = parsed
                    return taxname, gene, ebi_related_transcripts
            else:
                raise
        except Exception as e:
            raise RuntimeError("Failed to parse error : {e}")

    if ebi_lookup_json and not ebi_lookup_json.get("error"):
        parsed = _parse_ebi_lookup_json(ebi_lookup_json)
        if parsed:
            taxname, gene, ebi_related_transcripts = parsed
        else:
            raise ValueError(f"Failed to retrieve related data from ENSEMBL for {accession_base}")
    return taxname, gene, ebi_related_transcripts


def filter_related(ID_base, related):
    filtered_genes = []

    for gene in related.get("genes", []):
        filtered_transcripts = []
        for transcript in gene.get("transcripts", []):
            has_tag = "tag" in transcript
            matching_providers = [
                provider for provider in transcript.get("providers", [])
                if ID_base == (provider.get("transcript_accession") or "").split(".")[0]
                or ID_base == (provider.get("protein_accession") or "").split(".")[0]
            ]

            if has_tag or matching_providers:
                filtered_transcript = dict(transcript) # shallow copy
                if has_tag:
                    filtered_transcript["providers"] = transcript.get("providers", [])
                else:
                    filtered_transcript["providers"] = matching_providers
                filtered_transcripts.append(filtered_transcript)
        if filtered_transcripts:
            filtered_gene = dict(gene)
            filtered_gene["transcripts"] = filtered_transcripts
            filtered_genes.append(filtered_gene)

    return {"assemblies": related.get("assemblies"), "genes": filtered_genes}



def _fetch_related_from_ncbi_dataset_report(gene_ids, timeout):
    base_url = NCBI_URLs("Datasets_gene")
    gene_id_str = quote(",".join(map(str, gene_ids)))
    url = f"{base_url}/id/{gene_id_str}/dataset_report"
    genes_dataset_report_json = json.loads(request(url=url, timeout=timeout))
    taxname, dataset_related = _parse_dataset_report(genes_dataset_report_json)
    return taxname, dataset_related


def _fetch_related_from_ncbi_product_report(gene_ids, taxname, timeout):
    base_url = NCBI_URLs("Datasets_gene")
    gene_id_str = quote(",".join(map(str, gene_ids)))
    url = f"{base_url}/id/{gene_id_str}/product_report"
    genes_product_report_json = json.loads(request(url=url, timeout=timeout))
    if genes_product_report_json:
        return _parse_product_report(genes_product_report_json)


def _get_gene_related(gene_ids, timeout):
    taxname, genomic_related = _fetch_related_from_ncbi_dataset_report(gene_ids, timeout)
    
    product_related = _fetch_related_from_ncbi_product_report(gene_ids, taxname, timeout)
    return genomic_related, product_related


def _get_related_by_chr_location(accession, locations, timeout):
    gene_ids = list()
    related = dict()
    assembly_accession = _get_assembly_accession(accession, timeout=timeout)
    if not assembly_accession:
        raise ValueError(f"Assembly accession could not be determined for {accession}")

    location_params = [f"{accession}:{start}-{end}" for start, end in locations]
    encoded_locations = [quote(loc) for loc in location_params]
    location_string = "&".join([f"locations={loc}" for loc in encoded_locations])
    url = f"{NCBI_URLs("Datasets_genome")}/accession/{assembly_accession}/annotation_report?{location_string}"
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


def _get_related_by_accession_from_NCBI(accessions, timeout):
    """
    Input: A sequence (transcripts or proteins) accession.
    Output: A set of related sequence accessions.
    """
    related = dict()
    genomic_related = None
    product_related = None

    base_url = NCBI_URLs("Datasets_gene")
    accessions_without_versions = accessions.split(".")[0]
    dataset_report_url = f"{base_url}/accession/{accessions_without_versions}/dataset_report"

    dataset_json = json.loads(request(url=dataset_report_url, timeout=timeout))
    if dataset_json:
        taxname, genomic_related = _parse_dataset_report(dataset_json)

    product_related_url = f"{base_url}/accession/{accessions_without_versions}/product_report"
    product_json = json.loads(request(url=product_related_url, timeout=timeout))
    if product_json:
        product_related = _parse_product_report(product_json)
    if genomic_related or product_related:
        related = _merge_related(genomic_related, product_related)
        return taxname,  related
    return None, None, None


def _get_related_ensembl(accession_base):
    taxname, gene_sysmbol, genes_ebi = _get_related_by_accession_from_EBI(accession_base)
    _, related_ncbi = _get_related_by_gene_symbol(gene_sysmbol, taxname)
    if genes_ebi:
        related = merge_ncbi_ebi(genes_ebi, related_ncbi)
        if taxname.upper() == "HOMO SAPIENS":
            related = filter_related(accession_base, related)
        return related


def get_new_related(ID, locations=[0,0], timeout=30):
    """
    Input: A refseq identifier or human gene name.
    Output: A dictionary of related.
    """

    ID_base = ID.split(".")[0]

    related = {}   

    # Ensembl declares its identifiers should be in the form of 
    # ENS[species prefix][feature type prefix][a unique eleven digit number]
    # See at https://www.ensembl.org/info/genome/stable_ids/index.html
    if ID.startswith("ENS"):
        related = _get_related_ensembl(ID_base)


    for db in ["nucleotide", "protein"]:
        summary = _get_summary_result_one(_fetch_ncbi_esummary(db, ID_base, timeout))
        if summary:
            break
    
    if summary:
        moltype = summary.get("moltype", "")
        try: 
            if "RNA" or "AA" in moltype.upper():
                taxname, related = _get_related_by_accession_from_NCBI(ID, timeout)
            elif "DNA" in moltype.upper() and "NC_" in ID:
                taxname, related = _get_related_by_chr_location(ID, locations, timeout)
            else:
                raise NameError(f"Could not retrieve related for {ID}.")
        except Exception as e:
            raise RuntimeError(f"Error fetching related accessions: {e}")            
    else:
        # not a valid NCBI/EBI accession , try with gene name
        try:
            _, related = _get_related_by_gene_symbol(ID, timeout)
            if not related:
                raise NameError(f"Could not retrieve related accessions for {ID}")
        except Exception as e:
            raise RuntimeError(f"Error fetching related accessions: {e}")
        
  





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get related sequences")
    parser.add_argument("ID", help="A sequence ID, support gene name, NCBI or EBI accession.")
    parser.add_argument("locations", nargs="?", help="A lis of locations")

    args = parser.parse_args()
