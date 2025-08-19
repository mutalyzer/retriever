# A script to get related references using Datasets or Entrez

import json
import requests
import jsonschema
from urllib.parse import quote
from jsonschema import validate
from mutalyzer_retriever.related import _get_summary_result_one
from mutalyzer_retriever.request import Http400, RequestErrors, request



def _get_current_NCBI_version(summary):
    if summary.get("accessionversion"):
        return summary.get("accessionversion")


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


def _get_current_EBI_version(accessions_without_version, timeout=10):
    url = f"https://rest.ensembl.org/lookup/id/{accessions_without_version}?content-type=application/json"
    response = json.loads(request(url=url, timeout=timeout))
    if response:
        return f"{accessions_without_version}.{response.get('version')}"


def NCBI_URLs(endpoint: str):
    URLs = {
        "ELink": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi",
        "ESummary": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        "EFetch": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        "Datasets_gene": "https://api.ncbi.nlm.nih.gov/datasets/v2/gene",
        "Datasets_genome": "https://api.ncbi.nlm.nih.gov/datasets/v2/genome"
    }
    return URLs.get(endpoint)


def get_uids(linkset:dict):
    uids_from_linkset = []
    for linksetdb in linkset["linksetdbs"]:
        uids_from_linkset.extend(linksetdb["links"])
    return uids_from_linkset


def _get_grch37_chr_accession(chr, timeout):
    url = f"{NCBI_URLs("Datasets_genome")}/accession/GCF_000001405.25/sequence_reports?chromosomes={chr}&role_filters=assembled-molecule" 
    response = json.loads(request(url=url, timeout=timeout))
    return response.get("reports", [{}])[0].get("refseq_accession")


def _parse_dataset_report(json_report):
    """
    Groups unique genomic_accession_version and reference standard accession_version values by gene symbol.
    Converts genomic accession versions into a list of dicts with 'assembly_name' and 'accession'.
    Also returns the taxname as a dictionary with an 'organism' key.

    Args:
        json_report (dict): JSON data returned from the NCBI API.

    Returns:
        dict: A dictionary with 'taxname' as a dictionary and 'assemblies' as a list of dicts.
    """
    assemblies = []
    refseqgenes = []
    taxname = None

    # Iterate through the reports to extract assembly and accession version information
    for report in json_report.get("reports", []):
        gene = report.get("gene", {})
        symbol = gene.get("symbol", [])
        
        # Extract taxname (organism name)
        if not taxname:
            taxname = gene.get("taxname", "Unknown")

        # Get genomic accessions and convert them into the desired format
        for annotation in gene.get("annotations", []):
            for loc in annotation.get("genomic_locations", []):
                genomic_acc = loc.get("genomic_accession_version")
                if genomic_acc:
                    assembly_name = annotation.get("assembly_name")
                    assemblies.append({
                        "assembly_name": assembly_name,
                        "accession": genomic_acc
                    })

        # Optionally, handle reference standard accessions if needed (for completeness)
        ref_acc = gene.get("reference_standards", [])
        hgnc_id_str = gene.get('nomenclature_authority', {}).get('identifier', None)
        hgnc_id = hgnc_id_str.split(':')[1] if hgnc_id_str else None
        for ref in ref_acc:
            ref_accession = ref.get("gene_range", {}).get("accession_version", None)
            
            if ref_accession:
                refseqgenes.append({
                    "name": symbol,
                    "accession": ref_accession,
                    "hgnc_id": hgnc_id
                })

    return taxname, {"assemblies":assemblies, "refseqgene":refseqgenes}




def _parse_product_report(data):
    """
    Returns a dictionary mapping gene symbol to a list of transcript mappings:
    Each transcript mapping is a dictionary:
    {
        "transcript": transcript_acc,
        "ensembl_transcript": ensembl_transcript,
        "protein": protein_acc,
        "ensembl_protein": ensembl_protein,
        "Tag": "mane_select" (if present)
    }

    Args:
        data (dict): JSON data from NCBI API.

    Returns:
        dict: Mapping gene_symbol -> list of transcript info dicts
    """
    genes_dict = {}

    for report in data.get("reports", []):
        product = report.get("product", {})
        transcripts = product.get("transcripts", [])
        symbol = product.get("symbol", "UNKNOWN")
        gene_products = []

        for transcript in transcripts:
            
            gene_transcript = []
            transcript_acc = transcript.get("accession_version")
            protein = transcript.get("protein", {})

            if not transcript_acc:
                break

            # if transcript.get("select_category") == "MANE_SELECT":
            #     gene_product["tag"] = "mane_select"

            gene_transcript.append({"name":"NCBI", "transcript_id":transcript.get("accession_version"), "protein_id":protein.get("accession_version")})
            gene_transcript.append({"name":"ENSEMBL", "transcript_id":transcript.get("ensembl_transcript"), "protein_id":protein.get("ensembl_protein")})
            # gene_product["providers"] = gene_transcript
            gene_products.append({"providers":gene_transcript, "tag":transcript.get("select_category")})


        genes_dict[symbol] = gene_products


    return genes_dict



def _fetch_related_from_ncbi_dataset_report(gene_ids, timeout):
    base_url = NCBI_URLs("Datasets_gene")
    gene_id_str = quote(",".join(map(str, gene_ids)))
    url = f"{base_url}/id/{gene_id_str}/dataset_report"
    genes_dataset_report_json = json.loads(request(url=url, timeout=timeout))
    taxname, dataset_relted = _parse_dataset_report(genes_dataset_report_json)
    return taxname, dataset_relted


def _fetch_related_from_ncbi_product_report(gene_ids, taxname, timeout):
    base_url = NCBI_URLs("Datasets_gene")
    gene_id_str = quote(",".join(map(str, gene_ids)))
    url = f"{base_url}/id/{gene_id_str}/product_report"
    genes_product_report_json = json.loads(request(url=url, timeout=timeout))
    if genes_product_report_json:
        return _get_closest_products(_parse_product_report(genes_product_report_json),None, taxname=="Homo sapiens")   


def _get_gene_related(gene_ids, timeout):
    taxname, genomic_related = _fetch_related_from_ncbi_dataset_report(gene_ids, timeout)
    
    product_related = _fetch_related_from_ncbi_product_report(gene_ids, taxname, timeout)
    return genomic_related, product_related


def _get_closest_products(product_related, accession_without_version, human):

    """
    Filters product data for each gene to include only entries where:
    - the transcript matches the given accession (ignoring version), or
    - the product is marked as "mane_select".

    Removes duplicate product entries.

    Args:
        product_related (dict): Mapping of gene -> list of product dicts.
        accession_without_version (str): Transcript accession without version (e.g., "NM_003002").

    Returns:
        dict: Filtered mapping with only relevant, deduplicated product entries.
    """
    if not human:
        return product_related
    else:
        closest_related = {}
 
        for gene, products in product_related.items():
            filtered_products = []

            for product in products:
                transcript_ids = [item['transcript_id'].split(".")[0] for item in product.get("providers") if item.get('transcript_id') is not None]

                select = product.get("tag")
                filters = (
                    select != "None" or
                    (transcript_ids and accession_without_version in transcript_ids)
                )

                if filters:
                    filtered_products.append(product)
            if filtered_products:
                closest_related[gene] = filtered_products


    return closest_related


def sort_related(related):
    #TODO:  sort by key gene_name
    assemblies = related.pop("Assemblies")
    sorted_related = {"Assemblies": assemblies}
    for key in sorted(related.keys()):
        sorted_related[key] = related[key]
    return sorted_related



def _merge_related(genomic_related, product_related):
    related = {}

    related["assemblies"] = genomic_related["assemblies"]
    related["genes"] = []

    # Add gene-related info and corresponding products
    for gene in genomic_related["refseqgene"]:
        # Prepare the gene information for the final result
        gene_info = {
            "name": gene["name"],
            "hgnc_id": gene["hgnc_id"],
            "transcripts": product_related[gene["name"]],
            "providers": [{"name":"NCBI", "accession": gene.get("accession")}]
        }
        related["genes"].append(gene_info)

    # Add taxname (if it exists in genomic data)
    taxname = genomic_related.get("taxname", None)
    if taxname:
        related["taxname"] = {"organism": taxname}
    return related
    # Sort the final dictionary (if needed)
    return sort_related(related)




def _get_related_by_accession(accessions, timeout):
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
        product_related = _get_closest_products(_parse_product_report(product_json), accessions_without_versions, taxname=="Homo sapiens")
    if genomic_related or product_related:
        related = _merge_related(genomic_related, product_related)
        return taxname, related
    return None, None


def _get_related_by_gene_symbol(gene_symbol, timeout):
    """
    Input: A sequence (transcripts or proteins) accession.
    Output: A set of related sequence accessions.
    """
    related = dict()
    genomic_related = None
    product_related = None


    # Build request to NCBI ELink
    base_url = NCBI_URLs("Datasets_gene")
    
    dataset_report_url = f"{base_url}/symbol/{gene_symbol.upper()}/taxon/9606/dataset_report"


    dataset_json = json.loads(request(url=dataset_report_url, timeout=timeout))
    if dataset_json:
        taxname, genomic_related = _parse_dataset_report(dataset_json)

    product_related_url = f"{base_url}/symbol/{gene_symbol.upper()}/taxon/9606/product_report"
    product_json = json.loads(request(url=product_related_url, timeout=timeout))
    if product_json:
        product_related = _get_closest_products(_parse_product_report(product_json), None, human=True)
    if product_related or genomic_related:
        related = _merge_related(genomic_related, product_related)
        return taxname, related   
    return None, None


def _get_assembly_accession(accession, timeout):
    url = f"{NCBI_URLs("Datasets_genome")}/sequence_accession/{accession}/sequence_assemblies"
    assembly_json = json.loads(request(url=url, timeout=timeout))
    accessions = assembly_json.get("accessions")
    if isinstance(accessions, list) and accessions:
        return accessions[0]
    return None


def _get_related_by_accession_and_location(accession, locations, timeout):
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


from schema import SchemaError
from mutalyzer_retriever.new_related_schema import query_schema

def get_new_related(accession, locations=None, timeout=30):
    if locations is None:
        locations = [0, 0]

    accession_base = accession.split(".")[0]

    query_output = {
        "query": f"{accession}+{locations}",
        "organism": "Unknown",
        "moltype": "Unknown",
        "accession": accession,
        "related": {
            "assemblies": [],
            "genes": []
        }
    }

    summary = _get_summary_result_one(
        _fetch_ncbi_esummary("nucleotide", accession_base, timeout)
    )

    if summary:
        moltype = summary.get("moltype", "")
        organism = summary.get("organism", "Unknown")

        query_output["moltype"] = moltype
        query_output["organism"] = organism
        query_output["accession"] = summary.get("accessionversion", accession)

        try:
            if "RNA" in moltype.upper():
                _, related = _get_related_by_accession(accession, timeout=timeout)
            elif "DNA" in moltype.upper() and "NC_" in accession:
                _, related = _get_related_by_accession_and_location(accession, locations, timeout)
            else:
                raise NameError(f"Could not retrieve related accessions for {accession}")
        except Exception as e:
            raise RuntimeError(f"Error fetching related accessions: {e}")
    else:
        # Fallbacks: EBI or gene symbol search
        organism, related = _get_related_by_accession(accession, timeout=timeout)
        if related:
            query_output["organism"] = organism
            query_output["accession"] = _get_current_EBI_version(accession_base, timeout=timeout)
        else:
            organism, related = _get_related_by_gene_symbol(accession, timeout)
            if related:
                query_output["organism"] = organism
                query_output["accession"] = accession
            else:
                raise NameError(f"Could not retrieve related accessions for {accession}")

    query_output["related"] = related


    # Validate against schema
    try:
        query_schema.validate(query_output)
    except SchemaError as e:
        print("Schema validation error:", e)
        raise NameError("Output did not conform to schema")

    return query_output

            


if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description="Get related from Entrez or Datasets")
    # parser.add_argument("accession", help="Sequence ID")
    # parser.add_argument("location", help="Location")

    # args = parser.parse_args
    test_list = [
        ["NC_000011.10",  [[112088970, 112088970], [100000000, 100000006]]], # hg38, t2t; two sets of mane select accessions
        ["NC_000011.9",  [[112088970, 112088970], [100000000, 100000006]]],  # hg38, t2t; two sets of mane select accessions
        # ["NC_000011.8",  [[112088970, 112088970], [100000000, 100000006]]],  # hg38, t2t; two sets of mane select accessions
        # # # ["NG_012337.3",  [[274, 277]]], # NG_ itself, since no way to retrieve for NG_, NW_, NT
        # # ["NT_033899.8",  [[1000, 1010]]], # current version, NT_033899.9
        # ["NM_003002.2", [[200, 202], [100, 110]]], # NC_ from hg38, T2T; two sets of Mane select of current versions
        # ["NM_003002.4", [[200, 202], [100, 110]]], # NC_ from hg38, T2T; one set of Mane select of current version
        ["NM_001276506.2", [[200, 202], [100, 110]]], # NC_ from hg38, T2T; two sets: one Mane select, one of NM_001276506
        ["SDHD", [[200, 202], [100, 110]]], # NC_ from hg38, T2T; one set of Mane select of current versions
        ["Sdhd", [[200, 202], [100, 110]]],  # NC_ from hg38, T2T; one set of Mane select of current versions     
        ["ENST00000375549.8", [[200, 202], [100, 110]]], # NC_ from hg38, T2T; one set of Mane select of current version
        ["ENST00000375549.7", [[200, 202], [100, 110]]], # NC_ from hg38, T2T; one set of Mane select of current version
        ["ENST00000375549", [[200, 202], [100, 110]]], # NC_ from hg38, T2T; one set of Mane select of current version
        ["ENSP00000364699.3", [[200, 202], [100, 110]]], # NC_ from hg38, T2T; one set of Mane select of current version
        ["NC_000075.7",  [[50507000,50510000]]], # mouse model, related 
    ]
    import pprint

    for test_sample in test_list:
        pprint.pprint(get_new_related(test_sample[0], test_sample[1]))
