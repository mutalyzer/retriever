"""
Module for ENSEMBL lookup endpoint response parsing.
https://rest.ensembl.org/lookup/id/ENST00000530458?expand=1;content-type=application/json
"""
from mutalyzer_retriever.util import DataSource


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