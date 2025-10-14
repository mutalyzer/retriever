"""
Module for ENSEMBL lookup endpoint response parsing.
https://rest.ensembl.org/lookup/id/ENST00000530458?expand=1;content-type=application/json
"""

def parse_ensembl_gene_lookup_json(response):
    """
    Parse response from Ensembl lookup endpoint.

    Args:
        response (dict): Response fetched from Ensembl.

    Returns:
        output (dict): A dictionary containing parsed information for one Ensembl gene, including:
            - "taxon_name": Species name, in uppercase.
            - "name": Gene symbol/name.
            - "accession": Ensembl gene id if available.
            - "transcripts": List of Ensembl transcripts with available protein accessions and description.
                - "trancript_accession": Ensembl transcript id.
                - "description": Display name of the transcript
                - "protein_accession": Ensembl protein id.
    """
    gene_id = response.get("id")
    gene_symbol = response.get("display_name")

    if not (gene_id and gene_symbol):
        return {}

    output = {"name": gene_symbol, "accession": gene_id}

    taxon_name = response.get("species")
    if taxon_name:
        output["taxon_name"] = taxon_name.replace("_", " ").upper()

    transcripts = []
    for t in response.get("Transcript", []):
        t_id = t.get("id")
        t_version = t.get("version")
        t_desc = t.get("display_name")

        # Build transcript entry, it accession and description
        transcript = {}
        if t_id and t_version:
            transcript["transcript_accession"] = f"{t_id}.{t_version}"
        if t_desc:
            transcript["description"] = t_desc

        # Get protein accession
        translation = t.get("Translation", {})
        p_id = translation.get("id")
        p_version = translation.get("version")
        if p_id and p_version:
            transcript["protein_accession"] = f"{p_id}.{p_version}"

        if transcript:
            transcripts.append(transcript)

    if transcripts:
        output["transcripts"] = transcripts

    return output
