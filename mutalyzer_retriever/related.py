from .request import Http400, RequestErrors, request
import json


def _fetch_ncbi_esummary(db, query_id, timeout=10):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {"db": db, "id": query_id, "retmode": "json"}
    return json.loads(request(url=url, params=params, timeout=timeout))


def _fetch_ncbi_elink(db, dbfrom, query_id, timeout=10):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
    params = {
        "db": db,
        "dbfrom": dbfrom,
        "id": query_id,
        "cmd": "neighbor",
        "retmode": "json",
    }
    return json.loads(request(url=url, params=params, timeout=timeout))


def _fetch_ncbi_datasets_gene_accession(accession_id, timeout=1):
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v1/gene/accession/{accession_id}"
    return json.loads(request(url=url, timeout=timeout))


def _extract_link_uids(links, genome):
    linksetdbs = _extract(links, ["linksets", 0, "linksetdbs"])
    uid_links = set()
    if linksetdbs:
        for linksetdb in linksetdbs:
            if linksetdb.get("links"):
                if genome != "chromosome" or (
                    genome == "chromosome"
                    and linksetdb.get("linkname")
                    in ["nuccore_nuccore_comp", "nuccore_nuccore_rsgb"]
                ):
                    uid_links.update(linksetdb["links"])
    return list(uid_links)


def _get_summary_result_one(summary):
    if (
        summary.get("error")
        or not summary["result"].get("uids")
        or len(summary["result"]["uids"]) != 1
    ):
        return {}
    return summary["result"][summary["result"]["uids"][0]]


def _get_summary_accession_versions(summary):
    output = set()
    uids = _extract(summary, ["result", "uids"])
    if uids:
        for uid in uids:
            accession_version = _extract(summary, ["result", uid, "accessionversion"])
            if accession_version:
                output.add(accession_version)
    return output


def _get_new_versions(summary, timeout):
    new_versions = set()
    if summary.get("replacedby"):
        new_versions.add(summary["replacedby"])
        new_versions.update(
            _get_new_versions(
                _get_summary_result_one(
                    _fetch_ncbi_esummary("nucleotide", summary["replacedby"], timeout)
                ),
                timeout,
            )
        )
    return new_versions


def _get_linked_references(reference_id, genome, timeout):
    links = _fetch_ncbi_elink("nucleotide", "nucleotide", reference_id, timeout)
    link_uids = _extract_link_uids(links, genome)
    if link_uids:
        summary = _fetch_ncbi_esummary("nucleotide", ",".join(link_uids), timeout)
        # TODO: Make sure that request uri is not too long (414)
        return _get_summary_accession_versions(summary)
    return set()


def _fetch_ncbi_entrez_eutils_esummary(gene_id, timeout=1):
    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
        f"db=gene&id={gene_id}&retmode=json"
    )
    try:
        response = request(url=url, timeout=timeout)
    except RequestErrors:
        raise ConnectionError
    except Http400 as e:
        if "Failed to understand id" in e.response.text:
            raise NameError
        else:
            raise ConnectionError
    return json.loads(response)


def _extract(d, path):
    s_d = d
    for p in path:
        if (isinstance(s_d, dict) and s_d.get(p)) or (
            isinstance(s_d, list) and len(s_d) > p
        ):
            s_d = s_d[p]
        else:
            return None
    return s_d


def _add(d, k, v, not_none=True):
    if not_none:
        for v_v in v:
            if v_v is None:
                return
    if d.get(k) is None:
        d[k] = set()
    d[k].add(v)


def _update(d, k, s):
    if s:
        if d.get(k) is None:
            d[k] = set()
        d[k].update(s)


def _to_model(d):
    output = {}
    for k in d:
        output[k] = []
        # print(sorted(d[k], key=lambda r_id: r_id['id'] if r_id.get("selector") is None else (r_id['id'], r_id["selector"]["id"])))
        for i in d[k]:
            if len(i) == 1:
                output[k].append({"id": i[0]})
            else:
                output[k].append({"id": i[0], "selector": {"id": i[1]}})
    return output


def _extract_datasets(gene):
    related = {}
    paths = [
        ["genomic_ranges", 0, "accession_version"],
        ["reference_standards", 0, "gene_range", "accession_version"],
    ]
    for p in paths:
        _add(related, "ncbi", (_extract(gene, p),))
    transcripts = _extract(gene, ["transcripts"])
    if transcripts and isinstance(transcripts, list):
        for transcript in transcripts:
            _add(related, "ncbi", (_extract(transcript, ["accession_version"]),))
            _add(related, "ensembl", (_extract(transcript, ["ensembl_transcript"]),))
            if transcript.get("genomic_range") and transcript["genomic_range"].get(
                "accession_version"
            ):
                _add(
                    related,
                    "ncbi",
                    (
                        _extract(transcript, ["genomic_range", "accession_version"]),
                        _extract(transcript, ["accession_version"]),
                    ),
                )
    if gene.get("ensembl_gene_ids"):
        ensembl_genes = set()
        for ensembl_gene_id in gene.get("ensembl_gene_ids"):
            ensembl_genes.add((ensembl_gene_id,))
        _update(related, "ensembl", ensembl_genes)
    return related


def _extract_gene_summary(gene_summary, gene_id):
    related = set()
    locationhist = _extract(gene_summary, ["result", gene_id, "locationhist"])
    if locationhist and isinstance(locationhist, list):
        for l_h in locationhist:
            if l_h.get("chraccver"):
                related.add((l_h["chraccver"],))
    return related


def _get_ncbi_datasets_non_chromosome_related(reference_id, timeout=10):
    ncbi = _fetch_ncbi_datasets_gene_accession(reference_id, timeout)

    if ncbi.get("genes") and len(ncbi["genes"]) == 1 and ncbi["genes"][0].get("gene"):
        gene = ncbi["genes"][0]["gene"]
        related = _extract_datasets(gene)
        if gene.get("gene_id"):
            gene_summary = _fetch_ncbi_esummary("gene", gene.get("gene_id"), timeout)
            _update(
                related, "ncbi", _extract_gene_summary(gene_summary, gene["gene_id"])
            )
        return related
    return {}


def _get_related_from_summary(summary):
    related = set()
    if summary.get("accessionversion"):
        related.add(summary["accessionversion"])
    if summary.get("assemblyacc"):
        related.add(summary["assemblyacc"])
    return related


def _get_related_ncbi(reference_id, timeout=1):
    related = set()
    summary = _get_summary_result_one(
        _fetch_ncbi_esummary("nucleotide", reference_id, timeout)
    )
    if not summary:
        return {}

    related.update(_get_related_from_summary(summary))
    related.update(_get_new_versions(summary, timeout))
    related.update(_get_linked_references(reference_id, summary.get("genome"), timeout))
    related = {"ncbi": set([(i,) for i in related if i != reference_id])}
    if summary.get("biomol") in ["mRNA", "peptide", "ncRNA|lncRNA"]:
        datasets_related = _get_ncbi_datasets_non_chromosome_related(reference_id)
        for k in datasets_related:
            _update(related, k, datasets_related[k])
    return related


def get_related(reference_id, timeout=1):
    ncbi = _get_related_ncbi(reference_id, timeout)
    return _to_model(ncbi)
