from .request import Http400, RequestErrors, request
import json


def _fetch_ncbi_datasets_gene_accession(accession_id, timeout=1):
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v1/gene/accession/{accession_id}"
    try:
        response = request(url=url, timeout=timeout)
    except RequestErrors:
        raise ConnectionError
    except Http400 as e:
        if "Failed to understand id" in e.response.text:
            raise NameError
        else:
            raise ConnectionError
    else:
        return response


def _fetch_ncbi_entrez_eutils_esummary(gene_id, timeout=1):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={gene_id}&retmode=json"
    try:
        response = request(url=url, timeout=timeout)
    except RequestErrors:
        raise ConnectionError
    except Http400 as e:
        if "Failed to understand id" in e.response.text:
            raise NameError
        else:
            raise ConnectionError
    else:
        return response


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


def _append(d, k, v, not_none=True):
    if not_none and v is not None:
        if d.get(k) is None:
            d[k] = []
        d[k].append(v)


def _extend(d, k, l):
    if l:
        if d.get(k) is None:
            d[k] = []
        d[k].extend(l)


def _merge(d1, d2):
    pass


def _extract_datasets(gene):
    related = {}
    paths = [
        ["genomic_ranges", 0, "accession_version"],
        ["reference_standards", 0, "gene_range", "accession_version"],
    ]
    for p in paths:
        _append(related, "ncbi", _extract(gene, p))
    transcripts = _extract(gene, ["transcripts"])
    if transcripts and isinstance(transcripts, list):
        for transcript in transcripts:
            _append(related, "ncbi", _extract(transcript, ["accession_version"]))
            _append(related, "ensembl", _extract(transcript, ["ensembl_transcript"]))
    if gene.get("ensembl_gene_ids"):
        _extend(related, "ensembl", gene.get("ensembl_gene_ids"))
    return related


def _extract_gene_summary(gene_summary, gene_id):
    related = set()
    locationhist = _extract(gene_summary, ["result", gene_id, "locationhist"])
    if locationhist and isinstance(locationhist, list):
        for l_h in locationhist:
            if l_h.get("chraccver"):
                related.add(l_h["chraccver"])
    return list(related)


def _clean(d):
    for k in d:
        d[k] = list(set(d[k]))


def get_related(reference_id, timeout=1):
    ncbi = json.loads(_fetch_ncbi_datasets_gene_accession(reference_id, timeout))
    if ncbi.get("genes") and len(ncbi["genes"]) == 1 and ncbi["genes"][0].get("gene"):
        gene = ncbi["genes"][0]["gene"]
        related = _extract_datasets(gene)
        if gene.get("gene_id"):
            gene_summary = json.loads(
                _fetch_ncbi_entrez_eutils_esummary(gene.get("gene_id"), timeout)
            )
            _extend(
                related, "ncbi", _extract_gene_summary(gene_summary, gene["gene_id"])
            )
        _clean(related)
    return related
