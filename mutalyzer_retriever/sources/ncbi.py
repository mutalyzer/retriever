from http.client import HTTPException, IncompleteRead
from urllib.error import HTTPError

from Bio import Entrez

from ..configuration import settings
from ..request import Http400, RequestErrors, request
from ..util import f_e

Entrez.email = settings.get("EMAIL")
Entrez.api_key = settings.get("NCBI_API_KEY")


class ReferenceToLong(Exception):
    """
    Raised when the reference length exceeds maximum size.
    """

    pass


def fetch_ncbi_databases(reference_id):
    """
    Queries NCBI to identify in what databases a specific reference appears.

    Note:
    Whenever a reference without the most recent version is employed it seems
    that not all the databases are returned (including the most important ones
    for us, i.e., nuccore and protein). Hence, we always strip the reference
    and employ the accession only. Example:
    https://eutils.ncbi.nlm.nih.gov/gquery?term=NC_000001.10&retmode=xml
    versus:
    https://eutils.ncbi.nlm.nih.gov/gquery?term=NC_000001&retmode=xml
    Strange enough, but this works:
    https://www.ncbi.nlm.nih.gov/search/all/?term=NC_000001.10
    I was not able to find any additional parameter that would expand the
    search.

    :arg str reference_id: The id of the reference.
    :returns set: Set with NCBI databases.
    """
    if "." in reference_id:
        reference_id = reference_id.rsplit(".")[0]
    try:
        handle = Entrez.egquery(term=reference_id)
    except (IOError, HTTPError, HTTPException) as e:
        raise ConnectionError

    result = Entrez.read(handle)
    databases = set()
    for item in result["eGQueryResult"]:
        if item["Status"].upper() == "OK" and int(item["Count"]) >= 1:
            databases.add(item["DbName"])
    return databases


def _get_database(reference_id):
    try:
        databases = fetch_ncbi_databases(reference_id)
    except ConnectionError as e:
        raise e

    if "nuccore" in databases:
        return "nuccore"
    elif "protein" in databases:
        return "protein"
    elif "nucest" in databases:
        return "nucest"
        # Todo: Pay attention to the following:
        # https://ncbiinsights.ncbi.nlm.nih.gov/2018/07/30/upcoming-changes-est-gss-databases/
    else:
        raise NameError


def _fetch_reference_summary(reference_id):
    """
    Retrieves the reference summary if available on the NCBI.

    :arg str reference_id: The id of the reference.

    :returns dict:
    """
    try:
        handle = Entrez.esummary(id=reference_id)
    except (IOError, HTTPError, HTTPException):
        raise ConnectionError
    else:
        try:
            record = Entrez.read(handle)
        except RuntimeError:
            raise NameError
        else:
            handle.close()

    return {
        "reference_id": record[0]["AccessionVersion"],
        "db": _get_database(reference_id),
        "length": int(record[0]["Length"]),
    }


def fetch_genbank(reference_id, size_on=True):
    """
    Retrieve a genbank reference from the NCBI.

    :arg str reference_id: The id of the reference.
    :arg bool size_on: Consider or not the maximum sequence length.
    :returns: Reference content.
    :rtype: str
    """
    reference_summary = _fetch_reference_summary(reference_id)

    if size_on and reference_summary["length"] > settings["MAX_FILE_SIZE"]:
        raise ReferenceToLong
    try:
        handle = Entrez.efetch(
            db=reference_summary["db"],
            id=reference_id,
            rettype="gbwithparts",
            retmode="text",
        )
    except (IOError, HTTPError, HTTPException):
        raise NameError
    else:
        raw_data = handle.read()
        handle.close()
        return raw_data


def fetch_fasta(reference_id, db):
    """
    Retrieve the sequence of corresponding reference ID.

    :arg str reference_id: The reference ID.

    :returns str: The sequence.
    """
    try:
        handle = Entrez.efetch(db=db, id=reference_id, rettype="fasta")
    except HTTPError as e:
        if e.code == 400:
            # TODO: Check whether in the response is mentioned that the
            #  reference_id was not found.
            raise NameError(f_e("fasta", e))
        else:
            raise ConnectionError(f_e("fasta", e))
    except (IOError, HTTPException) as e:
        raise ConnectionError(f_e("fasta", e))
    else:
        try:
            raw_data = handle.read()
        except IncompleteRead as e:
            raw_data = e.partial.decode()
        handle.close()
        return raw_data


def fetch_gff3(reference_id, db, timeout=1):
    """
    Retrieve the gff3 for the corresponding reference ID.

    :arg str reference_id: The reference ID.

    :returns str: gff3 content.
    """
    url = settings["NCBI_GFF3_URL"]
    params = {"db": db, "report": "gff3", "id": reference_id}
    try:
        response = request(url=url, params=params, timeout=timeout)
    except RequestErrors as e:
        raise ConnectionError(f"(gff3) Original: {str(e)}")
    except Http400 as e:
        if "Failed to understand id" in e.response.text:
            raise NameError(f_e("gff3", e))
        else:
            raise ConnectionError(f_e("gff3", e))
    else:
        if response.startswith("Error"):
            raise NameError(f_e("gff3"))
        return response


def fetch(reference_id, reference_type, size_on=True, timeout=1):
    """
    Fetch the raw annotation for the corresponding reference ID.

    :arg str reference_id: The reference ID.
    :arg str reference_type: The reference type ("gff3" - default, or "genbank").
    :arg bool size_on: Consider maximum file size.

    :returns tuple: raw annotations, type ("gff3" or "genbank")
    """
    db = "nuccore"
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4702849/
    # https://support.nlm.nih.gov/knowledgebase/article/KA-03437/
    # https://support.nlm.nih.gov/knowledgebase/article/KA-03434/
    # https://support.nlm.nih.gov/knowledgebase/article/KA-03389/
    # https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/
    if (
        reference_id.startswith("AP_")
        or reference_id.startswith("NP_")
        or reference_id.startswith("WP_")
        or reference_id.startswith("XP_")
        or reference_id.startswith("YP_")
        or reference_id.startswith("ZP_")
    ):
        db = "protein"
    if reference_type in [None, "gff3"]:
        return fetch_gff3(reference_id, db, timeout), "gff3"
    elif reference_type == "fasta":
        return fetch_fasta(reference_id, db), "fasta"
    elif reference_type == "genbank":
        return fetch_genbank(reference_id, size_on), "genbank"

    raise ValueError(
        "NCBI fetch does not support '{}' reference type.".format(reference_type)
    )

