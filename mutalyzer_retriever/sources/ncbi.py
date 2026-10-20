from http.client import HTTPException, IncompleteRead
from urllib.error import HTTPError

from Bio import Entrez

from ..configuration import DEFAULT_TIMEOUT, settings
from ..request import Http400, RequestErrors, request
from ..util import f_e

Entrez.email = settings.get("EMAIL")
Entrez.api_key = settings.get("NCBI_API_KEY")


class ReferenceToLong(Exception):
    """
    Raised when the reference length exceeds maximum size.
    """


def _fetch_reference_summary(reference_id, db):
    """
    Retrieves the reference summary if available on the NCBI.

    :arg str reference_id: The id of the reference.
    :arg str db: The NCBI database ("nuccore" or "protein").

    :returns dict:
    """
    try:
        handle = Entrez.esummary(db=db, id=reference_id)
    except (IOError, HTTPError, HTTPException) as e:
        raise ConnectionError from e
    try:
        record = Entrez.read(handle)
    except RuntimeError as e:
        # ESummary reports both an unknown ID and an ID in another database
        # as errors. Other errors (e.g. a backend failure) are not lookup misses.
        if str(e).startswith(("Invalid uid ", "Otherdb uid=")):
            raise NameError(str(e)) from e
        raise ConnectionError(str(e)) from e
    except (IOError, HTTPException) as e:
        raise ConnectionError from e
    finally:
        handle.close()

    return {
        "reference_id": record[0]["AccessionVersion"],
        "length": int(record[0]["Length"]),
    }


def fetch_genbank(reference_id, db, size_on=True):
    """
    Retrieve a genbank reference from the NCBI.

    :arg str reference_id: The id of the reference.
    :arg str db: The NCBI database ("nuccore" or "protein").
    :arg bool size_on: Consider or not the maximum sequence length.
    :returns: Reference content.
    :rtype: str
    """
    reference_summary = _fetch_reference_summary(reference_id, db)

    if size_on and reference_summary["length"] > settings["MAX_FILE_SIZE"]:
        raise ReferenceToLong
    try:
        handle = Entrez.efetch(
            db=db,
            id=reference_id,
            rettype="gbwithparts",
            retmode="text",
        )
        try:
            raw_data = handle.read()
        finally:
            handle.close()
    except (IOError, HTTPError, HTTPException) as e:
        raise ConnectionError from e
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
            raise NameError(f_e("fasta", e)) from e
        raise ConnectionError(f_e("fasta", e)) from e
    except (IOError, HTTPException) as e:
        raise ConnectionError(f_e("fasta", e)) from e
    try:
        raw_data = handle.read()
    except IncompleteRead as e:
        raw_data = e.partial.decode()
    handle.close()
    if not raw_data.startswith(">"):
        # Ids NCBI cannot parse come back as HTTP 200 with a (letter-spaced)
        # "Error: Failed to understand id" body. Other errors are not misses.
        if raw_data.replace(" ", "").startswith("Error:Failedtounderstandid"):
            raise NameError(f_e("fasta", extra=raw_data.strip()))
        raise ConnectionError(f_e("fasta", extra=raw_data[:100].strip()))
    return raw_data


def fetch_gff3(reference_id, db, timeout=DEFAULT_TIMEOUT):
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
        raise ConnectionError(f"(gff3) Original: {str(e)}") from e
    except Http400 as e:
        if "Failed to understand id" in e.response.text:
            raise NameError(f_e("gff3", e)) from e
        raise ConnectionError(f_e("gff3", e)) from e
    if response.startswith("Error"):
        raise NameError(f_e("gff3"))
    return response


def fetch(reference_id, reference_type, size_on=True, timeout=DEFAULT_TIMEOUT):
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

    def _fetch(db):
        if reference_type in [None, "gff3"]:
            return fetch_gff3(reference_id, db, timeout), "gff3"
        if reference_type == "fasta":
            return fetch_fasta(reference_id, db), "fasta"
        if reference_type == "genbank":
            return fetch_genbank(reference_id, db, size_on), "genbank"
        raise ValueError(f"NCBI fetch does not support '{reference_type}' reference type.")

    try:
        return _fetch(db)
    except NameError:
        # Not found: the prefixes above miss e.g. INSDC proteins (BFP69721.1).
        return _fetch("protein" if db == "nuccore" else "nuccore")
