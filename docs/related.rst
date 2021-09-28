Related references
==================

To retrieve the related reference IDs the following steps are performed:

- Get a summary of the reference id from the NCBI using the Eutils `ESummary
  <https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_ESummary_>`_
  endpoint (for example `LRG_303 <https://eutils.ncbi.nlm.nih.gov/entrez/eutil
  s/esummary.fcgi?db=nucleotide&id=LRG_303&retmode=json>`_).

    - We consider the `nucleotide` database in the query.
    - If in the response the `accessionversion` is different than the input
      reference id, we include it within the related reference ids
      (`NG_008376.4` for the previous `LRG_303` example).
    - The `assemblyacc` value is also considered as related (`AC254562.1` for
      `LRG_303`).
    - If a new version is available, the `replacedby` key is present in the
      response (see `NG_012337.1 <https://eutils.ncbi.nlm.nih.gov/entrez/eutil
      s/esummary.fcgi?db=nucleotide&id=NG_012337.1&retmode=json>`_). In this
      case we use the ESummary endpoint recursively to retrieve all the newer
      versions.

- Get the NCBI linked uids by using the Eutils `ELink
  <https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_ELink_>`_ endpoint
  (for example `LRG_303 <https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
  elink.fcgi?db=nucleotide&dbfrom=nucleotide&id=LRG_303&cmd=neighbor&retmode=
  json>`__).

    - We use `nucleotide` in both database parameters (`db` and `dbfrom`).
    - We extract all the NCBI uids from the `linksetdbs`.

      - Note that for chromosomes we do not want all the transcript ids. For
        this reason, if `genome` equals `chromosome` in the previous ESummary
        response (for example `NC_000022.11 <https://eutils.ncbi.nlm.nih.gov
        /entrez/eutils/esummary.fcgi?db=nucleotide&id=NC_000022.11&ret
        mode=json>`_), we consider only those for which `nuccore_nuccore_comp`
        and `nuccore_nuccore_rsgb` appear as `linknames` in the ELink
        `response <https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?db
        =nucleotide&dbfrom=nucleotide&id=NC_000022.11&cmd=neighbor&retmode=jso
        n>`_.

- Use the ESummary endpoint with the NCBI uids extracted in the previous step
  to obtain their reference ids (for example `LRG_303 <https://eutils.ncbi.nlm
  .nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id=509155882,1435110251,5
  94191220&retmode=json>`__).

- If the reference id is from a transcript we make use of the `NCBI Datasets
  REST API <https://api.ncbi.nlm.nih.gov/datasets/docs/reference-docs/rest-api
  />`_ to obtain further related references (see `NM_003002.2 <https://api.ncb
  i.nlm.nih.gov/datasets/v1/gene/accession/NM_003002.2>`_).
