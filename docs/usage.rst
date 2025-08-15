Usage
=====

This package provides a :doc:`command line interface <cli>`.

Configuration
-------------

If you plan to retrieve NCBI references, then create a configuration file where you
include the email address used to communicate with the NCBI:

.. code-block:: console

    $ echo EMAIL = your.email@address.com > config.txt

Optionally, include the NCBI API key:

.. code-block:: console

    $ echo NCBI_API_KEY = your_NCBI_key >> config.txt

Finally:

.. code-block:: console

    export MUTALYZER_SETTINGS="$(pwd)/config.txt"


Retrieve a reference
--------------------

To retrieve a reference mention its id with the ``--id`` option.

.. code-block:: console

    $ mutalyzer_retriever --id "NG_012337.1"
    ##sequence-region NG_012337.1 1 15948
    ...


Retrieve a reference model
--------------------------

To retrieve the reference model add ``--parse`` (``-p``). Optionally, choose the
preferred indentation with ``--indent``.

.. code-block:: console

    $ mutalyzer_retriever --id "NG_012337.1" -p --indent 2
    {
      "annotations": {
        "id": "NG_012337.1",
        "type": "record",
    ...


Output directory and split the model
------------------------------------

Specify an output directory with ``--output``.

.. code-block:: console

    $ mutalyzer_retriever --id "NG_012337.1" -p --indent 2 --output .
    $ less NG_012337.1
    {
      "annotations": {
        "id": "NG_012337.1",
        "type": "record",
    ...


Split the model between annotations and sequence with ``--split``.

.. code-block:: console

    $ mutalyzer_retriever --id "NG_012337.1" -p --indent 2 --output . --split
    $ less NG_012337.1.annotations
    {
      "id": "NG_012337.1",
      "type": "record",
    ...
    $ less NG_012337.1.sequence
    GGGCTTGGTTCTACCATATCTCTACTTTGTGTTTATGTTTGTGTATGCATGTACTCCAAAGTCTT
    ...


Choose the retrieval source
---------------------------

By default all the sources are accessed (in the following order: LRG, NCBI,
Esembl) and the reference is retrieved from the first one found. However,
a specific source can be specified with ``-source`` (``-s``).

.. code-block:: console

    $ mutalyzer_retriever --id "NG_012337.1" -s ncbi
    ...


Choose the retrieval file type
------------------------------

For NCBI and Ensembl the default retrieved reference is ``gff3``. However,
a ``fasta`` file can be also retrieved with ``--type`` (``-t``).

.. code-block:: console

    $ mutalyzer_retriever --id "NG_012337.1" -t fasta
    >NG_012337.1 Homo sapiens succinate dehydrogenase complex, ...
    GGGCTTGGTTCTACCATATCTCTACTTTGTGTTTATGTTTGTGTATGCATGTACTCCAA...
    ...

If ``--parse`` (``p``) is added to the previous command, the sequence model
is obtained (no annotations are included).

.. code-block:: console

    $ mutalyzer_retriever --id "NG_012337.1" -t fasta -p
    {"sequence": {"seq": "GGGCTTGGTTCTACCATATCTCTACTTT

For the moment, this is not the case when ``--parse`` (``p``) is used in
combination with ``-t gff3``.

Raw genbank files can be retrieved from NCBI with ``-t genbank``, but they
cannot be parsed to obtain a model.


Parse local files
-----------------

To obtain a model from local files (``gff3`` with ``fasta`` and ``lrg``) use
the ``from_file`` command.

.. code-block:: console

    $ mutalyzer_retriever from_file -h
    usage: mutalyzer_retriever from_file [-h]
                                         [--paths PATHS [PATHS ...]]
                                         [--is_lrg]

    optional arguments:
      -h, --help            show this help message and exit
      --paths PATHS [PATHS ...]
                            both gff3 and fasta paths or just an lrg
      --is_lrg              there is one file which is lrg

An example with ``gff3`` and ``fasta`` is as follows.

.. code-block:: console

    $ mutalyzer_retriever from_file --paths NG_012337.1.gff3 NG_012337.1.fasta
    {"annotations": {"id": "NG_012337.1", "type": "record", "location": ...
    ...

For an ``lrg`` file the ``--is_lrg`` flag needs to be added.

.. code-block:: console

    $ mutalyzer_retriever from_file --paths LRG_417 --is_lrg
    {"annotations": {"type": "record", "id": "LRG_417", "location": ...


Retrieve the NCBI reference models from FTP
-------------------------------------------

Starting from scratch, i.e., connect to the FTP location to retrieve the assembly
versions and to download the annotations files. Please note that the following
command will retrieve, besides the chromosomes (``NC_``), also the contigs
(``NT_``) and the scaffolds (``NW_``).

.. code-block:: console

    $ mutalyzer_retriever ncbi_assemblies
    Downloading assembly releases:
     - assembly: GRCh37
       - dir: ncbi_annotation_releases/GRCh37/20190906
       - dir: ncbi_annotation_releases/GRCh37/20220307
       - dir: ncbi_annotation_releases/GRCh37/20240902
     - assembly: GRCh38
       - dir: ncbi_annotation_releases/GRCh38/20180213
       ...
       - dir: ncbi_annotation_releases/GRCh38/20240823
     - assembly: T2T-CHM13v2
       - dir: ncbi_annotation_releases/T2T-CHM13v2/20230315
       - dir: ncbi_annotation_releases/T2T-CHM13v2/20231002
       - dir: ncbi_annotation_releases/T2T-CHM13v2/20240823
    Get annotation models:
    - get from: GRCh38, date: 20180213
      - NC_000001.11
      - NT_187361.1
      ...

To restrict only to specific reference ids and assuming that the input files are
already present in the ``./ncbi_annotation_releases`` (default) directory:

.. code-block:: console

    $ mutalyzer_retriever ncbi_assemblies --input ncbi_annotation_releases --ref_id_start NC_000023 --downloaded
    Using downloaded releases from:
     ./ncbi_annotation_releases
    Get annotation models:
    - get from: GRCh38, date: 20180213
      - NC_000023.11
    ...
    - get from: GRCh38, date: 20240823
      - NC_000023.11
    - get from: GRCh37, date: 20190906
      - NC_000023.10
    - get from: GRCh37, date: 20220307
      - NC_000023.10
    - get from: GRCh37, date: 20240902
      - NC_000023.10
    - get from: T2T-CHM13v2, date: 20230315
    - get from: T2T-CHM13v2, date: 20231002
    - get from: T2T-CHM13v2, date: 20240823
    - writing ./ncbi_annotation_models/NC_000023.11.annotations
    - writing ./ncbi_annotation_models/NC_000023.10.annotations


To restrict only to a specific reference id and an assembly id, with the input
files already present in the ``./ncbi_annotation_releases`` directory, and to
download also the sequences (``--include_sequence``) in the same directory:

.. code-block:: console

    $ mutalyzer_retriever ncbi_assemblies --ref_id_start NC_0 --assembly_id_start GRCh37 --downloaded --include_sequence
    Using downloaded releases from:
     ./ncbi_annotation_releases
    Get annotation models:
    - get from: GRCh37, date: 20190906
      - NC_000001.10
      ...
      - NC_000024.9
    - get from: GRCh37, date: 20220307
      - NC_000001.10
      ...
      - NC_000024.9
    - get from: GRCh37, date: 20240902
      - NC_000001.10
      ...
      - NC_000024.9
    - writing ./ncbi_annotation_models/NC_000001.10.annotations
      ...
    - writing ./ncbi_annotation_models/NC_000024.9.annotations
    Downloading the sequences:
    - get the sequence for NC_000001.10
    - writing ./ncbi_annotation_models/NC_000023.10.sequence
    ...
    - get the sequence for NC_000023.10
    - writing ./ncbi_annotation_models/NC_000023.10.sequence


Retrieve related reference ids
------------------------------

To obtain the related reference ids use the ``related`` flag.

.. code-block:: console

    $ mutalyzer_retriever --id LRG_303 --related --indent 2
    {
        "ncbi": [
            {
              "id": "NG_008376.4"
            },
            {
              "id": "AC254562.1"
            },
            {
              "id": "NM_000106.6"
            },
            {
              "id": "NR_034118.2"
            }
        ]
    }


