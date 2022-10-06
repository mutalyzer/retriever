Usage
=====

This package provides a :doc:`command line interface <cli>`.

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

Starting from scratch, i.e., connect to the FTP location to retrieve the
assembly versions and to download the annotations files.

.. code-block:: console

    $ mutalyzer_retriever ncbi_assemblies
    - local output directory set up to ./models
      done
      ...

Restrict only to specific reference ids and assuming that the input files are
already present in the ``downloads/`` directory.

.. code-block:: console

    $ mutalyzer_retriever ncbi_assemblies  --input downloads/ --ref_id_start NC_000023 --downloaded
    - local output directory set up to ./models
      done
    - processing 109 from 20180213, (GRCh38.p12, GCF_000001405.38)
      - NC_000023.11
    ...
    - processing 105.20220307 from 20220307, (GRCh37.p13, GCF_000001405.25)
      - NC_000023.10
    - writing ./models/NC_000023.10





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


