Usage
=====

The package provides a command line interface. Use ``-h`` to see the available options.

.. code-block:: console

    $ mutalyzer_retriever -h
    usage: mutalyzer_retriever [-h] [-v]
                               [--id ID]
                               [-s {ncbi,ensembl,lrg}]
                               [-t {gff3,genbank,json,fasta}]
                               [-p]
                               [-m {all,sequence,annotations}
                               [--timeout TIMEOUT]
                               [--indent INDENT]
                               [--sizeoff]
                               [-c CONFIGURATION]
                               {from_file} ...

    Mutalyzer genomic reference retriever.

    positional arguments:
      {from_file}           parse files to get the model

    optional arguments:
      -h, --help            show this help message and exit
      -v                    show program's version number and exit
      --id ID               the reference id
      -s {ncbi,ensembl,lrg}, --source {ncbi,ensembl,lrg}
                            retrieval source
      -t {gff3,genbank,json,fasta}, --type {gff3,genbank,json,fasta}
                            reference type
      -p, --parse           parse reference content
      -m {all,sequence,annotations}, --model_type {all,sequence,annotations}
                            include the complete model or parts of it
      --timeout TIMEOUT     timeout
      --indent INDENT       indentation spaces
      --sizeoff             do not consider file size
      -c CONFIGURATION, --configuration CONFIGURATION
                            configuration file path


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
