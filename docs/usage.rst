Usage
=====

This package provides a command line interface. To see the full options list,
use ``-h``.

.. code-block:: console

    $ mutalyzer_retriever -h
    usage: mutalyzer_retriever [-h] [-v]
                               [--id ID]
                               [-s {ncbi,ensembl,lrg}]
                               [-t {gff3,genbank,json,fasta}]
                               [-p]
                               [-m {all,sequence,annotations}
                               [--timeout TIMEOUT] [--indent INDENT] [--sizeoff] [-c CONFIGURATION]
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


Retrieve reference
------------------

.. code-block:: console

    $ mutalyzer_retriever --id "NG_012337.1"
    ##sequence-region NG_012337.1 1 15948
    ...

For NCBI and Ensembl references the gff3.


Retrieve reference model
------------------------

To retrieve the reference model add ``--parse`` (``-p``). Optionally, choose the
preferred indentation with ``--indent``.

.. code-block:: console

    $ mutalyzer_retriever --id "NG_012337.1" -p --indent 2
    {
      "annotations": {
        "id": "NG_012337.1",
        "type": "record",
    ...


Choose a retrieval file type
----------------------------



Choose a retrieval source
-------------------------

By default all the sources are accessed (in the following order: LRG, NCBI,
Esembl) and the reference is retrieved from the first one found. However,
a specific source can be specified with ``-source`` (``-s``).

.. code-block:: console

    $ mutalyzer_retriever --id "NG_012337.1" -s ncbi
    ...
