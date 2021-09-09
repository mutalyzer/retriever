Retrieval sources
=================

NCBI
----

NCBI references are retrieved using the `Eutils sviewer <https://www.ncbi.nlm.
nih.gov/tools/sviewer/>`_ endpoint (`example <https://eutils.ncbi.nlm.nih.gov/
sviewer/viewer.cgi?report=gff3;id=NG_012337.3>`_). Note that for GRCh37
references the annotations are incomplete and can be manually retrieved from
the following `location <https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refs
eq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/ref_GRCh37.p13_top_level.gf
f3.gz>`__. They can be merged later with the sequence to obtain the model.

Ensembl
-------

Ensembl offers an `API <https://rest.ensembl.org/>`_ from where the most
recent reference versions can be retrieved. Queries are not accepted with the
version included, e.g., `ENST00000383925.1 <https://rest.ensembl.org/overlap/
id/ENST00000383925.1?feature=gene;feature=transcript;feature=exon;feature=cds;
content-type=application/json>`_, the version being part of the response,
e.g., `ENST00000383925 <https://rest.ensembl.org/overlap/id/ENST00000383925?fe
ature=gene;feature=transcript;feature=exon;feature=cds;content-type=applicatio
n/json>`_. For this reason we check if the provided reference id includes the
version, case in which we use the following `endpoint <https://rest.ensembl.or
g/lookup/id/ENST00000383925>`_ to check if the most recent version equals the
provided one. If not, for humans we check if the version matches the `GRCh37
dedicated API <https://grch37.rest.ensembl.org/>`_. If the reference has the
same id in GRCh37 and GRCh38 the retrieved one is from GRCh38.  The
`transcript archive <http://dev-tark.ensembl.org/api/>`_ may be employed in
future to retrieve other versions, but currently the annotation provided is
not complete.

Note that the retriever accepts only `stable ensembl ids <https://www.ensembl.
org/info/genome/stable_ids/index.html>`_, which start with ``ENS``.

LRG
---

LRG references are retrieved from the following `location <http://ftp.ebi.ac.uk
/pub/databases/lrgex/>`__.
