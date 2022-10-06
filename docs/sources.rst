Retrieval sources
=================

NCBI
----

References annotations are retrieved using the `Eutils sviewer <https://www.nc
bi.nlm.nih.gov/tools/sviewer/>`_ endpoint (`example <https://eutils.ncbi.nlm.n
ih.gov/sviewer/viewer.cgi?report=gff3;id=NG_012337.3>`_), while the fasta
sequences are retrieved using the Entrez API (`example <https://eutils.ncbi.nl
m.nih.gov/entrez/eutils/efetch.fcgi?db=sequences&id=NG_012337.3&rettype=fasta&
retmode=text>`__).

Assemblies
^^^^^^^^^^

For human chromosomal references (NC\_) the following `FTP location <https://f
tp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotatio
n_releases/>`__ is used to manually retrieve the annotations making sure that
the history is taken into account.

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
