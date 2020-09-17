from Bio import SeqIO
from io import StringIO


def parse(fasta):
    records = []
    for record in SeqIO.parse(StringIO(fasta), "fasta"):
        records.append({"seq": str(record.seq), "description": record.description})
    if not records:
        raise ValueError
    if len(records) == 1:
        return records[0]
    # TODO: What to do when there are multiple records?
