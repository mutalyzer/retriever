import pytest

from mutalyzer_retriever.reference import (
    get_assembly_chromosome_accession,
    get_assembly_id,
)


@pytest.mark.parametrize(
    "assembly_id,chromosome,expected",
    [
        ("GRCh38", "chr1", "NC_000001.11"),
        ("hg19", "1", "NC_000001.10"),
        ("T2T-CHM13v2.0", "chr1", "NC_060925.1"),
        ("hs1", "chr1", "NC_060925.1"),
        ("GRCh38", "chrM", "NC_012920.1"),
        ("hg38", "MT", "NC_012920.1"),
        ("t2t", "chrM", "CP068254.1"),
        ("GRCh37", "MT", "NC_012920.1"),
        ("GRCh37", "chrM", "NC_012920.1"),
        ("hg19", "chrMT", "NC_012920.1"),
        ("hg19", "chrM", "NC_012920.1"),
    ],
)
def test_get_assembly_chromosome_accession(assembly_id, chromosome, expected):
    assert get_assembly_chromosome_accession(assembly_id, chromosome) == expected


@pytest.mark.parametrize(
    "assembly_id,expected",
    [("T2T-CHM13v2.0", "T2T"), ("hs1", "T2T"), ("hg38", "GRCH38"), ("hg18", None)],
)
def test_get_assembly_id(assembly_id, expected):
    assert get_assembly_id(assembly_id) == expected
