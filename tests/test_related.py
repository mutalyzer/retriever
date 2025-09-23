import pytest
import json
from unittest.mock import MagicMock, Mock
import json
import pytest
from pathlib import Path
from unittest.mock import MagicMock
from mutalyzer_retriever.request import Http400
from mutalyzer_retriever.reference import GRCH37
from mutalyzer_retriever.related_schema import related_schema
from mutalyzer_retriever.related import get_related, _merge_datasets, _parse_dataset_report, _parse_product_report, _parse_genome_annotation_report

class MockHttp400(Exception):
    def __init__(self, response):
        super().__init__("Bad Request")
        self.response = response



def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    with open(str(data_file), "r") as file:
        content = file.read()
    return content


def _get_grch37_chr_accession(chr):
    return GRCH37[chr]


def _get_assembly_accession(chr):
    asssembly_map = {
        "NC_000011.10": "GCF_000001405.40",
        "NC_000011.9" : "GCF_000001405.25",
        "NC_060935.1" : "GCF_009914755.1"
    }
    return asssembly_map.get(chr)


def _get_related_by_chr_location(accession, locations):
    taxname = None
    related = {}
    genome_annotation_json = json.loads(_get_content(f"data/NCBI_annotation_report_{accession}_{locations}.json"))
    taxname, related = _parse_genome_annotation_report(genome_annotation_json)
    return taxname, related


def _fetch_related_from_ncbi_product_report(gene_ids):
    genes_product_report_json = json.loads(_get_content(f"data/NCBI_product_report_{gene_ids}.json"))
    taxname, product_related = _parse_product_report(genes_product_report_json)
    return taxname, product_related


def _fetch_related_from_ncbi_dataset_report(gene_ids):
    genes_dataset_report_json = json.loads(_get_content(f"data/NCBI_dataset_report_{gene_ids}.json"))
    taxname, dataset_related = _parse_dataset_report(genes_dataset_report_json)
    return taxname, dataset_related

def _get_genomic_related_from_datasets(gene_ids):
    genes_dataset_report_json = json.loads(_get_content(f"data/NCBI_dataset_report_gene_id_{gene_ids}.json"))
    taxname, dataset_related = _parse_dataset_report(genes_dataset_report_json)
    return taxname, dataset_related    


def _get_product_related_from_datasets(gene_ids):
    genes_product_report_json = json.loads(_get_content(f"data/NCBI_product_report_gene_id_{gene_ids}.json"))
    taxname, product_related = _parse_product_report(genes_product_report_json)
    return taxname, product_related


def _get_related_by_accession_from_ncbi(accession):
    accession_without_versions = accession.split(".")[0]
    dataset_json = json.loads(_get_content(f"data/NCBI_dataset_report_{accession}.json"))
    product_json = json.loads(_get_content(f"data/NCBI_product_report_{accession}.json"))
    taxname, genomic_related= _parse_dataset_report(dataset_json)
    _, product_related = _parse_product_report(product_json)
    return taxname, _merge_datasets(genomic_related, product_related) 


def _get_related_by_gene_symbol_from_ncbi(gene_symbol, taxname):
    dataset_json = json.loads(_get_content(f"data/NCBI_dataset_report_{gene_symbol}.json"))
    product_json = json.loads(_get_content(f"data/NCBI_product_report_{gene_symbol}.json"))
    taxname, genomic_related= _parse_dataset_report(dataset_json)
    taxname, product_related = _parse_product_report(product_json)
    return taxname, _merge_datasets(genomic_related, product_related)


@pytest.fixture
def monkeypatch(monkeypatch):
    def _fetch_ebi_lookup_grch38(accession_base, expand=1):
        return json.loads(_get_content(f"data/EBI_lookup_expand_{expand}_{accession_base}.json"))

    monkeypatch.setattr("mutalyzer_retriever.related._get_assembly_accession", _get_assembly_accession)
    monkeypatch.setattr("mutalyzer_retriever.related._fetch_related_from_ncbi_dataset_report", _fetch_related_from_ncbi_dataset_report)
    monkeypatch.setattr("mutalyzer_retriever.related._fetch_related_from_ncbi_product_report", _fetch_related_from_ncbi_product_report)            
    monkeypatch.setattr("mutalyzer_retriever.related._get_grch37_chr_accession", _get_grch37_chr_accession)
    monkeypatch.setattr("mutalyzer_retriever.related._get_product_related_from_datasets", _get_product_related_from_datasets)
    monkeypatch.setattr("mutalyzer_retriever.related._get_genomic_related_from_datasets", _get_genomic_related_from_datasets)
    monkeypatch.setattr("mutalyzer_retriever.related._get_related_by_accession_from_ncbi", _get_related_by_accession_from_ncbi)
    monkeypatch.setattr("mutalyzer_retriever.related._get_related_by_gene_symbol_from_ncbi", _get_related_by_gene_symbol_from_ncbi)
    monkeypatch.setattr("mutalyzer_retriever.related._get_related_by_chr_location", _get_related_by_chr_location)
    monkeypatch.setattr("mutalyzer_retriever.related._fetch_ebi_lookup_grch38", _fetch_ebi_lookup_grch38)

    yield


@pytest.mark.parametrize("accession", ["ENST00000375549.8"])
def test_ensembl_mane_select_transcript(accession, monkeypatch):
    """
    A MANE select ENSEMBL transcript, with a NCBI match.
    Expect chr accessions on three assemblies and one set of transcripts
    One is itself, also MANE Select.
    """    
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENST00000646730.1"])
def test_ensembl_mane_plus_clinical_transcript(accession, monkeypatch):
    """
    A MANE Plus Clinical ENSEMBL transcript with a NCBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself and the other ones are from MANE Select.
    """
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENST00000714087.1"])
def test_ensembl_transcript_no_ncbi_match_transcript(accession, monkeypatch):
    """
    Not a MANE select ENSEMBL transcript, without NCBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ENSEMBL) and the other ones are from MANE Select.
    """
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENST00000528048.5"])
def test_ensembl_transcript_with_ncbi_match_transcript(accession, monkeypatch):
    """
    Not a MANE select ENSEMBL transcript, with a NCBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ENSEMBL and NCBI) and the other ones are from MANE Select.
    """
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSMUST00000000175.6"])
def test_ensembl_mouse_transcript_with_ncbi_match_accession(accession, monkeypatch):
    """
    An ENSEMBL mouse transcript, with a NCBI match.
    Expect a chr accessions on this mouse assemblies and one set of transcripts:
    From ENSEMBL and NCBI.
    """
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENST00000530923.6"])
def test_ensembl_non_coding_transcript_with_ncbi_match_accession(accession, monkeypatch):
    """
    Not a MANE select ENSEMBL transcript, non-coding, with a NCBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ENSEMBL and NCBI) and the other ones are from MANE Select.
    """
    related = get_related(accession)
    assert related_schema.validate(related)
    #TODO: this matched up information is not the same from two sources



@pytest.mark.parametrize("accession", ["ENST00000714091.5"])
def test_ensembl_non_coding_transcript_without_ncbi_match_accession(accession, monkeypatch):
    """
    Not a MANE select ENSEMBL transcript, non-coding, without a NCBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ENSEMBL) and the other ones are from MANE Select.
    """
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSG00000204370.14"])
def test_ensembl_gene(accession, monkeypatch):
    """
    An ENSEMBL gene,
    Expect chr accessions on three assemblies and one set of transcripts:
    One is MANE select from NCBI and ENSEMBL
    """
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSG00000204370.10"])
def test_ensembl_gene_older_version(accession, monkeypatch):
    """
    An ENSEMBL gene of an older version, 
    Expect chr accessions on three assemblies and one set of transcripts:
    One is MANE select from NCBI and ENSEMBL. The same as with the latest version.
    """
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_ENSG00000204370.14.json"))


@pytest.mark.parametrize("accession", ["ENSMUSG00000000171.6"])
def test_ensembl_mouse_gene(accession, monkeypatch):
    """
    An ENSEMBL mouse gene,
    Expect chr accessions on one mouse assembly and one set of transcripts:
    From NCBI and ENSEMBL
    """
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSP00000364699.3"])
def test_ensembl_mane_select_protein(accession, monkeypatch):
    """
    An ENSEMBL protein ID, MANE Select
    Expect chr accessions on three assemblies and one set of transcripts:
    From NCBI and ENSEMBL
    """    
    related = get_related(accession)
    assert related_schema.validate(related)    
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSP00000519382.1"])
def test_ensembl_not_mane_select_protein(accession, monkeypatch):
    """
    An ENSEMBL protein ID, non_MANE Select
    Expect chr accessions on three assemblies and two sets of transcripts:
    One is MANE Select from NCBI and ENSEMBL, the other one is itself
    """    
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSE00003479002.1"])
def test_ensembl_invalid_exon_id_raises(accession, monkeypatch):
    """
    An ENSEMBL exon ID, not support as input, check the error handling
    """
    with pytest.raises(ValueError, match=f"Unsupported molecule type: exon"):
        get_related(accession)


@pytest.mark.parametrize("accession", ["NM_003002.4"])
def test_ncbi_mane_select_transcript(accession, monkeypatch):
    """
    A MANE select ncbi transcript, with a EBI match.
    Expect chr accessions on three assemblies and one set of transcripts
    One is itself, also MANE Select.
    """    
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NM_001374258.1"])
def test_ncbi_mane_plus_clinical_transcript(accession, monkeypatch):
    """
    A MANE Plus Clinical ncbi transcript with a EBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself and the other ones are from MANE Select.
    """
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["XM_024454345.2"])
def test_ncbi_transcript_no_ensembl_match_transcript(accession, monkeypatch):
    """
    Not a MANE select ncbi transcript, without EBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ncbi) and the other ones are from MANE Select.
    """
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NM_001276506.2"])
def test_ncbi_transcript_with_ensembl_match_transcript(accession, monkeypatch):
    """
    Not a MANE select ncbi transcript, with a EBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ncbi and EBI) and the other ones are from MANE Select.
    """
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NM_025848.3"])
def test_ncbi_mouse_transcript_with_ensembl_match_accession(accession, monkeypatch):
    """
    An ncbi mouse transcript, with a EBI match.
    Expect a chr accessions on this mouse assemblies and one set of transcripts:
    From ncbi and EBI.
    """
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession, locations", [
    ("NC_000011.10", "112086000-112088000"),
])
def test_ncbi_two_genes_at_hg38_chr_location(accession, locations, monkeypatch):
    """
    A genomic range covers two genes .
    Expect chr accessions on three assemblies and two sets of MANE Select
    transcripts from NCBI and EBI.
    """
    related = get_related(accession, locations)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}_{locations}.json"))


@pytest.mark.parametrize("accession, locations", [
    ("NC_000011.10", "112088000-112088100"),
])
def test_ncbi_one_gene_at_hg38_chr_location(accession, locations, monkeypatch):
    """
    A genomic range covers two genes .
    Expect chr accessions on three assemblies and one set of MANE Select
    transcripts from NCBI and EBI.
    """
    related = get_related(accession, locations)
    assert related == json.loads(_get_content(f"data/related_{accession}_{locations}.json"))    


@pytest.mark.parametrize("accession, locations", [
    ("NC_000011.10", "112096000-112100000"),
])
def test_ncbi_no_gene_at_hg38_chr_location(accession, locations, monkeypatch):
    """
    A genomic range covers no genes .
    Expect chr accessions on three assemblies.
    """
    related = get_related(accession, locations)
    assert related == {}



@pytest.mark.parametrize("accession, locations", [
    ("NC_060935.1", "112097000-112100000"),
])
def test_ncbi_two_genes_at_t2t_chr_location(accession, locations, monkeypatch):
    """
    A genomic range covers two genes .
    Expect chr accessions on three assemblies and two sets of MANE Select
    transcripts from NCBI and EBI.
    """
    related = get_related(accession, locations)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}_{locations}.json"))


@pytest.mark.parametrize("accession, locations", [
    ("NC_000011.9", "111960000-111966000"),
])
def test_ncbi_one_gene_at_hg37_chr_location(accession, locations, monkeypatch):
    """
    A genomic range covers two genes on hg37.
    Expect chr accessions on three assemblies and one set of MANE Select
    transcripts from NCBI and EBI.
    """
    related = get_related(accession, locations)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}_{locations}.json"))


@pytest.mark.parametrize("accession", ["NR_077060.2"])
def test_ncbi_non_coding_transcript_with_ensembl_match_accession(accession, monkeypatch):
    """
    Not a MANE select ncbi transcript, non-coding, with a EBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ncbi and EBI) and the other ones are from MANE Select.
    """
    related = get_related(accession)
    assert related_schema.validate(related)
    #TODO: this matched up information is not the same from two sources


@pytest.mark.parametrize("accession", ["NR_157078.3"])
def test_ncbi_non_coding_transcript_without_ensembl_match_accession(accession, monkeypatch):
    """
    Not a MANE select ncbi transcript, non-coding, without a EBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ncbi) and the other ones are from MANE Select.
    """
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NP_002993.1"])
def test_ncbi_mane_select_protein(accession, monkeypatch):
    """
    An ncbi protein ID, MANE Select
    Expect chr accessions on three assemblies and one set of transcripts:
    From EBI and ncbi
    """    
    related = get_related(accession)
    assert related_schema.validate(related)    
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NP_001263435.1"])
def test_ncbi_not_mane_select_protein(accession, monkeypatch):
    """
    An ncbi protein ID, non_MANE Select
    Expect chr accessions on three assemblies and two sets of transcripts:
    One is MANE Select from EBI and ncbi, the other one is itself
    """    
    related = get_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))
