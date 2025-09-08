import pytest
import json
from pathlib import Path
from unittest.mock import MagicMock
from mutalyzer_retriever.request import Http400
from mutalyzer_retriever.reference import GRCH37
from schema import SchemaError
from mutalyzer_retriever.new_related_schema import related_schema
from mutalyzer_retriever.new_related import get_new_related, _merge_related, _parse_dataset_report, _parse_product_report

class MockHttp400(Exception):
    def __init__(self, response):
        super().__init__("Bad Request")
        self.response = response



def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    with open(str(data_file), "r") as file:
        content = file.read()
    return content


def _fetch_ncbi_esummary(db, query_id, timeout=10):
    return json.loads(_get_content(f"data/esummary_{db}_{query_id}.json"))


def _get_grch37_chr_accession(chr, timeout=10):
    return GRCH37[chr]


def _fetch_related_from_ncbi_dataset_report(gene_ids, taxname, timout=10):
    genes_dataset_report_json = json.loads(_get_content(f"data/NCBI_dataset_report_geneid_{gene_ids}.json"))
    taxname, dataset_related = _parse_dataset_report(genes_dataset_report_json)
    return taxname, dataset_related    

def _fetch_related_from_ncbi_product_report(gene_ids, taxname, timeout=10):
    genes_product_report_json = json.loads(_get_content(f"data/NCBI_product_report_geneid_{gene_ids}.json"))
    return _parse_product_report(genes_product_report_json)

def _get_related_by_accession(accession, timeout):
    accession_without_versions = accession.split(".")[0]
    dataset_json = json.loads(_get_content(f"data/NCBI_dataset_report_{accession}.json"))
    product_json = json.loads(_get_content(f"data/NCBI_product_report_{accession}.json"))
    taxname, genomic_related= _parse_dataset_report(dataset_json)
    _, product_related = _parse_product_report(product_json)
    return taxname, _merge_related(genomic_related, product_related) 

def _get_related_by_gene_symbol(gene_symbol, timeout):
    dataset_json = json.loads(_get_content(f"data/NCBI_dataset_report_{gene_symbol}.json"))
    product_json = json.loads(_get_content(f"data/NCBI_product_report_{gene_symbol}.json"))
    taxname, genomic_related= _parse_dataset_report(dataset_json)
    product_related = _parse_product_report(product_json)
    return taxname, _merge_related(genomic_related, product_related) 


from unittest.mock import MagicMock
import json
import pytest

# Mock Http400 exception class that mimics real one
class MockHttp400(Exception):
    def __init__(self, response):
        super().__init__("Bad Request")
        self.response = response

@pytest.fixture
def monkeypatch_expand_fail(monkeypatch):
    def _fetch_ebi_lookup_grch38(accession_base, expand=1, timeout=10):
        if expand == 1:
            # For expand=1, raise error for proteins (ENSP)
            if accession_base.startswith("ENSP"):
                mock_response = MagicMock()
                mock_response.json.return_value = {
                    "error": "Expand option only available for Genes and Transcripts"
                }
                mock_response.status_code = 400
                raise MockHttp400(mock_response)
            else:
                return json.loads(_get_content(f"data/EBI_lookup_expand_1_{accession_base}.json"))
        elif expand == 0:
            return json.loads(_get_content(f"data/EBI_lookup_expand_0_{accession_base}.json"))
        else:
            raise ValueError(f"Unsupported expand value {expand}")
        
    monkeypatch.setattr("mutalyzer_retriever.new_related._fetch_ebi_lookup_grch38", _fetch_ebi_lookup_grch38)
    monkeypatch.setattr("mutalyzer_retriever.new_related._fetch_ncbi_esummary", _fetch_ncbi_esummary)
    monkeypatch.setattr("mutalyzer_retriever.new_related._get_grch37_chr_accession", _get_grch37_chr_accession)
    monkeypatch.setattr("mutalyzer_retriever.new_related._get_related_by_gene_symbol", _get_related_by_gene_symbol)
    monkeypatch.setattr("mutalyzer_retriever.new_related.Http400", MockHttp400)

    yield


@pytest.fixture
def monkeypatch_expand_success(monkeypatch):
    def _fetch_ebi_lookup_grch38(accession_base, expand=1, timeout=10):
        return json.loads(_get_content(f"data/EBI_lookup_expand_{expand}_{accession_base}.json"))
    monkeypatch.setattr("mutalyzer_retriever.new_related._fetch_ncbi_esummary", _fetch_ncbi_esummary)
    monkeypatch.setattr("mutalyzer_retriever.new_related._get_grch37_chr_accession", _get_grch37_chr_accession)
    # monkeypatch.setattr("mutalyzer_retriever.new_related._fetch_related_from_ncbi_product_report", _fetch_related_from_ncbi_product_report)
    # monkeypatch.setattr("mutalyzer_retriever.new_related._fetch_related_from_ncbi_dataset_report", _fetch_related_from_ncbi_dataset_report)
    # monkeypatch.setattr("mutalyzer_retriever.new_related._get_related_by_accession", _get_related_by_accession)
    monkeypatch.setattr("mutalyzer_retriever.new_related._get_related_by_gene_symbol", _get_related_by_gene_symbol)
    monkeypatch.setattr("mutalyzer_retriever.new_related._fetch_ebi_lookup_grch38", _fetch_ebi_lookup_grch38)

    yield


@pytest.mark.parametrize("accession", ["ENST00000375549.8"])
def test_ensembl_mane_select_transcript(accession, monkeypatch_expand_success):
    """
    A MANE select ENSEMBL transcript, with a NCBI match.
    Expect chr accessions on three assemblies and one set of transcripts
    One is itself, also MANE Select.
    """    
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENST00000646730.1"])
def test_ensembl_mane_plus_clinical_transcript(accession, monkeypatch_expand_success):
    """
    A MANE Plus Clinical ENSEMBL transcript with a NCBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself and the other ones are from MANE Select.
    """
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENST00000714087.1"])
def test_ensembl_transcript_no_ncbi_match_transcript(accession, monkeypatch_expand_success):
    """
    Not a MANE select ENSEMBL transcript, without NCBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ENSEMBL) and the other ones are from MANE Select.
    """
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENST00000528048.8"])
def test_ensembl_transcript_with_ncbi_match_transcript(accession, monkeypatch_expand_success):
    """
    Not a MANE select ENSEMBL transcript, with a NCBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ENSEMBL and NCBI) and the other ones are from MANE Select.
    """
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSMUST00000000175.6"])
def test_ensembl_mouse_transcript_with_ncbi_match_accession(accession, monkeypatch_expand_success):
    """
    An ENSEMBL mouse transcript, with a NCBI match.
    Expect a chr accessions on this mouse assemblies and one set of transcripts:
    From ENSEMBL and NCBI.
    """
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENST00000530923.6"])
def test_ensembl_non_coding_transcript_with_ncbi_match_accession(accession, monkeypatch_expand_success):
    """
    Not a MANE select ENSEMBL transcript, non-coding, with a NCBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ENSEMBL and NCBI) and the other ones are from MANE Select.
    """
    related = get_new_related(accession)
    assert related_schema.validate(related)
    #TODO: this matched up information is not the same from two sources



@pytest.mark.parametrize("accession", ["ENST00000714091.5"])
def test_ensembl_non_coding_transcript_without_ncbi_match_accession(accession, monkeypatch_expand_success):
    """
    Not a MANE select ENSEMBL transcript, non-coding, without a NCBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ENSEMBL) and the other ones are from MANE Select.
    """
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSG00000204370.14"])
def test_ensembl_gene(accession, monkeypatch_expand_success):
    """
    An ENSEMBL gene,
    Expect chr accessions on three assemblies and one set of transcripts:
    One is MANE select from NCBI and ENSEMBL
    """
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSG00000204370.10"])
def test_ensembl_gene_older_version(accession, monkeypatch_expand_success):
    """
    An ENSEMBL gene of an older version, 
    Expect chr accessions on three assemblies and one set of transcripts:
    One is MANE select from NCBI and ENSEMBL. The same as with the latest version.
    """
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_ENSG00000204370.14.json"))


@pytest.mark.parametrize("accession", ["ENSMUSG00000000171.6"])
def test_ensembl_mouse_gene(accession, monkeypatch_expand_success):
    """
    An ENSEMBL mouse gene,
    Expect chr accessions on one mouse assembly and one set of transcripts:
    From NCBI and ENSEMBL
    """
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSP00000364699.3"])
def test_ensembl_mane_select_protein(accession, monkeypatch_expand_fail):
    """
    An ENSEMBL protein ID, MANE Select
    Expect chr accessions on three assemblies and one set of transcripts:
    From NCBI and ENSEMBL
    """    
    related = get_new_related(accession)
    assert related_schema.validate(related)    
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSP00000519382.1"])
def test_ensembl_not_mane_select_protein(accession, monkeypatch_expand_fail):
    """
    An ENSEMBL protein ID, non_MANE Select
    Expect chr accessions on three assemblies and two sets of transcripts:
    One is MANE Select from NCBI and ENSEMBL, the other one is itself
    """    
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSE00003479002.1"])
def test_ensembl_invalid_exon_id_raises(accession, monkeypatch_expand_success):
    """
    An ENSEMBL exon ID, not support as input, check the error handling
    """
    with pytest.raises(ValueError, match=f"Failed to retrieve related data from ENSEMBL for {accession.split('.')[0]}"):
        get_new_related(accession)


@pytest.mark.parametrize("accession", ["NM_003002.4"])
def test_ncbi_mane_select_transcript(accession, monkeypatch_expand_success):
    """
    A MANE select ncbi transcript, with a EBI match.
    Expect chr accessions on three assemblies and one set of transcripts
    One is itself, also MANE Select.
    """    
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NM_001374258.1"])
def test_ncbi_mane_plus_clinical_transcript(accession, monkeypatch_expand_success):
    """
    A MANE Plus Clinical ncbi transcript with a EBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself and the other ones are from MANE Select.
    """
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["XM_024454345.2"])
def test_ncbi_transcript_no_ensembl_match_transcript(accession, monkeypatch_expand_success):
    """
    Not a MANE select ncbi transcript, without EBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ncbi) and the other ones are from MANE Select.
    """
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NM_001276506.2"])
def test_ncbi_transcript_with_ensembl_match_transcript(accession, monkeypatch_expand_success):
    """
    Not a MANE select ncbi transcript, with a EBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ncbi and EBI) and the other ones are from MANE Select.
    """
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NM_025848.3"])
def test_ncbi_mouse_transcript_with_ensembl_match_accession(accession, monkeypatch_expand_success):
    """
    An ncbi mouse transcript, with a EBI match.
    Expect a chr accessions on this mouse assemblies and one set of transcripts:
    From ncbi and EBI.
    """
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NR_077060.2"])
def test_ncbi_non_coding_transcript_with_ensembl_match_accession(accession, monkeypatch_expand_success):
    """
    Not a MANE select ncbi transcript, non-coding, with a EBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ncbi and EBI) and the other ones are from MANE Select.
    """
    related = get_new_related(accession)
    assert related_schema.validate(related)
    #TODO: this matched up information is not the same from two sources


@pytest.mark.parametrize("accession", ["NR_157078.3"])
def test_ncbi_non_coding_transcript_without_ensembl_match_accession(accession, monkeypatch_expand_success):
    """
    Not a MANE select ncbi transcript, non-coding, without a EBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ncbi) and the other ones are from MANE Select.
    """
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NP_002993.1"])
def test_ncbi_mane_select_protein(accession, monkeypatch_expand_fail):
    """
    An ncbi protein ID, MANE Select
    Expect chr accessions on three assemblies and one set of transcripts:
    From EBI and ncbi
    """    
    related = get_new_related(accession)
    assert related_schema.validate(related)    
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NP_001263435.1"])
def test_ncbi_not_mane_select_protein(accession, monkeypatch_expand_fail):
    """
    An ncbi protein ID, non_MANE Select
    Expect chr accessions on three assemblies and two sets of transcripts:
    One is MANE Select from EBI and ncbi, the other one is itself
    """    
    related = get_new_related(accession)
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))
    

# @pytest.mark.parametrize("accession, locations", [
#     ("NC_000011.10", [[112088970, 112088970], [100000000, 100000006]]),
# ])
# def test_current_genomic_accessions_with_hg38(accession, locations):
#     related = get_new_related(accession, locations)
#     #check schema
#     assert related_schema.validate(related)
#     #check assemblies
#     assert "assemblies" in related, "Missing 'assemblies' key"
#     assert len({a["accession"] for a in related.get("assemblies", [])}) == 3
#     #check related genes and also the order
#     assert related["genes"][0]["name"] == "CNTN5"
#     assert related["genes"][1]["name"] == "SDHD"


# @pytest.mark.parametrize("accession, locations", [
#     ("NC_000011.9", [[112088970, 112088970], [100000000, 100000006]]),
# ])
# def test_hg19_current_genomic_accesions(accession, locations):
#     related = get_new_related(accession, locations)
#     assert related_schema.validate(related)



# def test_t2t_current_genomic_accessions(accession, locations):


# def test_genomic_accession_without_version(accession, locations):


# def test_hg19_older_genomic_accessions(accession, locations):


# def test_mouse_genomic_accesions(accession, locations):


# def test_unsupported_accessions(accession, locations):


# def test_older_unsupported_accesions(accession, locations):


# def test_mane_select_transcript_accessions(accession):
#     # here should be accession from a gene has mane select 
#     # and mane plus clinical transcripts. e,g. BRAF


# def test_mane_plus


# def test_not_mane_select_transcript_accessions(accession):


# def test_older_transcript_aebi_json.get("Parent")ccessions(accesion):


# def test_mouse_transcript_accessions(accesion):


# def test_protein_accessions(accession):


# def test_gene_name(gene_name):


# def test_ensembl_transcript_accession(accesion):


# def test_older_ensembl_transcript_accession(accession):
    

# def test_ensembl_transcript_no_ncbi_match_accession(accession):


