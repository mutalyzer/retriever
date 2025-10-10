import pytest
import json
from unittest.mock import MagicMock, Mock, patch
import json
import pytest
from pathlib import Path
from unittest.mock import MagicMock
from mutalyzer_retriever.request import Http400
from mutalyzer_retriever.reference import GRCH37
from mutalyzer_retriever.related_schema import related_schema
from mutalyzer_retriever.related import  get_related
from mutalyzer_retriever.util import DataSource, HUMAN_TAXON

class MockHttp400(Exception):
    def __init__(self, response):
        super().__init__("Bad Request")
        self.response = response

def normalize_enums(obj):
    # convert DataSource enums to its value
    if isinstance(obj, dict):
        return {
            k: normalize_enums(v)
            for k, v in obj.items()
        }
    elif isinstance(obj, list):
        return [normalize_enums(i) for i in obj]
    elif isinstance(obj, DataSource):
        return obj.value
    return obj


def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    with open(str(data_file), "r") as file:
        content = file.read()
    return content


def mock_get_grch37_chr_accession(chr):
    return GRCH37[chr]

def mock_get_accession_dataset_report(accession):
    return json.loads(_get_content(f"data/NCBI_dataset_report_{accession}.json"))

def mock_get_accession_product_report(accession):
    return json.loads(_get_content(f"data/NCBI_product_report_{accession}.json"))

def mock_get_gene_id_dataset_report(gene_ids):
    return json.loads(_get_content(f"data/NCBI_dataset_report_{gene_ids}.json"))

def mock_get_gene_id_product_report(gene_ids):
    return json.loads(_get_content(f"data/NCBI_product_report_{gene_ids}.json"))

def mock_get_gene_symbol_dataset_report(gene_symbol, taxon_name):
    return json.loads(_get_content(f"data/NCBI_dataset_report_{gene_symbol}.json"))

def mock_get_gene_symbol_product_report(gene_symbol, taxon_name):
    return json.loads(_get_content(f"data/NCBI_product_report_{gene_symbol}.json"))

def mock_get_genome_annotation_report(accession, locations):
    if isinstance(locations, list):
        locations = locations[0]
    return json.loads(_get_content(f"data/NCBI_annotation_report_{locations}.json"))

def mock_get_assembly_accession(chr):
    asssembly_map = {
        "NC_000011.10": "GCF_000001405.40",
        "NC_000011.9" : "GCF_000001405.25",
        "NC_060935.1" : "GCF_009914755.1"
    }
    return asssembly_map.get(chr)

def mock_lookup_symbol(gene_symbol, taxon_name=HUMAN_TAXON):
    return json.loads(_get_content(f"data/EBI_lookup_expand_1_{gene_symbol}.json"))

def mock_lookup_id(accession, expand):
    return json.loads(_get_content(f"data/EBI_lookup_expand_{expand}_{accession}.json"))


@pytest.fixture
def mock_ncbi_client(monkeypatch):
    with patch('mutalyzer_retriever.related.NCBIClient') as class_patch:
        mock_client_instance = Mock()
        class_patch.return_value = mock_client_instance

        mock_client_instance.get_accession_dataset_report = Mock(side_effect=mock_get_accession_dataset_report)
        mock_client_instance.get_accession_product_report = Mock(side_effect=mock_get_accession_product_report)
        mock_client_instance.get_gene_id_dataset_report = Mock(side_effect=mock_get_gene_id_dataset_report)
        mock_client_instance.get_gene_id_product_report = Mock(side_effect=mock_get_gene_id_product_report)
        mock_client_instance.get_gene_symbol_product_report = Mock(side_effect=mock_get_gene_symbol_product_report)
        mock_client_instance.get_gene_symbol_dataset_report = Mock(side_effect=mock_get_gene_symbol_dataset_report)
        mock_client_instance.get_assembly_accession = Mock(side_effect=mock_get_assembly_accession)
        mock_client_instance.get_genome_annotation_report = Mock(side_effect=mock_get_genome_annotation_report)

        yield mock_client_instance


@pytest.fixture
def mock_ensembl_client(monkeypatch):
    with patch('mutalyzer_retriever.related.EnsemblClient') as class_patch:
        mock_client_instance = Mock()
        class_patch.return_value = mock_client_instance

        mock_client_instance.lookup_symbol = Mock(side_effect=mock_lookup_symbol)
        mock_client_instance.lookup_id = Mock(side_effect=mock_lookup_id)

        yield mock_client_instance


@pytest.mark.parametrize("accession", ["ENST00000375549.8"])
def test_ensembl_mane_select_transcript(accession, mock_ncbi_client, mock_ensembl_client):
    """
    A MANE select ENSEMBL transcript, with a NCBI match.
    Expect chr accessions on three assemblies and one set of transcripts
    One is itself, also MANE Select.
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENST00000646730.1"])
def test_ensembl_mane_plus_clinical_transcript(accession, mock_ncbi_client, mock_ensembl_client):
    """
    A MANE Plus Clinical ENSEMBL transcript with a NCBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself and the other ones are from MANE Select.
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENST00000714087.1"])
def test_ensembl_transcript_no_ncbi_match_transcript(accession, mock_ncbi_client, mock_ensembl_client):
    """
    Not a MANE select ENSEMBL transcript, without NCBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ENSEMBL) and the other ones are from MANE Select.
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENST00000528048.5"])
def test_ensembl_transcript_with_ncbi_match_transcript(accession, mock_ncbi_client, mock_ensembl_client):
    """
    Not a MANE select ENSEMBL transcript, with a NCBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ENSEMBL and NCBI) and the other ones are from MANE Select.
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSMUST00000000175.6"])
def test_ensembl_mouse_transcript_with_ncbi_match_accession(accession, mock_ncbi_client, mock_ensembl_client):
    """
    An ENSEMBL mouse transcript, with a NCBI match.
    Expect a chr accessions on this mouse assemblies and one set of transcripts:
    From ENSEMBL and NCBI.
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENST00000530923.6"])
def test_ensembl_non_coding_transcript_with_ncbi_match_accession(accession, mock_ncbi_client, mock_ensembl_client):
    """
    Not a MANE select ENSEMBL transcript, non-coding, with a NCBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ENSEMBL and NCBI) and the other ones are from MANE Select.
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    #TODO: this matched up information is not the same from two sources



@pytest.mark.parametrize("accession", ["ENST00000714091.5"])
def test_ensembl_non_coding_transcript_without_ncbi_match_accession(accession, mock_ncbi_client, mock_ensembl_client):
    """
    Not a MANE select ENSEMBL transcript, non-coding, without a NCBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ENSEMBL) and the other ones are from MANE Select.
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSG00000204370.14"])
def test_ensembl_gene(accession, mock_ncbi_client, mock_ensembl_client):
    """
    An ENSEMBL gene,
    Expect chr accessions on three assemblies and one set of transcripts:
    One is MANE select from NCBI and ENSEMBL
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSG00000204370.10"])
def test_ensembl_gene_older_version(accession, mock_ncbi_client, mock_ensembl_client):
    """
    An ENSEMBL gene of an older version,
    Expect chr accessions on three assemblies and one set of transcripts:
    One is MANE select from NCBI and ENSEMBL. The same as with the latest version.
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_ENSG00000204370.14.json"))


@pytest.mark.parametrize("accession", ["ENSMUSG00000000171.6"])
def test_ensembl_mouse_gene(accession, mock_ncbi_client, mock_ensembl_client):
    """
    An ENSEMBL mouse gene,
    Expect chr accessions on one mouse assembly and one set of transcripts:
    From NCBI and ENSEMBL
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSDARG00000017744.10"])
def test_ensembl_zebrafish_gene(accession, mock_ncbi_client, mock_ensembl_client):
    """
    An ENSEMBL mouse gene,
    Expect chr accessions on one mouse assembly and one set of transcripts:
    From NCBI and ENSEMBL
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSP00000364699.3"])
def test_ensembl_mane_select_protein(accession, mock_ncbi_client, mock_ensembl_client):
    """
    An ENSEMBL protein ID, MANE Select
    Expect chr accessions on three assemblies and one set of transcripts:
    From NCBI and ENSEMBL
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSP00000519382.1"])
def test_ensembl_not_mane_select_protein(accession, mock_ncbi_client, mock_ensembl_client):
    """
    An ENSEMBL protein ID, non_MANE Select
    Expect chr accessions on three assemblies and two sets of transcripts:
    One is MANE Select from NCBI and ENSEMBL, the other one is itself
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["ENSE00003479002.1"])
def test_ensembl_invalid_exon_id_raises(accession, mock_ncbi_client, mock_ensembl_client):
    """
    An ENSEMBL exon ID, not support as input, check the error handling
    """
    with pytest.raises(ValueError, match=f"Unsupported molecule type: exon"):
        get_related(accession)


@pytest.mark.parametrize("accession", ["NM_003002.4"])
def test_ncbi_mane_select_transcript(accession, mock_ncbi_client, mock_ensembl_client):
    """
    A MANE select ncbi transcript, with a EBI match.
    Expect chr accessions on three assemblies and one set of transcripts
    One is itself, also MANE Select.
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NM_001374258.1"])
def test_ncbi_mane_plus_clinical_transcript(accession, mock_ncbi_client, mock_ensembl_client):
    """
    A MANE Plus Clinical ncbi transcript with a EBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself and the other ones are from MANE Select.
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["XM_024454345.2"])
def test_ncbi_transcript_no_ensembl_match_transcript(accession, mock_ncbi_client, mock_ensembl_client):
    """
    Not a MANE select ncbi transcript, without EBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ncbi) and the other ones are from MANE Select.
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NM_001276506.2"])
def test_ncbi_transcript_with_ensembl_match_transcript(accession, mock_ncbi_client, mock_ensembl_client):
    """
    Not a MANE select ncbi transcript, with a EBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ncbi and EBI) and the other ones are from MANE Select.
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NM_025848.3"])
def test_ncbi_mouse_transcript_with_ensembl_match_accession(accession, mock_ncbi_client, mock_ensembl_client):
    """
    An ncbi mouse transcript, with a EBI match.
    Expect a chr accessions on this mouse assemblies and one set of transcripts:
    From ncbi and EBI.
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession, locations", [
    ("NC_000011.10", "112086000_112088000"),
])
def test_ncbi_two_genes_at_hg38_chr_location(accession, locations, mock_ncbi_client, mock_ensembl_client):
    """
    A genomic range covers two genes .
    Expect chr accessions on three assemblies and two sets of MANE Select
    transcripts from NCBI and EBI.
    """
    related = get_related(accession, locations)
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}_{locations}.json"))


@pytest.mark.parametrize("accession, locations", [
    ("NC_000011.10", "112088000_112088100"),
])
def test_ncbi_one_gene_at_hg38_chr_location(accession, locations, mock_ncbi_client, mock_ensembl_client):
    """
    A genomic range covers two genes .
    Expect chr accessions on three assemblies and one set of MANE Select
    transcripts from NCBI and EBI.
    """
    related = get_related(accession, locations)
    assert related ==  json.loads(_get_content(f"data/related_{accession}_{locations}.json"))


@pytest.mark.parametrize("accession, locations", [
    ("NC_000011.10", "112096000_112100000"),
])
def test_ncbi_no_gene_at_hg38_chr_location(accession, locations, mock_ncbi_client, mock_ensembl_client):
    """
    A genomic range covers no genes .
    Expect chr accessions on three assemblies.
    """
    related = get_related(accession, locations)
    assert related ==  {}

@pytest.mark.parametrize("accession, locations", [
    ("NC_060935.1", "112097000_112100000"),
])
def test_ncbi_two_genes_at_t2t_chr_location(accession, locations, mock_ncbi_client, mock_ensembl_client):
    """
    A genomic range covers two genes .
    Expect chr accessions on three assemblies and two sets of MANE Select
    transcripts from NCBI and EBI.
    """
    related = get_related(accession, locations)
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}_{locations}.json"))


@pytest.mark.parametrize("accession, locations", [
    ("NC_000011.9", "111960000_111966000"),
])
def test_ncbi_one_gene_at_hg37_chr_location(accession, locations, mock_ncbi_client, mock_ensembl_client):
    """
    A genomic range covers two genes on hg37.
    Expect chr accessions on three assemblies and one set of MANE Select
    transcripts from NCBI and EBI.
    """
    related = get_related(accession, locations)
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}_{locations}.json"))


@pytest.mark.parametrize("accession", ["NR_077060.2"])
def test_ncbi_non_coding_transcript_with_ensembl_match_accession(accession, mock_ncbi_client, mock_ensembl_client):
    """
    Not a MANE select ncbi transcript, non-coding, with a EBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ncbi and EBI) and the other ones are from MANE Select.
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    #TODO: this matched up information is not the same from two sources


@pytest.mark.parametrize("accession", ["NR_157078.3"])
def test_ncbi_non_coding_transcript_without_ensembl_match_accession(accession, mock_ncbi_client, mock_ensembl_client):
    """
    Not a MANE select ncbi transcript, non-coding, without a EBI match.
    Expect chr accessions on three assemblies and multiple sets of transcripts:
    One is itself (from ncbi) and the other ones are from MANE Select.
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related ==  json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NP_002993.1"])
def test_ncbi_mane_select_protein(accession, mock_ncbi_client, mock_ensembl_client):
    """
    An ncbi protein ID, MANE Select
    Expect chr accessions on three assemblies and one set of transcripts:
    From EBI and ncbi
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NP_001263435.1"])
def test_ncbi_not_mane_select_protein(accession, mock_ncbi_client, mock_ensembl_client):
    """
    An ncbi protein ID, non_MANE Select
    Expect chr accessions on three assemblies and two sets of transcripts:
    One is MANE Select from EBI and ncbi, the other one is itself
    """
    related = normalize_enums(get_related(accession))
    assert related_schema.validate(related)
    assert related == json.loads(_get_content(f"data/related_{accession}.json"))


@pytest.mark.parametrize("accession", ["NP_000000.0"])
def test_invalid_ncbi_accession(accession, mock_ncbi_client, mock_ensembl_client):
    assert get_related(accession) == {}


@pytest.mark.parametrize("accession", ["ENSG00000000000.0"])
def test_invalid_ensembl_accession(accession, mock_ncbi_client, mock_ensembl_client):
    mock_error = Mock()
    mock_error.response = Mock()
    mock_error.response = Mock(status_code=400, text="Bad Request")

    mock_ensembl_client.look_id.side_effect = Http400(mock_error)
    assert get_related(accession) == {}


@pytest.mark.parametrize("accession", ["abcde"])
def test_invalid_gene_symbol(accession, mock_ncbi_client, mock_ensembl_client):
    assert get_related(accession) == {}


@pytest.mark.parametrize("accession, locations", [("NC_060935.1", "112097000_112100000_112100000")])
def test_invalid_chr_range(accession, locations, mock_ncbi_client, mock_ensembl_client):
    with pytest.raises(NameError, match = f"Invalid location format: '{locations}'. Expected format: point or range 'point;start_end'."):
        get_related(accession, locations)
