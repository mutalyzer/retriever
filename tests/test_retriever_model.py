import json
from pathlib import Path

import pytest

from mutalyzer_retriever import retrieve_model
from mutalyzer_retriever.parsers.gff3 import parse
from mutalyzer_retriever.retriever import (
    NoReferenceError,
    ReferenceParseError,
    retrieve_model_from_file,
)

from .commons import _get_content, _retrieve_raw, patch_retriever, references


def get_tests(references):

    tests = []

    for r_source in references:
        for r_type in references[r_source]:
            for r_id in references[r_source][r_type]:
                if r_type == "json":
                    p = Path(__file__).parent / "data" / str(r_id + ".tark.model.json")
                else:
                    p = Path(__file__).parent / "data" / str(r_id + ".model.json")
                with p.open() as f:
                    r_model = json.loads(f.read())
                tests.append(
                    pytest.param(
                        r_id,
                        r_source,
                        r_type,
                        r_model,
                        id=f"{r_source}-{r_type}-{r_id}",
                    )
                )

    return tests


def _seq_from_rest(r_id):
    return _get_content("data/" + str(r_id) + ".sequence")


def load_expected_model(r_id):
    with (Path(__file__).parent / "data" / f"{r_id}.model.json").open() as f:
        return json.load(f)


@pytest.mark.parametrize(
    "r_id, r_source, r_type, expected_model", get_tests(references)
)
def test_model(r_id, r_source, r_type, expected_model, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.parsers.json_ensembl._seq_from_rest",
        lambda _0, _1, _2, _3, _4, _5: _seq_from_rest(r_id),
    )
    assert retrieve_model(r_id, r_source, r_type) == expected_model


def test_model_type_sequence():
    # model_type only changes what retrieve_model() returns, not per-reference
    # behaviour, so a single fixture is enough to cover the "gff3" branch.
    r_id = "NM_000077.4"
    expected_model = load_expected_model(r_id)

    assert retrieve_model(r_id, "ncbi", "gff3", model_type="sequence") == {
        "sequence": expected_model["sequence"]
    }


def test_model_type_annotations():
    r_id = "NM_000077.4"
    expected_model = load_expected_model(r_id)

    assert (
        retrieve_model(r_id, "ncbi", "gff3", model_type="annotations")
        == expected_model["annotations"]
    )


def test_retrieve_model_from_file_gff3_and_fasta():
    data_dir = Path(__file__).parent / "data"
    model = retrieve_model_from_file(
        paths=[str(data_dir / "AA010203.1.gff3"), str(data_dir / "AA010203.1.fasta")]
    )
    with (data_dir / "AA010203.1.model.json").open() as f:
        expected_model = json.load(f)
    assert model == expected_model


def test_retrieve_model_from_file_lrg():
    data_dir = Path(__file__).parent / "data"
    model = retrieve_model_from_file(paths=[str(data_dir / "LRG_11.lrg")], is_lrg=True)
    with (data_dir / "LRG_11.model.json").open() as f:
        expected_model = json.load(f)
    assert model == expected_model


def test_unparsable_gff3():
    # NC_037304.1: CDS lines outside the lines of their trans-spliced nad5 gene.
    with pytest.raises(ReferenceParseError):
        retrieve_model("NC_037304.1", "ncbi", "gff3")


def test_translation_exception_split_codon():
    # A codon split by an intron is given as join(...), which crashed the parsing
    # of whole chromosomes in production (Xenopus, zebrafish, chicken, cow, etc.).
    # Trimmed lines of NC_141022.1 (zebrafish), XP_087097748.1.
    content = (
        "##gff-version 3\n"
        "##sequence-region NC_141022.1 1 63929756\n"
        "NC_141022.1\tRefSeq\tregion\t1\t63929756\t.\t+\t.\tID=NC_141022.1:1..63929756;gbkey=Src\n"
        "NC_141022.1\tGnomon\tgene\t53243316\t53268470\t.\t+\t.\tID=gene-LOC101882419;Name=LOC101882419;gbkey=Gene;gene=LOC101882419;gene_biotype=protein_coding\n"
        "NC_141022.1\tGnomon\tmRNA\t53243316\t53268470\t.\t+\t.\tID=rna-XM_087241622.1;Parent=gene-LOC101882419;Name=XM_087241622.1;gbkey=mRNA;gene=LOC101882419;transcript_id=XM_087241622.1\n"
        "NC_141022.1\tGnomon\tCDS\t53245300\t53245374\t.\t+\t0\tID=cds-XP_087097748.1;Parent=rna-XM_087241622.1;Name=XP_087097748.1;gbkey=CDS;gene=LOC101882419;protein_id=XP_087097748.1;transl_except=(pos:join(53263975%2C53265999..53266000)%2Caa:Other)\n"
        "NC_141022.1\tGnomon\tCDS\t53246650\t53246752\t.\t+\t0\tID=cds-XP_087097748.1;Parent=rna-XM_087241622.1;Name=XP_087097748.1;gbkey=CDS;gene=LOC101882419;protein_id=XP_087097748.1;transl_except=(pos:join(53263975%2C53265999..53266000)%2Caa:Other)\n"
        "NC_141022.1\tGnomon\tCDS\t53246866\t53247444\t.\t+\t2\tID=cds-XP_087097748.1;Parent=rna-XM_087241622.1;Name=XP_087097748.1;gbkey=CDS;gene=LOC101882419;protein_id=XP_087097748.1;transl_except=(pos:join(53263975%2C53265999..53266000)%2Caa:Other)\n"
        "NC_141022.1\tGnomon\tCDS\t53257673\t53257752\t.\t+\t2\tID=cds-XP_087097748.1;Parent=rna-XM_087241622.1;Name=XP_087097748.1;gbkey=CDS;gene=LOC101882419;protein_id=XP_087097748.1;transl_except=(pos:join(53263975%2C53265999..53266000)%2Caa:Other)\n"
        "NC_141022.1\tGnomon\tCDS\t53259929\t53260667\t.\t+\t0\tID=cds-XP_087097748.1;Parent=rna-XM_087241622.1;Name=XP_087097748.1;gbkey=CDS;gene=LOC101882419;protein_id=XP_087097748.1;transl_except=(pos:join(53263975%2C53265999..53266000)%2Caa:Other)\n"
        "NC_141022.1\tGnomon\tCDS\t53260670\t53263975\t.\t+\t2\tID=cds-XP_087097748.1;Parent=rna-XM_087241622.1;Name=XP_087097748.1;gbkey=CDS;gene=LOC101882419;protein_id=XP_087097748.1;transl_except=(pos:join(53263975%2C53265999..53266000)%2Caa:Other)\n"
        "NC_141022.1\tGnomon\tCDS\t53265999\t53266013\t.\t+\t2\tID=cds-XP_087097748.1;Parent=rna-XM_087241622.1;Name=XP_087097748.1;gbkey=CDS;gene=LOC101882419;protein_id=XP_087097748.1;transl_except=(pos:join(53263975%2C53265999..53266000)%2Caa:Other)\n"
        "NC_141022.1\tGnomon\tCDS\t53267737\t53268470\t.\t+\t2\tID=cds-XP_087097748.1;Parent=rna-XM_087241622.1;Name=XP_087097748.1;gbkey=CDS;gene=LOC101882419;protein_id=XP_087097748.1;transl_except=(pos:join(53263975%2C53265999..53266000)%2Caa:Other)\n"
    )
    cds = parse(content, "ncbi")["features"][0]["features"][0]["features"][0]
    assert cds["qualifiers"]["translation_exception"] == {
        "exceptions": [
            {
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": 53263974},
                    "end": {"type": "point", "position": 53266000},
                    "strand": 1,
                },
                "amino_acid": "X",
            }
        ]
    }


def test_sequence_not_retrieved_tark(monkeypatch):
    # The transcript is found in Tark, but the Ensembl sequence download fails.
    def seq_from_rest(*args):
        raise ConnectionError("(ensembl sequence)(Read timed out.)")

    monkeypatch.setattr(
        "mutalyzer_retriever.parsers.json_ensembl._seq_from_rest", seq_from_rest
    )

    with pytest.raises(NoReferenceError) as e:
        retrieve_model("ENST00000383925.1", "ensembl_tark", "json")
    assert e.value.uncertain_sources == ["ensembl_tark"]


def test_sequence_not_retrieved_fasta(monkeypatch):
    # The GFF3 is retrieved, but the FASTA download fails.
    def retrieve_raw(r_id, r_source=None, r_type=None, size_off=True, timeout=1):
        if r_type == "fasta":
            raise ConnectionError("(fasta)(Read timed out.)")
        return _retrieve_raw(r_id, r_source, r_type, size_off, timeout=timeout)

    monkeypatch.setattr("mutalyzer_retriever.retriever.retrieve_raw", retrieve_raw)

    with pytest.raises(NoReferenceError) as e:
        retrieve_model("NM_000077.4", "ncbi", "gff3")
    assert e.value.uncertain_sources == ["ncbi"]
