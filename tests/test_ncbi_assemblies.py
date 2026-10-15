import json
from pathlib import Path

from mutalyzer_retriever.sources.ncbi_assemblies import (
    get_annotation_models,
    load_metadata,
)

from .commons import _get_content, patch_retriever, references


def _expected_model(ref_id):
    p = Path(Path(__file__).parent) / "data" / str(ref_id + ".annotations")
    return json.load(p.open())


def test_ncbi_assemblies_general():
    metadata = load_metadata(
        str(Path(__file__).parent.joinpath("data/annotation_releases/"))
    )
    models = get_annotation_models(metadata)
    for ref_id in ["NT_113901.1", "NT_167208.1", "NT_167222.1"]:
        assert models["GRCh37"][ref_id] == _expected_model(ref_id)
    assert models["GRCh38"]["NT_167214.1"] == _expected_model("NT_167214.1")
    # tRNA genes in 26.70-26.82 Mb only: two TRNAM-CAU genes (one id, two
    # loci) renamed after 20190607, both kept as history.
    assert models["GRCh38"]["NC_000006.12"] == _expected_model("NC_000006.12")
    # RNU2-2 only: its snRNA NR_199791.1 appears in 20250801, as an ncRNA.
    assert models["GRCh38"]["NC_000011.10"] == _expected_model("NC_000011.10")
    # RMRP only: RNase_MRP_RNA NR_003051.3, replaced by NR_003051.4 in 20250801.
    assert models["GRCh38"]["NC_000009.12"] == _expected_model("NC_000009.12")
    # One gene per ncRNA type not found elsewhere (snoRNA, scaRNA, Y_RNA, ...);
    # NC_000007.14 also has the CFAP144P1 transcribed pseudogene (NR_028347.1).
    for ref_id in ["NC_000003.12", "NC_000014.9", "NC_000005.10", "NC_000007.14"]:
        assert models["GRCh38"][ref_id] == _expected_model(ref_id)


def test_ncbi_assemblies_restrict_ref_id_grch38():
    metadata = load_metadata(
        str(Path(__file__).parent.joinpath("data/annotation_releases/")), "GRCh38"
    )
    ref_id = "NT_167214.1"
    models = get_annotation_models(metadata, ref_id_start=ref_id)
    assert models["GRCh38"][ref_id] == _expected_model(ref_id)


def test_ncbi_assemblies_restrict_ref_id():
    metadata = load_metadata(
        str(Path(__file__).parent.joinpath("data/annotation_releases/"))
    )
    ref_id = "NT_113901.1"
    models = get_annotation_models(metadata, ref_id_start=ref_id)
    assert models["GRCh37"][ref_id] == _expected_model(ref_id)


def test_ncbi_assemblies_restrict_ref_id_and_assembly_id():
    assembly_id = "GRCh37"
    metadata = load_metadata(
        str(Path(__file__).parent.joinpath("data/annotation_releases/")), assembly_id
    )
    ref_id = "NT_113901.1"
    models = get_annotation_models(metadata, ref_id_start=ref_id)
    assert models[assembly_id][ref_id] == _expected_model(ref_id)
