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
