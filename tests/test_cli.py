from mutalyzer_retriever.cli import (
    _endpoint,
    _from_file,
    _parse_args,
    _related,
    _retrieve_assemblies,
    _retrieve_model,
    _retrieve_raw,
)


def test_retrieve_raw():
    command = "--id ref_id"
    args = _parse_args(command.split())

    assert args.id == "ref_id"
    assert args.parse is False

    assert args.command is None
    assert args.indent is None
    assert args.related is False
    assert args.sizeoff is False
    assert args.source is None
    assert args.timeout is None
    assert args.type is None

    assert _endpoint(args) == _retrieve_raw


def test_retrieve_model():
    command = "--id ref_id --parse"
    args = _parse_args(command.split())

    assert args.id == "ref_id"
    assert args.parse is True

    assert args.command is None
    assert args.indent is None
    assert args.related is False
    assert args.sizeoff is False
    assert args.source is None
    assert args.timeout is None
    assert args.type is None

    assert _endpoint(args) == _retrieve_model


def test_retrieve_model_indent():
    command = "--id ref_id --parse --indent 2"
    args = _parse_args(command.split())

    assert args.id == "ref_id"
    assert args.parse is True
    assert args.indent == 2

    assert args.command is None
    assert args.related is False
    assert args.sizeoff is False
    assert args.source is None
    assert args.timeout is None
    assert args.type is None

    assert _endpoint(args) == _retrieve_model


def test_related():
    command = "--id ref_id --related"
    args = _parse_args(command.split())

    assert args.id == "ref_id"
    assert args.related is True

    assert args.parse is False
    assert args.indent is None
    assert args.command is None
    assert args.sizeoff is False
    assert args.source is None
    assert args.timeout is None
    assert args.type is None

    assert _endpoint(args) == _related


def test_from_file():
    command = "from_file --paths some.gff3 some.fasta"
    args = _parse_args(command.split())

    assert args.command == "from_file"
    assert args.paths == ["some.gff3", "some.fasta"]

    assert args.id is None
    assert args.indent is None
    assert args.parse is False
    assert args.related is False
    assert args.sizeoff is False
    assert args.source is None
    assert args.timeout is None
    assert args.type is None

    assert _endpoint(args) == _from_file


def test_from_file_indent():
    command = "--indent 2 from_file --paths some.gff3 some.fasta"
    args = _parse_args(command.split())

    assert args.command == "from_file"
    assert args.paths == ["some.gff3", "some.fasta"]

    assert args.id is None
    assert args.indent == 2
    assert args.parse is False
    assert args.related is False
    assert args.sizeoff is False
    assert args.source is None
    assert args.timeout is None
    assert args.type is None

    assert _endpoint(args) == _from_file


def ncbi_assemblies():
    command = "ncbi_assemblies --input downloads --output models"
    args = _parse_args(command.split())

    assert args.command == "ncbi_assemblies"
    assert args.downloaded is False
    assert args.input == "downloads"
    assert args.output == "models"

    assert args.id is None
    assert args.indent is None
    assert args.parse is False
    assert args.related is False
    assert args.sizeoff is False
    assert args.source is None
    assert args.timeout is None
    assert args.type is None

    assert _endpoint(args) == _retrieve_assemblies


def test_ncbi_assemblies_downloaded():
    command = "ncbi_assemblies --input downloads --output models --downloaded"
    args = _parse_args(command.split())

    assert args.command == "ncbi_assemblies"
    assert args.downloaded is True
    assert args.input == "downloads"
    assert args.output == "models"

    assert args.id is None
    assert args.indent is None
    assert args.parse is False
    assert args.related is False
    assert args.sizeoff is False
    assert args.source is None
    assert args.timeout is None
    assert args.type is None

    assert _endpoint(args) == _retrieve_assemblies
