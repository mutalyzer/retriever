import pytest
import requests

from mutalyzer_retriever import parser
from mutalyzer_retriever.parsers.json_ensembl import _seq_from_rest


def test_json_parse_does_not_fetch_sequence(monkeypatch):
    def seq_from_rest(*args, **kwargs):
        raise AssertionError("json_ensembl.parse() should not fetch sequence")

    monkeypatch.setattr(
        "mutalyzer_retriever.parsers.json_ensembl._seq_from_rest", seq_from_rest
    )

    tark_results = {
        "results": [
            {
                "stable_id": "ENST00000000000",
                "stable_id_version": 1,
                "assembly": "GRCh38",
                "loc_region": "1",
                "loc_start": 10,
                "loc_end": 20,
                "loc_strand": 1,
                "exons": [],
                "translations": [],
                "genes": [
                    {
                        "stable_id": "ENSG00000000000",
                        "stable_id_version": 1,
                        "assembly": "GRCh38",
                        "loc_start": 10,
                        "loc_end": 20,
                        "loc_strand": 1,
                        "name": "GENE",
                    }
                ],
                "biotype": "protein_coding",
            }
        ]
    }

    model = parser.parse(tark_results, "json")

    assert model["type"] == "record"


def _unavailable(*args, **kwargs):
    response = requests.Response()
    response.status_code = 503
    return response


def _timeout(*args, **kwargs):
    raise requests.exceptions.ReadTimeout("Read timed out.")


@pytest.mark.parametrize("get", [_unavailable, _timeout])
def test_seq_from_rest_connection_error(get, monkeypatch):
    monkeypatch.setattr("mutalyzer_retriever.parsers.json_ensembl.requests.get", get)

    with pytest.raises(ConnectionError):
        _seq_from_rest("GRCh38", "1", 1, 10, 20, timeout=1)
