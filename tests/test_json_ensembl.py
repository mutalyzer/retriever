from mutalyzer_retriever import parser


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
