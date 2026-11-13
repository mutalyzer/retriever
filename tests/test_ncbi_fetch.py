from http.client import IncompleteRead
from unittest.mock import Mock
from urllib.error import HTTPError

import pytest

from mutalyzer_retriever.sources import ncbi


@pytest.mark.parametrize("reference_type", [None, "gff3", "fasta", "genbank"])
@pytest.mark.parametrize(
    "reference_id,databases",
    [
        ("BFP69721.1", ["nuccore", "protein"]),
        ("NP_000068.1", ["protein", "nuccore"]),
    ],
)
def test_fetch_falls_back_on_lookup_miss(
    monkeypatch, reference_type, reference_id, databases
):
    fetch = Mock(side_effect=[NameError("not found"), "content"])
    result_type = reference_type or "gff3"
    monkeypatch.setattr(ncbi, f"fetch_{result_type}", fetch)

    assert ncbi.fetch(reference_id, reference_type, False, 42) == (
        "content",
        result_type,
    )
    extra_args = {"gff3": (42,), "fasta": (), "genbank": (False,)}[result_type]
    assert [call.args for call in fetch.call_args_list] == [
        (reference_id, db, *extra_args) for db in databases
    ]


@pytest.mark.parametrize("reference_type", ["gff3", "fasta", "genbank"])
def test_fetch_stops_after_both_databases_miss(monkeypatch, reference_type):
    errors = [NameError("nuccore miss"), NameError("protein miss")]
    fetch = Mock(side_effect=errors)
    monkeypatch.setattr(ncbi, f"fetch_{reference_type}", fetch)

    with pytest.raises(NameError) as exc:
        ncbi.fetch("BFP69721.1", reference_type)

    assert exc.value is errors[1]
    assert fetch.call_count == 2


@pytest.mark.parametrize("reference_type", ["gff3", "fasta", "genbank"])
def test_fetch_does_not_fall_back_on_connection_error(monkeypatch, reference_type):
    error = ConnectionError("service unavailable")
    fetch = Mock(side_effect=error)
    monkeypatch.setattr(ncbi, f"fetch_{reference_type}", fetch)

    with pytest.raises(ConnectionError) as exc:
        ncbi.fetch("BFP69721.1", reference_type)

    assert exc.value is error
    assert fetch.call_count == 1


@pytest.mark.parametrize(
    "message,exception",
    [
        ('Otherdb uid="2785095217" db="protein" term="2785095217"', NameError),
        ("Invalid uid INVALID_ACCESSION_999999999 at position= 0", NameError),
        ("Unable to obtain query #1", ConnectionError),
        ("Service unavailable", ConnectionError),
    ],
)
def test_summary_distinguishes_lookup_misses(monkeypatch, message, exception):
    handle = Mock()
    error = RuntimeError(message)
    monkeypatch.setattr(ncbi.Entrez, "esummary", Mock(return_value=handle))
    monkeypatch.setattr(ncbi.Entrez, "read", Mock(side_effect=error))

    with pytest.raises(exception) as exc:
        ncbi._fetch_reference_summary("BFP69721.1", "nuccore")

    assert exc.value.__cause__ is error
    handle.close.assert_called_once()


def test_genbank_falls_back_after_summary_reports_other_database(monkeypatch):
    handles = [Mock(), Mock()]
    summary = Mock(side_effect=handles)
    monkeypatch.setattr(ncbi.Entrez, "esummary", summary)
    monkeypatch.setattr(
        ncbi.Entrez,
        "read",
        Mock(
            side_effect=[
                RuntimeError(
                    'Otherdb uid="2785095217" db="protein" term="2785095217"'
                ),
                [{"AccessionVersion": "BFP69721.1", "Length": "100"}],
            ]
        ),
    )
    handle = Mock()
    handle.read.return_value = "genbank content"
    fetch = Mock(return_value=handle)
    monkeypatch.setattr(ncbi.Entrez, "efetch", fetch)

    assert ncbi.fetch("BFP69721.1", "genbank") == ("genbank content", "genbank")
    assert [call.kwargs["db"] for call in summary.call_args_list] == [
        "nuccore",
        "protein",
    ]
    fetch.assert_called_once_with(
        db="protein", id="BFP69721.1", rettype="gbwithparts", retmode="text"
    )
    for summary_handle in handles:
        summary_handle.close.assert_called_once()
    handle.close.assert_called_once()


def test_genbank_summary_service_failure_does_not_fall_back(monkeypatch):
    handle = Mock()
    summary = Mock(return_value=handle)
    monkeypatch.setattr(ncbi.Entrez, "esummary", summary)
    monkeypatch.setattr(
        ncbi.Entrez, "read", Mock(side_effect=RuntimeError("Service unavailable"))
    )
    fetch = Mock()
    monkeypatch.setattr(ncbi.Entrez, "efetch", fetch)

    with pytest.raises(ConnectionError):
        ncbi.fetch("BFP69721.1", "genbank")

    summary.assert_called_once_with(db="nuccore", id="BFP69721.1")
    fetch.assert_not_called()
    handle.close.assert_called_once()


@pytest.mark.parametrize("stage", ["open", "read"])
@pytest.mark.parametrize(
    "error",
    [
        HTTPError("url", 503, "Unavailable", None, None),
        OSError("connection lost"),
        IncompleteRead(b"partial"),
    ],
    # Explicit ids: on Python 3.9, pytest's id lookup on HTTPError raises KeyError.
    ids=["http_error", "os_error", "incomplete_read"],
)
def test_genbank_transport_failure_does_not_fall_back(monkeypatch, stage, error):
    summary = Mock(return_value={"length": 100})
    monkeypatch.setattr(ncbi, "_fetch_reference_summary", summary)
    handle = Mock()
    fetch = Mock(return_value=handle)
    if stage == "open":
        fetch.side_effect = error
    else:
        handle.read.side_effect = error
    monkeypatch.setattr(ncbi.Entrez, "efetch", fetch)

    with pytest.raises(ConnectionError) as exc:
        ncbi.fetch("NM_000077.4", "genbank")

    assert exc.value.__cause__ is error
    summary.assert_called_once_with("NM_000077.4", "nuccore")
    assert fetch.call_count == 1
    if stage == "read":
        handle.close.assert_called_once()


@pytest.mark.parametrize(
    "body,exception",
    [
        ("Error: F a i l e d  t o  u n d e r s t a n d  i d :  X Y Z \n\n\n", NameError),
        ("Error: Failed to understand id: XYZ", NameError),
        ("Error: Service unavailable", ConnectionError),
        ("<html>Service unavailable</html>", ConnectionError),
    ],
)
def test_fasta_rejects_non_fasta_body(monkeypatch, body, exception):
    handle = Mock()
    handle.read.return_value = body
    monkeypatch.setattr(ncbi.Entrez, "efetch", Mock(return_value=handle))

    with pytest.raises(exception):
        ncbi.fetch_fasta("XYZ", "nuccore")
