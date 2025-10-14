from enum import Enum

class StrEnum(str, Enum):
    """
    Enum where members are also (and must be) strings
    """
    pass


# Constants
DEFAULT_TIMEOUT: int = 10
HUMAN_TAXON: str = "HOMO SAPIENS"
EMPTY_VALUES = (None, "", [])

class DataSource(StrEnum):
    """Data source enumeration"""
    ENSEMBL = "ENSEMBL"
    NCBI = "NCBI"
    OTHER = "OTHER"

class MoleculeType(StrEnum):
    """Molecule type enumeration"""
    DNA = "dna"
    RNA = "rna"
    PROTEIN = "protein"
    UNKNOWN = "unknown"
    OTHER = "other"


def f_e(source, e=None, extra=None):
    output = f"({source})"
    if e is not None:
        output += f"({str(e)})"
    if extra is not None:
        output += f" ({extra})"
    return output


def make_location(start, end=None, strand=None):
    if end is not None:
        location = {
            "type": "range",
            "start": {"type": "point", "position": int(start)},
            "end": {"type": "point", "position": int(end)},
        }
    else:
        location = {"type": "point", "position": int(start)}
    if strand is not None:
        location["strand"] = int(strand)
    return location
