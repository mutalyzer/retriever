from __future__ import annotations
import json
import time
import requests
from urllib.parse import quote
from mutalyzer_retriever.request import Http400, request
from mutalyzer_retriever.configuration import settings
from mutalyzer_retriever.util import HUMAN_TAXON, DEFAULT_TIMEOUT


class BaseAPIClient:
    """Base class for API clients"""

    def __init__(self, base_url: str, timeout: int = DEFAULT_TIMEOUT):
        self.base_url = base_url
        self.timeout = timeout

    def make_request(self, url: str, params: dict | None = None):
        """Make HTTP request"""
        response = request(url=url, params=params, timeout=self.timeout)
        return json.loads(response)

class NCBIClient(BaseAPIClient):
    """Client for NCBI API operations"""

    def __init__(self, timeout: int = DEFAULT_TIMEOUT):
        base_url = settings.get("NCBI_DATASETS_API")
        super().__init__(base_url, timeout)

    def get_accession_dataset_report(self, accession: str):
        """Fetch dataset report for given accession"""
        url = f"{self.base_url}/gene/accession/{accession}/dataset_report"
        return self.make_request(url)

    def get_accession_product_report(self, accession: str):
        """Fetch product report for given accession"""
        url = f"{self.base_url}/gene/accession/{accession}/product_report"
        return self.make_request(url)

    def get_gene_id_dataset_report(self, gene_ids: list[str]):
        """Fetch dataset report for gene IDs"""
        gene_id_str = quote(",".join(map(str, gene_ids)))
        url = f"{self.base_url}/gene/id/{gene_id_str}/dataset_report"
        return self.make_request(url)

    def get_gene_id_product_report(self, gene_ids: list[str]):
        """Fetch product report for gene IDs"""
        gene_id_str = quote(",".join(map(str, gene_ids)))
        url = f"{self.base_url}/gene/id/{gene_id_str}/product_report"
        return self.make_request(url)

    def get_gene_symbol_dataset_report(self, gene_symbol: str, taxon_name: str = HUMAN_TAXON):
        """Fetch dataset report for gene symbol"""
        taxon_name_url_str = quote(str(taxon_name), safe="")
        url = f"{self.base_url}/gene/symbol/{gene_symbol}/taxon/{taxon_name_url_str}/dataset_report"
        return self.make_request(url)

    def get_gene_symbol_product_report(self, gene_symbol: str, taxon_name: str = HUMAN_TAXON):
        """Fetch product report for gene symbol"""
        taxon_name_url_str = quote(taxon_name, safe="")
        url = f"{self.base_url}/gene/symbol/{gene_symbol}/taxon/{taxon_name_url_str}/product_report"
        return self.make_request(url)

    def get_assembly_accession(self, accession: str):
        """Get assembly accession for sequence accession"""
        url = f"{self.base_url}/genome/sequence_accession/{accession}/sequence_assemblies"
        response = self.make_request(url)
        accessions = response.get("accessions")
        if isinstance(accessions, list) and accessions:
            return accessions[0]
        return None

    def get_genome_annotation_report(self, assembly_accession: str, locations: list[str]):
        """Get genome annotation report for assembly and locations"""
        url = f"{self.base_url}/genome/accession/{assembly_accession}/annotation_report"
        params = [("locations", loc) for loc in locations]
        return self.make_request(url, params)


class NCBIEutilsClient(BaseAPIClient):
    """Client for NCBI E-utilities API operations"""

    def __init__(self, timeout: int = DEFAULT_TIMEOUT, max_retries: int = 3):
        super().__init__("https://eutils.ncbi.nlm.nih.gov/entrez/eutils", timeout)
        self.max_retries = max_retries

    def make_request(self, url: str, params: dict | None = None):
        """Make HTTP request with retry logic for transient NCBI errors"""
        last_exception = None

        for attempt in range(self.max_retries):
            try:
                r = requests.get(url, params=params, timeout=self.timeout)
                r.raise_for_status()
                data = r.json()

                if "ERROR" in data:
                    raise RuntimeError(f"NCBI E-utilities error: {data['ERROR']}")

                return data

            except (json.JSONDecodeError, RuntimeError) as e:
                last_exception = e
                if attempt < self.max_retries - 1:
                    time.sleep(1 * (attempt + 1))  # Backoff: 1s, 2s, 3s
                    continue

        raise last_exception

    def elink(self, accession: str, dbfrom: str, db: str) -> dict:
        """Fetch links between NCBI databases"""
        url = f"{self.base_url}/elink.fcgi"
        params = {
            "id": accession,
            "dbfrom": dbfrom,
            "db": db,
            "retmode": "json",
            "api_key": settings.get("NCBI_API_KEY"),
            "email": settings.get("EMAIL")
        }
        return self.make_request(url, params)


class EnsemblClient(BaseAPIClient):
    """Client for Ensembl API operations"""

    def __init__(self, timeout: int = DEFAULT_TIMEOUT):
        super().__init__(settings.get("ENSEMBL_API"), timeout)

    def lookup_symbol(self, gene_symbol: str, taxon_name: str = "homo_sapiens"):
        """Lookup gene by symbol"""
        url = f"{self.base_url}/lookup/symbol/{taxon_name}/{gene_symbol}?content-type=application/json;expand=1"
        try:
            return self.make_request(url)
        except Http400:
            return {}

    def lookup_id(self, accession_base: str, expand: int = 1):
        """Lookup by Ensembl ID"""
        url = f"{self.base_url}/lookup/id/{accession_base}?content-type=application/json;expand={expand}"
        try:
            return self.make_request(url)
        except Http400:
            return {}
