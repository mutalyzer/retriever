import json
import requests
from ..configuration import settings
from ..request import Http400, RequestErrors, request
from ..util import f_e


def fetch_json(feature_id, api_base, timeout=1):
    url = f"{api_base}/lookup/id/{feature_id}"
    params = {"feature": ["gene", "transcript", "cds"], "expand": 1}
    headers = {"Content-Type": "application/json"}
    try:
        response = request(url, params, headers, timeout=timeout)
    except RequestErrors as e:
        raise ConnectionError(f"(json) {str(e)}")
    except Http400 as e:
        response_json = e.response.json()
        if response_json and response_json.get("error") == "ID '{}' not found".format(
            feature_id
        ):
            raise NameError(f"(json) {str(e)}")
        else:
            raise e
    else:
        return response


def fetch_fasta(feature_id, api_base, timeout=1):
    url = f"{api_base}/sequence/id/{feature_id}"
    params = {"format": "fasta", "type": "genomic"}
    headers = {"Content-Type": "text/x-fasta"}

    try:
        response = request(url, params, headers, timeout=timeout)
    except RequestErrors as e:
        raise ConnectionError(f_e("gff3", e))
    except Http400 as e:
        response_json = e.response.json()
        if response_json and response_json.get("error") == "ID '{}' not found".format(
            feature_id
        ):
            raise NameError(f_e("fasta", e, response_json.get("error")))
        else:
            raise e
    else:
        return response


def fetch_gff3(feature_id, api_base, timeout=1):
    url = f"{api_base}/overlap/id/{feature_id}"
    params = {"feature": ["gene", "transcript", "cds", "exon"]}
    headers = {"Content-Type": "text/x-gff3"}

    try:
        response = request(url, params, headers, timeout=timeout)
    except RequestErrors as e:
        raise ConnectionError(f_e("gff3", e))
    except Http400 as e:
        response_json = e.response.json()
        if response_json and response_json.get("error") == "ID '{}' not found".format(
            feature_id
        ):
            raise NameError(f_e("gff3", e, response_json.get("error")))
        else:
            raise e
    else:
        return response


def _get_tark_versions(reference_id, api_base, timeout=4):
    endpoint = "transcript"
    
    
    params = {"stable_id": reference_id}
    tark_req = json.loads(request(url=f"{api_base}/{endpoint}", params=params))
    tark_versions_38 = []
    tark_versions_37 = []
    if tark_req["results"]:
        for r in tark_req["results"]:
            if r["assembly"] == "GRCh37":
                tark_versions_37.append(int(r["stable_id_version"]))
            elif r["assembly"] == "GRCh38":
                tark_versions_38.append(int(r["stable_id_version"]))
            
    return tark_versions_38, tark_versions_37

def _get_most_recent_version(reference_id, api_base, timeout=1):
    return int(_get_reference_information(reference_id, api_base, timeout)["version"])


def _get_reference_information(reference_id, api_base, timeout=1):
    url = f"{api_base}/lookup/id/{reference_id}"
    headers = {"Content-Type": "application/json"}
    return json.loads(request(url, headers=headers, timeout=timeout))


def _get_id_and_version(reference_id):
    r_id = None
    r_version = None
    if reference_id.startswith("ENS"):
        if (
            "." in reference_id
            and len(reference_id.split(".")) == 2
            and reference_id.split(".")[1].isdigit()
        ):
            r_id, r_version = reference_id.split(".")
            r_version = int(r_version)
        else:
            r_id = reference_id
    return r_id, r_version


def _in_grch37(r_id, r_version, r_info, timeout):
    api_base = settings.get("ENSEMBL_API_GRCH37")
    if r_info["species"] == "homo_sapiens" and int(r_info["version"]) > r_version:
        grch37_version = _get_most_recent_version(r_id, api_base, timeout)
        if grch37_version and grch37_version == r_version:
            return True
    return False




def fetch_tark(reference_id, api_base, assembly= "GRCh38"):
    endpoint = "transcript"
    reference_id, reference_version =  _get_id_and_version(reference_id)
    
    if reference_id is None:
        raise NameError
    if reference_version:
        params = {"stable_id": reference_id,
                "assembly_name": assembly,
                "stable_id_version":reference_version,
                "expand": "sequence, translations, genes, exons"}
    else:
        params = {"stable_id": reference_id,
                "assembly_name": assembly,
                "expand": "sequence, translations, genes, exons"}
    req = requests.request(method="get", url=f"{api_base}/{endpoint}", params=params)
 
    return req.json()




def get_api_base(r_id, r_version, transcript = False):
    rest_version_38 = _get_most_recent_version(r_id,settings.get("ENSEMBL_API"))
    rest_version_37 = _get_most_recent_version(r_id,settings.get("ENSEMBL_API_GRCH37"))
    if not transcript:
        if r_version in [None, rest_version_38]:
            return settings.get("ENSEMBL_API")
        elif r_version == rest_version_37:
            return settings.get("ENSEMBL_API_GRCH37") 

        
    if transcript:
        tark_versions_38, tark_versions_37 = _get_tark_versions(r_id,settings.get("ENSEMBL_TARK_API"))
        if not r_version:
            return settings.get("ENSEMBL_TARK_API"), "GRCh38"
        elif r_version in tark_versions_38:
            return settings.get("ENSEMBL_TARK_API"), "GRCh38"
        elif r_version in tark_versions_37:
            return settings.get("ENSEMBL_TARK_API"), "GRCh37"
    raise ValueError(f"Cannot fetch {r_id} with version {r_version} from Ensembl")



def fetch(reference_id, reference_type=None, timeout=1):
    r_id, r_version = _get_id_and_version(reference_id)
    if r_id is None:
        raise NameError

    if reference_type == "gff3":
        api_base = get_api_base(r_id, r_version)
        return fetch_gff3(r_id, api_base, timeout), "gff3"
    elif reference_type == "fasta":
        api_base = get_api_base(r_id, r_version)
        return fetch_fasta(r_id,api_base), "fasta"
    elif reference_type in [None,"json"] and "ENST" not in r_id:
        api_base = get_api_base(r_id, r_version)
        return fetch_json(r_id,api_base), "json"  
    elif reference_type in [None,"json"] and "ENST" in r_id:
        api_base, assembly = get_api_base(r_id, r_version, transcript=True)
        print(api_base,assembly)
        return fetch_tark(r_id,api_base,assembly), "json"
   
    elif reference_type == "genbank":
        return None, "genbank"

    raise ValueError(
        "Ensembl fetch does not support '{}' reference type.".format(reference_type)
    )

