"""
Retriever configuration.
"""
import configparser
import os

DEFAULT_SETTINGS = {
    "NCBI_GFF3_URL": "https://eutils.ncbi.nlm.nih.gov/sviewer/viewer.cgi",
    "LRG_URL": "http://ftp.ebi.ac.uk/pub/databases/lrgex/",
    "MAX_FILE_SIZE": 10 * 1048576,
    "ENSEMBL_API": "https://rest.ensembl.org",
    "ENSEMBL_API_GRCH37": "https://grch37.rest.ensembl.org",
}


def setup_settings():
    """
    Setting up the configuration either from a file path (specified directly
    or via the MUTALYZER_RETRIEVER_SETTINGS environment variable) or from the
    default dictionary.

    :arg str configuration_path: Path towards a configuration file.
    :returns dict: Configuration dictionary.
    """
    if os.environ.get("MUTALYZER_SETTINGS"):
        configuration_path = os.environ["MUTALYZER_SETTINGS"]
        with open(configuration_path) as f:
            configuration_content = "[config]\n" + f.read()
        loaded_settings = configparser.ConfigParser()
        loaded_settings.optionxform = str
        loaded_settings.read_string(configuration_content)
        loaded_settings = {
            sect: dict(loaded_settings.items(sect))
            for sect in loaded_settings.sections()
        }["config"]
    else:
        loaded_settings = DEFAULT_SETTINGS

    if loaded_settings.get("MAX_FILE_SIZE") and isinstance(
        loaded_settings["MAX_FILE_SIZE"], str
    ):
        loaded_settings["MAX_FILE_SIZE"] = eval(loaded_settings["MAX_FILE_SIZE"])

    return loaded_settings


settings = setup_settings()
