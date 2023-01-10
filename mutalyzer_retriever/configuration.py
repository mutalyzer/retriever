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
    Setting up the configuration from the default dictionary above or (/ond
    updated) from a file path specified via the MUTALYZER_SETTINGS
    environment variable.

    :returns dict: Configuration dictionary.
    """
    settings = DEFAULT_SETTINGS
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
        for k in loaded_settings:
            if loaded_settings[k] in {"yes", "true", "1"}:
                loaded_settings[k] = True
            elif loaded_settings[k] in {"no", "false", "0"}:
                loaded_settings[k] = False
            elif loaded_settings[k].isnumeric():
                loaded_settings[k] = int(loaded_settings[k])
        settings.update(loaded_settings)

    return settings


settings = setup_settings()


def cache_dir():
    return settings.get("MUTALYZER_CACHE_DIR")


def cache_url():
    return settings.get("MUTALYZER_API_URL")


def lru_cache_maxsize():
    return settings.get("MUTALYZER_LRU_CACHE_MAXSIZE")


def cache_add():
    return settings.get("MUTALYZER_FILE_CACHE_ADD")
