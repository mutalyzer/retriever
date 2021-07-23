"""
Retriever configuration.
"""
import configparser
import os

settings = {
    "EMAIL": "make_sure_to_configure_this@email.com",
    "NCBI_API_KEY": None,
    "LRG_URL": "http://ftp.ebi.ac.uk/pub/databases/lrgex/",
    "MAX_FILE_SIZE": 10 * 1048576,
}


def setup_settings(configuration_path=None):
    """
    Setting up the configuration either from a file path (specified directly
    or via the MUTALYZER_RETRIEVER_SETTINGS environment variable) or from the
    default dictionary.

    :arg str configuration_path: Path towards a configuration file.
    :returns dict: Configuration dictionary.
    """
    if configuration_path is None and os.environ.get("MUTALYZER_RETRIEVER_SETTINGS"):
        configuration_path = os.environ["MUTALYZER_RETRIEVER_SETTINGS"]
    if configuration_path:
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
        loaded_settings = settings

    if loaded_settings.get("MAX_FILE_SIZE") and isinstance(
        loaded_settings["MAX_FILE_SIZE"], str
    ):
        loaded_settings["MAX_FILE_SIZE"] = eval(loaded_settings["MAX_FILE_SIZE"])

    return loaded_settings
