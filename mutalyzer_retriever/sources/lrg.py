from urllib.error import URLError
from urllib.request import urlopen

from .. import configuration


class NoLrgSettings(Exception):
    """
    Raised when there are no LRG settings.
    """

    pass


class NoLrgUrlSet(Exception):
    """
    Raised when there is no LRG path specified in the settings.
    """

    pass


class LrgUrlAccessError(Exception):
    """
    Raised when there is some error when accessing the LRG url.
    """

    pass


class ReferenceToLong(Exception):
    """
    Raised when the reference length exceeds maximum size.
    """

    pass


class NotLrg(Exception):
    """
    Raised when the reference is not LRG.
    """

    pass


class NoSizeRetrieved(Exception):
    """
    Raised when the size of the LRG cannot be retrieved.
    """

    pass


def fetch_lrg(reference_id, size_on=True):
    """
    Fetch the LRG file content.

    :arg str reference_id: the name of the LRG file to fetch
    :arg bool size_on: flag for the maximum sequence length
    :returns: the file content or None when the file was not retrieved
    """
    settings = configuration.settings

    if settings.get("LRG_URL"):
        url = "{}/{}.xml".format(settings["LRG_URL"], reference_id)
    else:
        raise NoLrgUrlSet()

    try:
        handle = urlopen(url)
    except URLError:
        raise NameError

    info = handle.info()

    if info["Content-Type"] == "application/xml":
        if "Content-length" in info:
            if size_on:
                length = int(info["Content-Length"])
                if 512 > length or length > settings["MAX_FILE_SIZE"]:
                    handle.close()
                    raise ReferenceToLong(
                        "Filesize '{}' is not within the allowed boundaries "
                        "(512 < filesize < {} ) for {}.".format(
                            length, settings["MAX_FILE_SIZE"] // 1048576, reference_id,
                        )
                    )
        else:
            raise NoSizeRetrieved()
        raw_data = handle.read()
        handle.close()
        return raw_data.decode()
    else:
        handle.close()
        raise NotLrg()
