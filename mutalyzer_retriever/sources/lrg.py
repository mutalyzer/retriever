import socket
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from ..configuration import settings
from ..util import f_e


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


def fetch_lrg(reference_id, size_on=True, timeout=1):
    """
    Fetch the LRG file content.

    :arg str reference_id: the name of the LRG file to fetch
    :arg bool size_on: flag for the maximum sequence length
    :returns: the file content or None when the file was not retrieved
    """
    if settings.get("LRG_URL"):
        url = "{}/{}.xml".format(settings["LRG_URL"], reference_id)
    else:
        raise NoLrgUrlSet()

    try:
        handle = urlopen(url, timeout=timeout)
    except socket.timeout as e:
        raise ConnectionError(f_e("socket.timeout", e, f"timeout set to {timeout} s"))
    except HTTPError as e:
        if e.code == 404:
            raise NameError(f_e("HTTPError", e))
        else:
            raise ConnectionError(f_e("HTTPError", e))
    except URLError as e:
        if hasattr(e, "reason"):
            raise ConnectionError(f_e("URLError", e))

    info = handle.info()

    if info["Content-Type"].strip().endswith("/xml"):
        if "Content-length" in info:
            if size_on:
                length = int(info["Content-Length"])
                if 512 > length or length > settings["MAX_FILE_SIZE"]:
                    handle.close()
                    raise ReferenceToLong(
                        "Filesize '{}' is not within the allowed boundaries "
                        "(512 < filesize < {} ) for {}.".format(
                            length,
                            settings["MAX_FILE_SIZE"] // 1048576,
                            reference_id,
                        )
                    )
        else:
            raise NoSizeRetrieved()
        try:
            raw_data = handle.read()
        except socket.timeout as e:
            raise ConnectionError(f_e("socket.timeout", e, f"timeout set to {timeout} s"))
        handle.close()
        return raw_data.decode()
    else:
        handle.close()
        raise NotLrg()
