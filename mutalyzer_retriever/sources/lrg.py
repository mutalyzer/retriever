import socket
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from ..configuration import DEFAULT_TIMEOUT, settings
from ..util import f_e


class NoLrgSettings(Exception):
    """
    Raised when there are no LRG settings.
    """


class NoLrgUrlSet(Exception):
    """
    Raised when there is no LRG path specified in the settings.
    """


class LrgUrlAccessError(Exception):
    """
    Raised when there is some error when accessing the LRG url.
    """


class ReferenceToLong(Exception):
    """
    Raised when the reference length exceeds maximum size.
    """


class NotLrg(Exception):
    """
    Raised when the reference is not LRG.
    """


class NoSizeRetrieved(Exception):
    """
    Raised when the size of the LRG cannot be retrieved.
    """


def fetch_lrg(reference_id, size_on=True, timeout=DEFAULT_TIMEOUT):
    """
    Fetch the LRG file content.

    :arg str reference_id: the name of the LRG file to fetch
    :arg bool size_on: flag for the maximum sequence length
    :returns: the file content or None when the file was not retrieved
    """
    if settings.get("LRG_URL"):
        url = f"{settings['LRG_URL']}/{reference_id}.xml"
    else:
        raise NoLrgUrlSet()

    try:
        handle = urlopen(url, timeout=timeout)
    except socket.timeout as e:
        raise ConnectionError(
            f_e("socket.timeout", e, f"timeout set to {timeout} s")
        ) from e
    except HTTPError as e:
        if e.code == 404:
            raise NameError(f_e("HTTPError", e)) from e
        raise ConnectionError(f_e("HTTPError", e)) from e
    except URLError as e:
        if hasattr(e, "reason"):
            raise ConnectionError(f_e("URLError", e)) from e

    info = handle.info()

    if info["Content-Type"].strip().endswith("/xml"):
        if "Content-length" in info:
            if size_on:
                length = int(info["Content-Length"])
                if 512 > length or length > settings["MAX_FILE_SIZE"]:
                    handle.close()
                    raise ReferenceToLong(
                        f"Filesize '{length}' is not within the allowed boundaries "
                        f"(512 < filesize < {settings['MAX_FILE_SIZE'] // 1048576} ) "
                        f"for {reference_id}."
                    )
        else:
            raise NoSizeRetrieved()
        try:
            raw_data = handle.read()
        except socket.timeout as e:
            raise ConnectionError(
                f_e("socket.timeout", e, f"timeout set to {timeout} s")
            ) from e
        handle.close()
        return raw_data.decode()
    handle.close()
    raise NotLrg()
