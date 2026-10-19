from importlib.metadata import metadata

from .retriever import retrieve_model, retrieve_raw


def _get_metadata(name):
    meta = metadata(__package__)
    return meta.get(name, "")


_copyright_notice = (
    f"Copyright (c) {_get_metadata('Author')} <{_get_metadata('Author-email')}>"
)

usage = [_get_metadata("Summary"), _copyright_notice]


def doc_split(func):
    return func.__doc__.split("\n\n")[0]


def version(name):
    return (
        f"{_get_metadata('Name')} version {_get_metadata('Version')}\n\n"
        f"{_copyright_notice}\nHomepage: {_get_metadata('Home-page')}"
    )
