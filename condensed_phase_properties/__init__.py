"""Condensed phase properties namespace."""

from importlib.metadata import PackageNotFoundError, version

__author__ = "Anna Picha"
__email__ = "anna.picha@gmail.com"

# Used to automatically set version number from github actions
# as well as not break when being tested locally
try:
    __version__ = version(__package__)  # type: ignore
except PackageNotFoundError:  # pragma: no cover
    __version__ = "0.0.0"
