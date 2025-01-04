from datetime import datetime

__author__ = (
    "Haibao Tang",
    "Vivek Krishnakumar",
    "Adam Taranto",
    "Xingtan Zhang",
    "Won Cheol Yim",
)
__copyright__ = f"Copyright (c) 2010-{datetime.now().year}, Haibao Tang"
__email__ = "tanghaibao@gmail.com"
__license__ = "BSD"
__status__ = "Development"

try:
    from ._version import __version__  # noqa
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "Failed to find (autogenerated) _version.py. "
        "This might be because you are installing from GitHub's tarballs, "
        "use the PyPI ones."
    ) from exc
