__author__ = ("Haibao Tang", "Vivek Krishnakumar", "Jingping Li")
__copyright__ = "Copyright (c) 2010-2019, Haibao Tang"
__email__ = "tanghaibao@gmail.com"
__license__ = "BSD"
__status__ = "Development"

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions
