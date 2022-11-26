from datetime import datetime
from importlib.metadata import version, PackageNotFoundError


__author__ = ("Haibao Tang", "Vivek Krishnakumar", "Jingping Li")
__copyright__ = "Copyright (c) 2010-{}, Haibao Tang".format(datetime.now().year)
__email__ = "tanghaibao@gmail.com"
__license__ = "BSD"
__status__ = "Development"

try:
    __version__ = version("jcvi")
except PackageNotFoundError:
    # package is not installed
    pass
