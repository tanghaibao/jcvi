import logging


class BaseFile (object):
    def __init__(self, filename):

        self.filename = filename
        logging.debug("Load file %s" % filename)


class LineFile (BaseFile, list):
    """
    Generic file parser for line-based files
    """
    def __init__(self, filename):

        super(LineFile, self).__init__(filename)
