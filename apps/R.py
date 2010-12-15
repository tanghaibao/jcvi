"""
uses R for statistics and graphics
"""

import os
import logging
import tempfile
from string import Template

from jcvi.apps.base import sh


class RTemplate (object):
    """
    Creates a R script and runs it
    """
    def __init__(self, template, parameters):

        self.template = Template(template)
        self.parameters = parameters

    def run(self, clean=True):
        """
        Create a temporary file and run it
        """
        template = self.template
        parameters = self.parameters
        # write to a temporary R script
        fd, path = tempfile.mkstemp()
        fw = os.fdopen(fd, "w")

        fw.write(template.substitute(**parameters))
        fw.close()

        sh("Rscript %s" % path)
        if clean:
            os.remove(path)


def main():
    template = """
    png("$pngfile")
    data(cars)
    plot(cars, xlab="Speed (mph)", ylab="Stopping distance (ft)", \
            las=1, main="cars data")
    dev.off()
    """
    logging.basicConfig(level=logging.DEBUG)
    rtemplate = RTemplate(template, dict(pngfile="t.png"))
    rtemplate.run()


if __name__ == '__main__':
    main()
