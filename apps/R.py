"""
uses R for statistics and graphics
"""

import os
import os.path as op
import tempfile

from string import Template

from jcvi.formats.base import must_open
from jcvi.apps.base import sh, debug
debug()


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
        fw = must_open("tmp", "w")
        path = fw.name

        fw.write(template.safe_substitute(**parameters))
        fw.close()

        sh("Rscript %s" % path)
        if clean:
            os.remove(path)
            # I have no idea why using ggsave, there is one extra image
            # generated, but here I remove it
            rplotspdf = "Rplots.pdf"
            if op.exists(rplotspdf):
                os.remove(rplotspdf)


def main():
    template = """
    png("$pngfile")
    data(cars)
    plot(cars, xlab="Speed (mph)", ylab="Stopping distance (ft)", \
            las=1, main="cars data")
    dev.off()
    """
    rtemplate = RTemplate(template, dict(pngfile="t.png"))
    rtemplate.run()


if __name__ == '__main__':
    main()
