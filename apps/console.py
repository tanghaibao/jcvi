#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import sys
import time
import threading
import unicodedata


def tabular(listOfStuff, rulersize=80):
    """
    Print a tabular output, with horizontal separators
    """
    table_edge = "=" * rulersize + "\n"
    table_sep = "-" * rulersize + "\n"
    contents = table_sep.join(str(x) + "\n" for x in listOfStuff)
    return "".join((table_edge, contents, table_edge))

"""
An ASCII text progress bar. See __main__ for command line use (using \r to
move the cursor back to the start of the current line is the key, on
terminals that do not support this functionality the progress bar will
not work as well).

<http://code.google.com/p/python-progressbar/>
"""

from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
         RotatingMarker, ReverseBar, SimpleProgress


class SpinCursor(threading.Thread):
    """
    A console spin cursor class based on:
    <http://code.activestate.com/recipes/534142-spin-cursor/>
    """
    def __init__(self, msg='', maxspin=0, minspin=10, speed=5):
        # Count of a spin
        self.count = 0
        self.out = sys.stdout
        self.flag = False
        self.max = maxspin
        self.min = minspin
        # Any message to print first ?
        self.msg = msg
        # Complete printed string
        self.string = ''
        # Speed is given as number of spins a second
        # Use it to calculate spin wait time
        self.waittime = 1.0 / float(speed * 4)
        if os.name == 'posix':
            self.spinchars = (unicodedata.lookup('FIGURE DASH'), \
                    u'\\ ', u'| ', u'/ ')
        else:
            # The unicode dash character does not show
            # up properly in Windows console.
            self.spinchars = (u'-', u'\\ ', u'| ', u'/ ')
        threading.Thread.__init__(self, None, None, "Spin Thread")

    def spin(self):
        """ Perform a single spin """

        for x in self.spinchars:
            self.string = self.msg + "...\t" + x + "\r"
            self.out.write(self.string.encode('utf-8'))
            self.out.flush()
            time.sleep(self.waittime)

    def run(self):

        while (not self.flag) and \
              ((self.count < self.min) or (self.count < self.max)):
            self.spin()
            self.count += 1

        # Clean up display...
        self.out.write(" " * (len(self.string) + 1))

    def stop(self):
        self.flag = True


"""
ANSI Color formatting for output in terminal based on termcolor
<http://pypi.python.org/pypi/termcolor>

Copyright (C) 2008-2009 Konstantin Lepa <konstantin.lepa@gmail.com>.
"""
ATTRIBUTES = dict(zip(
    ('bold', 'dark', '', 'underline', 'blink', '', 'reverse', 'concealed'),
    range(1, 9))
)
del ATTRIBUTES['']

HIGHLIGHTS = dict(zip(
    ('on_grey', 'on_red', 'on_green', 'on_yellow', 'on_blue', 'on_magenta',
        'on_cyan', 'on_white'),
    range(40, 48))
)
ATTRIBUTES.update(HIGHLIGHTS)

COLORS = dict(zip(
    ('grey', 'red', 'green', 'yellow', 'blue', 'magenta', 'cyan', 'white', ),
    range(30, 38))
)
ATTRIBUTES.update(COLORS)

RESET = '\033[0m'


class ColoredText:

    def __init__(self, text, attrs=None):
        self.text = text
        attrs = [x.strip() for x in attrs.strip().split('|')]
        self.attrs = [x for x in attrs if x in ATTRIBUTES]

    def __str__(self):
        """Colorize text.

        Available text colors:
            red, green, yellow, blue, magenta, cyan, white.

        Available text highlights:
            on_red, on_green, on_yellow, on_blue, on_magenta, on_cyan, on_white.

        Available attributes:
            bold, dark, underline, blink, reverse, concealed.

        Example:
            ColoredText('Hello, World!', 'red|on_grey|blue|blink')
            ColoredText('Hello, World!', 'green')
        """
        ctext = None
        if os.getenv('ANSI_COLORS_DISABLED') is None:
            fmt_str = '\033[%dm%s'

            if self.attrs:
                for attr in self.attrs:
                    ctext = fmt_str % (ATTRIBUTES[attr], self.text)
                ctext += RESET

        return ctext or self.text

    __repr__ = __str__


def print_red(s):
    print(ColoredText(s, "red"))


def print_green(s):
    print(ColoredText(s, "green"))


if __name__ == '__main__':

    """
    test cases for progressbar taken from:
    <http://python-progressbar.googlecode.com/hg/progressbar/examples.py>
    """

    def example0():
        pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=300).start()
        for i in range(300):
            time.sleep(0.01)
            pbar.update(i + 1)
        pbar.finish()
        sys.stdout.write('\n')

    def example1():
        widgets = ['Test: ', Percentage(), ' ', Bar(marker=RotatingMarker()),
                   ' ', ETA(), ' ', FileTransferSpeed()]
        pbar = ProgressBar(widgets=widgets, maxval=10000000).start()
        for i in range(1000000):
            # do something
            pbar.update(10 * i + 1)
        pbar.finish()
        sys.stdout.write('\n')

    def example2():
        class CrazyFileTransferSpeed(FileTransferSpeed):
            "It's bigger between 45 and 80 percent"
            def update(self, pbar):
                if 45 < pbar.percentage() < 80:
                    return 'Bigger Now ' + FileTransferSpeed.update(self, pbar)
                else:
                    return FileTransferSpeed.update(self, pbar)

        widgets = [CrazyFileTransferSpeed(), ' <<<', Bar(), '>>> ',
                   Percentage(), ' ', ETA()]
        pbar = ProgressBar(widgets=widgets, maxval=10000000)
        # maybe do something
        pbar.start()
        for i in range(2000000):
            # do something
            pbar.update(5 * i + 1)
        pbar.finish()
        sys.stdout.write('\n')

    def example3():
        widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
        pbar = ProgressBar(widgets=widgets, maxval=10000000).start()
        for i in range(1000000):
            # do something
            pbar.update(10 * i + 1)
        pbar.finish()
        sys.stdout.write('\n')

    def example4():
        widgets = ['Test: ', Percentage(), ' ',
                   Bar(marker='0', left='[', right=']'),
                   ' ', ETA(), ' ', FileTransferSpeed()]
        pbar = ProgressBar(widgets=widgets, maxval=500)
        pbar.start()
        for i in range(100, 500 + 1, 50):
            time.sleep(0.2)
            pbar.update(i)
        pbar.finish()
        sys.stdout.write('\n')

    def example5():
        pbar = ProgressBar(widgets=[SimpleProgress()], maxval=17).start()
        for i in range(17):
            time.sleep(0.2)
            pbar.update(i + 1)
        pbar.finish()
        sys.stdout.write('\n')

    for f in (example0, example1, example2, example3, example4, example5):
        f()

    # test spin cursor
    spin = SpinCursor(msg="Spinning...", minspin=5, speed=5)
    spin.start()
    spin.join()
    print

    # testing ANSI colors and text
    print('Current terminal type: %s' % os.getenv('TERM'))
    print('Test basic colors:')
    print(ColoredText('Grey color', 'grey'))
    print(ColoredText('Red color', 'red'))
    print(ColoredText('Green color', 'green'))
    print(ColoredText('Yellow color', 'yellow'))
    print(ColoredText('Blue color', 'blue'))
    print(ColoredText('Magenta color', 'magenta'))
    print(ColoredText('Cyan color', 'cyan'))
    print(ColoredText('White color', 'white'))
    print('-' * 78)

    print('Test highlights:')
    print(ColoredText('On grey color', 'on_grey'))
    print(ColoredText('On red color', 'on_red'))
    print(ColoredText('On green color', 'on_green'))
    print(ColoredText('On yellow color', 'on_yellow'))
    print(ColoredText('On blue color', 'on_blue'))
    print(ColoredText('On magenta color', 'on_magenta'))
    print(ColoredText('On cyan color', 'on_cyan'))
    print(ColoredText('On white color', 'grey|on_white'))
    print('-' * 78)

    print('Test attributes:')
    print(ColoredText('Bold grey color', 'grey|bold'))
    print(ColoredText('Dark red color', 'red|dark'))
    print(ColoredText('Underline green color', 'green|underline'))
    print(ColoredText('Blink yellow color', 'yellow|blink'))
    print(ColoredText('Reversed blue color', 'blue|reverse'))
    print(ColoredText('Concealed Magenta color', 'magenta|concealed'))
    print(ColoredText('Bold underline reverse cyan color',
        'cyan|bold|underline|reverse'))
    print(ColoredText('Dark blink concealed white color',
        'white|dark|blink|concealed'))
    print('-' * 78)

    print('Test mixing:')
    print(ColoredText('Underline red on grey color',
        'red|on_grey|underline'))
    print(ColoredText('Reversed green on red color',
        'green|on_red|reverse'))
