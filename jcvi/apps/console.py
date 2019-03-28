#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import print_function
import os
import sys
import time
import threading
import unicodedata


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
    ["0;%d" % x for x in range(30, 38)])
)
ATTRIBUTES.update(COLORS)

DARKCOLORS = dict(zip(
    ('black', 'darkred', 'darkgreen', 'darkyellow', 'darkblue', 'darkmagenta', \
     'darkcyan', 'silver'),
    ["1;%d" % x for x in range(30, 38)])
)
ATTRIBUTES.update(DARKCOLORS)

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
            fmt_str = '\033[%sm%s'

            if self.attrs:
                for attr in self.attrs:
                    ctext = fmt_str % (ATTRIBUTES[attr], self.text)
                ctext += RESET

        return ctext or self.text

    __repr__ = __str__


grey = gray = lambda s: str(ColoredText(s, "grey"))
red = lambda s: str(ColoredText(s, "red"))
green = lambda s: str(ColoredText(s, "green"))
yellow = lambda s: str(ColoredText(s, "yellow"))
blue = lambda s: str(ColoredText(s, "blue"))
magenta = lambda s: str(ColoredText(s, "magenta"))
cyan = lambda s: str(ColoredText(s, "cyan"))
white = lambda s: str(ColoredText(s, "white"))
dark = lambda s: str(ColoredText(s, "dark"))


def test():
    # test spin cursor
    spin = SpinCursor(msg="Spinning...", minspin=5, speed=5)
    spin.start()
    spin.join()
    print()

    # test ANSI colors and text
    print('Current terminal type: %s' % os.getenv('TERM'))
    print('Test basic colors:')
    for c in COLORS.keys():
        print(ColoredText("{0} color".format(c.capitalize()), c))
    print('-' * 78)

    for c in DARKCOLORS.keys():
        print(ColoredText("{0} color".format(c.capitalize()), c))
    print('-' * 78)

    print('Test highlights:')
    for c in HIGHLIGHTS.keys():
        print(ColoredText("{0} color".format(c.capitalize()), c))
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


if __name__ == '__main__':
    test()
