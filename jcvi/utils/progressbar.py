#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# progressbar  - Text progress bar library for Python.
# Copyright (c) 2005 Nilton Volpato
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


'''Text progress bar library for Python.

A text progress bar is typically used to display the progress of a long
running operation, providing a visual cue that processing is underway.

The ProgressBar class manages the current progress, and the format of the line
is given by a number of widgets. A widget is an object that may display
differently depending on the state of the progress bar. There are three types
of widgets:
 - a string, which always shows itself

 - a ProgressBarWidget, which may return a different value every time its
   update method is called

 - a ProgressBarWidgetHFill, which is like ProgressBarWidget, except it
   expands to fill the remaining width of the line.

The progressbar module is very easy to use, yet very powerful. It will also
automatically enable features like auto-resizing when the system supports it.
'''

from __future__ import division

import math
import os
import signal
import sys
import time
import datetime

try:
    from fcntl import ioctl
    from array import array
    import termios
except ImportError:
    pass

try:
    from abc import ABCMeta, abstractmethod
except ImportError:
    AbstractWidget = object
    abstractmethod = lambda fn: fn
else:
    AbstractWidget = ABCMeta('AbstractWidget', (object,), {})


__author__ = 'Nilton Volpato'
__author_email__ = 'first-name dot last-name @ gmail.com'
__date__ = '2011-05-14'
__version__ = '2.3'


'''Default ProgressBar widgets'''


def format_updatable(updatable, pbar):
    if hasattr(updatable, 'update'): return updatable.update(pbar)
    else: return updatable


class Widget(AbstractWidget):
    '''The base class for all widgets

    The ProgressBar will call the widget's update value when the widget should
    be updated. The widget's size may change between calls, but the widget may
    display incorrectly if the size changes drastically and repeatedly.

    The boolean TIME_SENSITIVE informs the ProgressBar that it should be
    updated more often because it is time sensitive.
    '''

    TIME_SENSITIVE = False
    __slots__ = ()

    @abstractmethod
    def update(self, pbar):
        '''Updates the widget.

        pbar - a reference to the calling ProgressBar
        '''


class WidgetHFill(Widget):
    '''The base class for all variable width widgets.

    This widget is much like the \\hfill command in TeX, it will expand to
    fill the line. You can use more than one in the same line, and they will
    all have the same width, and together will fill the line.
    '''

    @abstractmethod
    def update(self, pbar, width):
        '''Updates the widget providing the total width the widget must fill.

        pbar - a reference to the calling ProgressBar
        width - The total width the widget must fill
        '''


class Timer(Widget):
    'Widget which displays the elapsed seconds.'

    __slots__ = ('format',)
    TIME_SENSITIVE = True

    def __init__(self, format='Elapsed Time: %s'):
        self.format = format

    @staticmethod
    def format_time(seconds):
        'Formats time as the string "HH:MM:SS".'

        return str(datetime.timedelta(seconds=int(seconds)))


    def update(self, pbar):
        'Updates the widget to show the elapsed time.'

        return self.format % self.format_time(pbar.seconds_elapsed)


class ETA(Timer):
    'Widget which attempts to estimate the time of arrival.'

    TIME_SENSITIVE = True

    def update(self, pbar):
        'Updates the widget to show the ETA or total time when finished.'

        if pbar.currval == 0:
            return 'ETA:  --:--:--'
        elif pbar.finished:
            return 'Time: %s' % self.format_time(pbar.seconds_elapsed)
        else:
            elapsed = pbar.seconds_elapsed
            eta = elapsed * pbar.maxval / pbar.currval - elapsed
            return 'ETA:  %s' % self.format_time(eta)


class FileTransferSpeed(Widget):
    'Widget for showing the transfer speed (useful for file transfers).'

    format = '%6.2f %s%s/s'
    prefixes = ' kMGTPEZY'
    __slots__ = ('unit', 'format')

    def __init__(self, unit='B'):
        self.unit = unit

    def update(self, pbar):
        'Updates the widget with the current SI prefixed speed.'

        if pbar.seconds_elapsed < 2e-6 or pbar.currval < 2e-6: # =~ 0
            scaled = power = 0
        else:
            speed = pbar.currval / pbar.seconds_elapsed
            power = int(math.log(speed, 1000))
            scaled = speed / 1000.**power

        return self.format % (scaled, self.prefixes[power], self.unit)


class AnimatedMarker(Widget):
    '''An animated marker for the progress bar which defaults to appear as if
    it were rotating.
    '''

    __slots__ = ('markers', 'curmark')

    def __init__(self, markers='|/-\\'):
        self.markers = markers
        self.curmark = -1

    def update(self, pbar):
        '''Updates the widget to show the next marker or the first marker when
        finished'''

        if pbar.finished: return self.markers[0]

        self.curmark = (self.curmark + 1) % len(self.markers)
        return self.markers[self.curmark]

# Alias for backwards compatibility
RotatingMarker = AnimatedMarker


class Counter(Widget):
    'Displays the current count'

    __slots__ = ('format',)

    def __init__(self, format='%d'):
        self.format = format

    def update(self, pbar):
        return self.format % pbar.currval


class Percentage(Widget):
    'Displays the current percentage as a number with a percent sign.'

    def update(self, pbar):
        return '%3d%%' % pbar.percentage()


class FormatLabel(Timer):
    'Displays a formatted label'

    mapping = {
        'elapsed': ('seconds_elapsed', Timer.format_time),
        'finished': ('finished', None),
        'last_update': ('last_update_time', None),
        'max': ('maxval', None),
        'seconds': ('seconds_elapsed', None),
        'start': ('start_time', None),
        'value': ('currval', None)
    }

    __slots__ = ('format',)
    def __init__(self, format):
        self.format = format

    def update(self, pbar):
        context = {}
        for name, (key, transform) in self.mapping.items():
            try:
                value = getattr(pbar, key)

                if transform is None:
                   context[name] = value
                else:
                   context[name] = transform(value)
            except: pass

        return self.format % context


class SimpleProgress(Widget):
    'Returns progress as a count of the total (e.g.: "5 of 47")'

    __slots__ = ('sep',)

    def __init__(self, sep=' of '):
        self.sep = sep

    def update(self, pbar):
        return '%d%s%d' % (pbar.currval, self.sep, pbar.maxval)


class Bar(WidgetHFill):
    'A progress bar which stretches to fill the line.'

    __slots__ = ('marker', 'left', 'right', 'fill', 'fill_left')

    def __init__(self, marker='#', left='|', right='|', fill=' ',
                 fill_left=True):
        '''Creates a customizable progress bar.

        marker - string or updatable object to use as a marker
        left - string or updatable object to use as a left border
        right - string or updatable object to use as a right border
        fill - character to use for the empty part of the progress bar
        fill_left - whether to fill from the left or the right
        '''
        self.marker = marker
        self.left = left
        self.right = right
        self.fill = fill
        self.fill_left = fill_left


    def update(self, pbar, width):
        'Updates the progress bar and its subcomponents'

        left, marker, right = (format_updatable(i, pbar) for i in
                               (self.left, self.marker, self.right))

        width -= len(left) + len(right)
        # Marker must *always* have length of 1
        marker *= int(pbar.currval / pbar.maxval * width)

        if self.fill_left:
            return '%s%s%s' % (left, marker.ljust(width, self.fill), right)
        else:
            return '%s%s%s' % (left, marker.rjust(width, self.fill), right)


class ReverseBar(Bar):
    'A bar which has a marker which bounces from side to side.'

    def __init__(self, marker='#', left='|', right='|', fill=' ',
                 fill_left=False):
        '''Creates a customizable progress bar.

        marker - string or updatable object to use as a marker
        left - string or updatable object to use as a left border
        right - string or updatable object to use as a right border
        fill - character to use for the empty part of the progress bar
        fill_left - whether to fill from the left or the right
        '''
        self.marker = marker
        self.left = left
        self.right = right
        self.fill = fill
        self.fill_left = fill_left


class BouncingBar(Bar):
    def update(self, pbar, width):
        'Updates the progress bar and its subcomponents'

        left, marker, right = (format_updatable(i, pbar) for i in
                               (self.left, self.marker, self.right))

        width -= len(left) + len(right)

        if pbar.finished: return '%s%s%s' % (left, width * marker, right)

        position = int(pbar.currval % (width * 2 - 1))
        if position > width: position = width * 2 - position
        lpad = self.fill * (position - 1)
        rpad = self.fill * (width - len(marker) - len(lpad))

        # Swap if we want to bounce the other way
        if not self.fill_left: rpad, lpad = lpad, rpad

        return '%s%s%s%s%s' % (left, lpad, marker, rpad, right)


'''Compatability methods and classes for the progressbar module'''
# Python 3.x (and backports) use a modified iterator syntax
# This will allow 2.x to behave with 3.x iterators
if not hasattr(__builtins__, 'next'):
    def next(iter):
        try:
            # Try new style iterators
            return iter.__next__()
        except AttributeError:
            # Fallback in case of a "native" iterator
            return iter.next()


# Python < 2.5 does not have "any"
if not hasattr(__builtins__, 'any'):
   def any(iterator):
      for item in iterator:
         if item: return True

      return False


class UnknownLength: pass


class ProgressBar(object):
    '''The ProgressBar class which updates and prints the bar.

    A common way of using it is like:
    >>> pbar = ProgressBar().start()
    >>> for i in range(100):
    ...    # do something
    ...    pbar.update(i+1)
    ...
    >>> pbar.finish()

    You can also use a ProgressBar as an iterator:
    >>> progress = ProgressBar()
    >>> for i in progress(some_iterable):
    ...    # do something
    ...

    Since the progress bar is incredibly customizable you can specify
    different widgets of any type in any order. You can even write your own
    widgets! However, since there are already a good number of widgets you
    should probably play around with them before moving on to create your own
    widgets.

    The term_width parameter represents the current terminal width. If the
    parameter is set to an integer then the progress bar will use that,
    otherwise it will attempt to determine the terminal width falling back to
    80 columns if the width cannot be determined.

    When implementing a widget's update method you are passed a reference to
    the current progress bar. As a result, you have access to the
    ProgressBar's methods and attributes. Although there is nothing preventing
    you from changing the ProgressBar you should treat it as read only.

    Useful methods and attributes include (Public API):
     - currval: current progress (0 <= currval <= maxval)
     - maxval: maximum (and final) value
     - finished: True if the bar has finished (reached 100%)
     - start_time: the time when start() method of ProgressBar was called
     - seconds_elapsed: seconds elapsed since start_time and last call to
                        update
     - percentage(): progress in percent [0..100]
    '''

    __slots__ = ('currval', 'fd', 'finished', 'last_update_time',
                 'left_justify', 'maxval', 'next_update', 'num_intervals',
                 'poll', 'seconds_elapsed', 'signal_set', 'start_time',
                 'term_width', 'update_interval', 'widgets', '_time_sensitive',
                 '__iterable')

    _DEFAULT_MAXVAL = 100
    _DEFAULT_TERMSIZE = 80
    _DEFAULT_WIDGETS = [Percentage(), ' ', Bar()]

    def __init__(self, maxval=None, widgets=None, term_width=None, poll=1,
                 left_justify=True, fd=sys.stderr):
        '''Initializes a progress bar with sane defaults'''

        # Don't share a reference with any other progress bars
        if widgets is None:
            widgets = list(self._DEFAULT_WIDGETS)

        self.maxval = maxval
        self.widgets = widgets
        self.fd = fd
        self.left_justify = left_justify

        self.signal_set = False
        if term_width is not None:
            self.term_width = term_width
        else:
            try:
                self._handle_resize()
                signal.signal(signal.SIGWINCH, self._handle_resize)
                self.signal_set = True
            except (SystemExit, KeyboardInterrupt): raise
            except:
                self.term_width = self._env_size()

        self.__iterable = None
        self._update_widgets()
        self.currval = 0
        self.finished = False
        self.last_update_time = None
        self.poll = poll
        self.seconds_elapsed = 0
        self.start_time = None
        self.update_interval = 1


    def __call__(self, iterable):
        'Use a ProgressBar to iterate through an iterable'

        try:
            self.maxval = len(iterable)
        except:
            if self.maxval is None:
                self.maxval = UnknownLength

        self.__iterable = iter(iterable)
        return self


    def __iter__(self):
        return self


    def __next__(self):
        try:
            value = next(self.__iterable)
            if self.start_time is None: self.start()
            else: self.update(self.currval + 1)
            return value
        except StopIteration:
            self.finish()
            raise


    # Create an alias so that Python 2.x won't complain about not being
    # an iterator.
    next = __next__


    def _env_size(self):
        'Tries to find the term_width from the environment.'

        return int(os.environ.get('COLUMNS', self._DEFAULT_TERMSIZE)) - 1


    def _handle_resize(self, signum=None, frame=None):
        'Tries to catch resize signals sent from the terminal.'

        h, w = array('h', ioctl(self.fd, termios.TIOCGWINSZ, '\0' * 8))[:2]
        self.term_width = w


    def percentage(self):
        'Returns the progress as a percentage.'
        return self.currval * 100.0 / self.maxval

    percent = property(percentage)


    def _format_widgets(self):
        result = []
        expanding = []
        width = self.term_width

        for index, widget in enumerate(self.widgets):
            if isinstance(widget, WidgetHFill):
                result.append(widget)
                expanding.insert(0, index)
            else:
                widget = format_updatable(widget, self)
                result.append(widget)
                width -= len(widget)

        count = len(expanding)
        while count:
            portion = max(int(math.ceil(width * 1. / count)), 0)
            index = expanding.pop()
            count -= 1

            widget = result[index].update(self, portion)
            width -= len(widget)
            result[index] = widget

        return result


    def _format_line(self):
        'Joins the widgets and justifies the line'

        widgets = ''.join(self._format_widgets())

        if self.left_justify: return widgets.ljust(self.term_width)
        else: return widgets.rjust(self.term_width)


    def _need_update(self):
        'Returns whether the ProgressBar should redraw the line.'
        if self.currval >= self.next_update or self.finished: return True

        delta = time.time() - self.last_update_time
        return self._time_sensitive and delta > self.poll


    def _update_widgets(self):
        'Checks all widgets for the time sensitive bit'

        self._time_sensitive = any(getattr(w, 'TIME_SENSITIVE', False)
                                    for w in self.widgets)


    def update(self, value=None):
        'Updates the ProgressBar to a new value.'

        if value is not None and value is not UnknownLength:
            if (self.maxval is not UnknownLength
                and not 0 <= value <= self.maxval):

                raise ValueError('Value out of range')

            self.currval = value


        if not self._need_update(): return
        if self.start_time is None:
            raise RuntimeError('You must call "start" before calling "update"')

        now = time.time()
        self.seconds_elapsed = now - self.start_time
        self.next_update = self.currval + self.update_interval
        self.fd.write(self._format_line() + '\r')
        self.last_update_time = now


    def start(self):
        '''Starts measuring time, and prints the bar at 0%.

        It returns self so you can use it like this:
        >>> pbar = ProgressBar().start()
        >>> for i in range(100):
        ...    # do something
        ...    pbar.update(i+1)
        ...
        >>> pbar.finish()
        '''

        if self.maxval is None:
            self.maxval = self._DEFAULT_MAXVAL

        self.num_intervals = max(100, self.term_width)
        self.next_update = 0

        if self.maxval is not UnknownLength:
            if self.maxval < 0: raise ValueError('Value out of range')
            self.update_interval = self.maxval / self.num_intervals


        self.start_time = self.last_update_time = time.time()
        self.update(0)

        return self


    def finish(self):
        'Puts the ProgressBar bar in the finished state.'

        self.finished = True
        self.update(self.maxval)
        self.fd.write('\n')
        if self.signal_set:
            signal.signal(signal.SIGWINCH, signal.SIG_DFL)


examples = []
def example(fn):
    try: name = 'Example %d' % int(fn.__name__[7:])
    except: name = fn.__name__

    def wrapped():
        try:
            sys.stdout.write('Running: %s\n' % name)
            fn()
            sys.stdout.write('\n')
        except KeyboardInterrupt:
            sys.stdout.write('\nSkipping example.\n\n')

    examples.append(wrapped)
    return wrapped


@example
def example0():
    pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=300).start()
    for i in range(300):
        time.sleep(0.01)
        pbar.update(i+1)
    pbar.finish()


@example
def example1():
    widgets = ['Test: ', Percentage(), ' ', Bar(marker=RotatingMarker()),
               ' ', ETA(), ' ', FileTransferSpeed()]
    pbar = ProgressBar(widgets=widgets, maxval=10000000).start()
    for i in range(1000000):
        # do something
        pbar.update(10*i+1)
    pbar.finish()


@example
def example2():
    class CrazyFileTransferSpeed(FileTransferSpeed):
        """It's bigger between 45 and 80 percent."""
        def update(self, pbar):
            if 45 < pbar.percentage() < 80:
                return 'Bigger Now ' + FileTransferSpeed.update(self,pbar)
            else:
                return FileTransferSpeed.update(self,pbar)

    widgets = [CrazyFileTransferSpeed(),' <<<', Bar(), '>>> ',
               Percentage(),' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=10000000)
    # maybe do something
    pbar.start()
    for i in range(2000000):
        # do something
        pbar.update(5*i+1)
    pbar.finish()


@example
def example3():
    widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
    pbar = ProgressBar(widgets=widgets, maxval=10000000).start()
    for i in range(1000000):
        # do something
        pbar.update(10*i+1)
    pbar.finish()


@example
def example4():
    widgets = ['Test: ', Percentage(), ' ',
               Bar(marker='0',left='[',right=']'),
               ' ', ETA(), ' ', FileTransferSpeed()]
    pbar = ProgressBar(widgets=widgets, maxval=500)
    pbar.start()
    for i in range(100,500+1,50):
        time.sleep(0.2)
        pbar.update(i)
    pbar.finish()


@example
def example5():
    pbar = ProgressBar(widgets=[SimpleProgress()], maxval=17).start()
    for i in range(17):
        time.sleep(0.2)
        pbar.update(i + 1)
    pbar.finish()


@example
def example6():
    pbar = ProgressBar().start()
    for i in range(100):
        time.sleep(0.01)
        pbar.update(i + 1)
    pbar.finish()


@example
def example7():
    pbar = ProgressBar()  # Progressbar can guess maxval automatically.
    for i in pbar(range(80)):
        time.sleep(0.01)


@example
def example8():
    pbar = ProgressBar(maxval=80)  # Progressbar can't guess maxval.
    for i in pbar((i for i in range(80))):
        time.sleep(0.01)


@example
def example9():
    pbar = ProgressBar(widgets=['Working: ', AnimatedMarker()])
    for i in pbar((i for i in range(50))):
        time.sleep(.08)


@example
def example10():
    widgets = ['Processed: ', Counter(), ' lines (', Timer(), ')']
    pbar = ProgressBar(widgets=widgets)
    for i in pbar((i for i in range(150))):
        time.sleep(0.1)


@example
def example11():
    widgets = [FormatLabel('Processed: %(value)d lines (in: %(elapsed)s)')]
    pbar = ProgressBar(widgets=widgets)
    for i in pbar((i for i in range(150))):
        time.sleep(0.1)


@example
def example12():
    widgets = ['Balloon: ', AnimatedMarker(markers='.oO@* ')]
    pbar = ProgressBar(widgets=widgets)
    for i in pbar((i for i in range(24))):
        time.sleep(0.3)


@example
def example17():
    widgets = [FormatLabel('Animated Bouncer: value %(value)d - '),
               BouncingBar(marker=RotatingMarker())]

    pbar = ProgressBar(widgets=widgets)
    for i in pbar((i for i in range(180))):
        time.sleep(0.05)


if __name__ == '__main__':
    try:
        for example in examples: example()
    except KeyboardInterrupt:
        sys.stdout('\nQuitting examples.\n')
