from __future__ import print_function

import os.path
import json
import webbrowser
from jcvi.apps.base import datadir as _DATADIR

try:
    from matplotlib.colors import LinearSegmentedColormap
except ImportError:  # pragma: no cover
    HAVE_MPL = False
else:
    HAVE_MPL = True


__all__ = (
    "COLOR_MAPS",
    "print_maps",
    "print_all_maps",
    "print_maps_by_type",
    "get_map",
    "get_maps",
    "MAP_TYPES",
    "BrewerMap",
)

_DATAFILE = os.path.join(_DATADIR, "colorbrewer_all_schemes.json")

with open(_DATAFILE, "r") as f:
    COLOR_MAPS = json.load(f)

MAP_TYPES = ("Sequential", "Diverging", "Qualitative")


def print_maps(map_type=None, number=None):
    """
    Print maps by type and/or number of defined colors.

    Parameters
    ----------
    map_type : {'Sequential', 'Diverging', 'Qualitative'}, optional
        Filter output by map type. By default all maps are printed.
    number : int, optional
        Filter output by number of defined colors. By default there is
        no numeric filtering.
    """
    if not map_type and not number:
        print_all_maps()

    elif map_type:
        print_maps_by_type(map_type, number)

    else:
        s = (
            "Invalid parameter combination. "
            "number without map_type is not supported."
        )
        raise ValueError(s)


def print_all_maps():
    """
    Print the name and number of defined colors of all available color maps.
    """
    for t in MAP_TYPES:
        print_maps_by_type(t)


def print_maps_by_type(map_type, number=None):
    """
    Print all available maps of a given type.

    Parameters
    ----------
    map_type : {'Sequential', 'Diverging', 'Qualitative'}
        Select map type to print.
    number : int, optional
        Filter output by number of defined colors. By default there is
        no numeric filtering.
    """
    map_type = map_type.lower().capitalize()
    if map_type not in MAP_TYPES:
        s = "Invalid map type, must be one of {0}".format(MAP_TYPES)
        raise ValueError(s)

    print(map_type)

    map_keys = sorted(COLOR_MAPS[map_type].keys())

    format_str = "{0:8}  :  {1}"

    for mk in map_keys:
        num_keys = sorted(COLOR_MAPS[map_type][mk].keys(), key=int)

        if not number or str(number) in num_keys:
            num_str = "{" + ", ".join(num_keys) + "}"
            print(format_str.format(mk, num_str))


class _ColorMap(object):
    """
    Representation of a color map with matplotlib compatible
    views of the map.

    Parameters
    ----------
    name : str
    map_type : str
    colors : list
        Colors as list of 0-255 RGB triplets.

    Attributes
    ----------
    name : str
    map_type : str
    number : int
        Number of colors in color map.
    colors : list
        Colors as list of 0-255 RGB triplets.
    hex_colors : list
    mpl_colors : list
    mpl_colormap : matplotlib LinearSegmentedColormap
    """

    def __init__(self, name, map_type, colors):
        self.name = name
        self.type = map_type
        self.number = len(colors)
        self.colors = colors

    @property
    def hex_colors(self):
        """
        Colors as a tuple of hex strings. (e.g. '#A912F4')

        """
        hc = []

        for color in self.colors:
            h = "#" + "".join("{0:>02}".format(hex(c)[2:].upper()) for c in color)
            hc.append(h)

        return hc

    @property
    def mpl_colors(self):
        """
        Colors expressed on the range 0-1 as used by matplotlib.

        """
        mc = []

        for color in self.colors:
            mc.append(tuple([x / 255.0 for x in color]))

        return mc

    @property
    def mpl_colormap(self):
        """
        A basic matplotlib color map. If you want to specify keyword arguments
        use the `get_mpl_colormap` method.

        """
        return self.get_mpl_colormap()

    def get_mpl_colormap(self, **kwargs):
        """
        A color map that can be used in matplotlib plots. Requires matplotlib
        to be importable. Keyword arguments are passed to
        `matplotlib.colors.LinearSegmentedColormap.from_list`.

        """
        if not HAVE_MPL:  # pragma: no cover
            raise RuntimeError("matplotlib not available.")

        cmap = LinearSegmentedColormap.from_list(self.name, self.mpl_colors, **kwargs)

        return cmap

    def show_as_blocks(self, block_size=100):
        """
        Show colors in the IPython Notebook using ipythonblocks.

        Parameters
        ----------
        block_size : int, optional
            Size of displayed blocks.

        """
        from ipythonblocks import BlockGrid

        grid = BlockGrid(self.number, 1, block_size=block_size)

        for block, color in zip(grid, self.colors):
            block.rgb = color

        grid.show()


class BrewerMap(_ColorMap):
    """
    Representation of a colorbrewer2 color map with matplotlib compatible
    views of the map.

    Parameters
    ----------
    name : str
    map_type : str
    colors : list
        Colors as list of 0-255 RGB triplets.

    Attributes
    ----------
    name : str
    map_type : str
    number : int
        Number of colors in color map.
    colors : list
        Colors as list of 0-255 RGB triplets.
    colorbrewer2_url : str
    hex_colors : list
    mpl_colors : list
    mpl_colormap : matplotlib LinearSegmentedColormap
    """

    @property
    def colorbrewer2_url(self):
        """
        URL that can be used to view the color map at colorbrewer2.org.

        """
        url = "http://colorbrewer2.org/index.html?type={0}&scheme={1}&n={2}"
        return url.format(self.type.lower(), self.name, self.number)

    def colorbrewer2(self):
        """
        View this color map at colorbrewer2.org. This will open
        colorbrewer2.org in your default web browser.

        """
        webbrowser.open_new_tab(self.colorbrewer2_url)  # pragma: no cover


def get_map(name, map_type, number, reverse=False):
    """
    Return a `BrewerMap` representation of the specified color map.

    Parameters
    ----------
    name : str
        Name of color map. Use `print_maps` to see available color maps.
    map_type : {'Sequential', 'Diverging', 'Qualitative'}
        Select color map type.
    number : int
        Number of defined colors in color map.
    reverse : bool, optional
        Set to True to get the reversed color map.
    """
    number = str(number)
    map_type = map_type.lower().capitalize()

    # check for valid type
    if map_type not in MAP_TYPES:
        s = "Invalid map type, must be one of {0}".format(MAP_TYPES)
        raise ValueError(s)

    # make a dict of lower case map name to map name so this can be
    # insensitive to case.
    # this would be a perfect spot for a dict comprehension but going to
    # wait on that to preserve 2.6 compatibility.
    # map_names = {k.lower(): k for k in COLOR_MAPS[map_type].iterkeys()}
    map_names = dict((k.lower(), k) for k in COLOR_MAPS[map_type].keys())

    # check for valid name
    if name.lower() not in map_names:
        s = "Invalid color map name {0!r} for type {1!r}.\n"
        s = s.format(name, map_type)
        valid_names = [str(k) for k in COLOR_MAPS[map_type].keys()]
        valid_names.sort()
        s += "Valid names are: {0}".format(valid_names)
        raise ValueError(s)

    name = map_names[name.lower()]

    # check for valid number
    if number not in COLOR_MAPS[map_type][name]:
        s = "Invalid number for map type {0!r} and name {1!r}.\n"
        s = s.format(map_type, str(name))
        valid_numbers = [int(k) for k in COLOR_MAPS[map_type][name].keys()]
        valid_numbers.sort()
        s += "Valid numbers are : {0}".format(valid_numbers)
        raise ValueError(s)

    colors = COLOR_MAPS[map_type][name][number]["Colors"]

    if reverse:
        name += "_r"
        colors = [x for x in reversed(colors)]

    return BrewerMap(name, map_type, colors)


def get_maps(map_type, number=None):
    """
    Get all available maps of a given type.

    Parameters
    ----------
    map_type : {'Sequential', 'Diverging', 'Qualitative'}
        Select map type to print.
    number : int, optional
        Filter output by number of defined colors. By default there is
        no numeric filtering.

    Returns
    -------
    maps: name of the maps
    """
    map_type = map_type.lower().capitalize()
    if map_type not in MAP_TYPES:
        s = "Invalid map type, must be one of {0}".format(MAP_TYPES)
        raise ValueError(s)

    map_keys = sorted(COLOR_MAPS[map_type].keys())
    maps = []

    for mk in map_keys:
        num_keys = sorted(COLOR_MAPS[map_type][mk].keys(), key=int)

        if not number or str(number) in num_keys:
            maps.append(get_map(mk, map_type, number))

    return maps


def _load_maps_by_type(map_type):
    """
    Load all maps of a given type into a dictionary.

    Color maps are loaded as BrewerMap objects. Dictionary is
    keyed by map name and then integer numbers of defined
    colors. There is an additional 'max' key that points to the
    color map with the largest number of defined colors.

    Parameters
    ----------
    map_type : {'Sequential', 'Diverging', 'Qualitative'}

    Returns
    -------
    maps : dict of BrewerMap
    """
    seq_maps = COLOR_MAPS[map_type]

    loaded_maps = {}

    for map_name in seq_maps:
        loaded_maps[map_name] = {}

        for num in seq_maps[map_name]:
            inum = int(num)
            colors = seq_maps[map_name][num]["Colors"]

            bmap = BrewerMap(map_name, map_type, colors)

            loaded_maps[map_name][inum] = bmap

        max_num = int(max(seq_maps[map_name].keys(), key=int))
        loaded_maps[map_name]["max"] = loaded_maps[map_name][max_num]

    return loaded_maps


if __name__ == "__main__":
    print_maps()
