#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Disjoint set data structure <http://code.activestate.com/recipes/387776/>
Author: Michael Droettboom
"""


class Grouper(object):
    """
    This class provides a lightweight way to group arbitrary objects
    together into disjoint sets when a full-blown graph data structure
    would be overkill.

    Objects can be joined using .join(), tested for connectedness
    using .joined(), and all disjoint sets can be retrieved using list(g)
    The objects being joined must be hashable.

    >>> g = Grouper()
    >>> g.join('a', 'b')
    >>> g.join('b', 'c')
    >>> g.join('d', 'e')
    >>> list(g)
    [['a', 'b', 'c'], ['d', 'e']]
    >>> g.joined('a', 'b')
    True
    >>> g.joined('a', 'c')
    True
    >>> 'f' in g
    False
    >>> g.joined('a', 'd')
    False
    >>> del g['b']
    >>> list(g)
    [['a', 'c'], ['d', 'e']]
    """

    def __init__(self, init=[]):
        mapping = self._mapping = {}
        for x in init:
            mapping[x] = [x]

    def join(self, a, *args):
        """
        Join given arguments into the same set. Accepts one or more arguments.
        """
        mapping = self._mapping
        set_a = mapping.setdefault(a, [a])

        for arg in args:
            set_b = mapping.get(arg)
            if set_b is None:
                set_a.append(arg)
                mapping[arg] = set_a
            elif set_b is not set_a:
                if len(set_b) > len(set_a):
                    set_a, set_b = set_b, set_a
                set_a.extend(set_b)
                for elem in set_b:
                    mapping[elem] = set_a

    def joined(self, a, b):
        """
        Returns True if a and b are members of the same set.
        """
        mapping = self._mapping
        try:
            return mapping[a] is mapping[b]
        except KeyError:
            return False

    def __iter__(self):
        """
        Returns an iterator returning each of the disjoint sets as a list.
        """
        seen = set()
        for elem, group in self._mapping.items():
            if elem not in seen:
                yield group
                seen.update(group)

    def __getitem__(self, key):
        """
        Returns the set that a certain key belongs.
        """
        return tuple(self._mapping[key])

    def __contains__(self, key):
        return key in self._mapping

    def __len__(self):
        group = set()
        for v in self._mapping.values():
            group.update([tuple(v)])
        return len(group)

    def __delitem__(self, key):
        group = self._mapping[key]
        group.remove(key)
        del self._mapping[key]

    @property
    def num_members(self):
        return sum(len(x) for x in self)

    def keys(self):
        return self._mapping.keys()


if __name__ == "__main__":
    import doctest

    doctest.testmod()
