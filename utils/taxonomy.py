"""
From my blog post:
<http://tanghaibao.blogspot.com/2010/02/getting-phylogeny-from-list-of.html>

Example:
>>> mylist = [3702, 3649, 3694, 3880]
>>> t = TaxIDTree(mylist)
>>> print t
(((Carica_papaya,Arabidopsis_thaliana)Brassicales,(Medicago_truncatula,Populus_trichocarpa)fabids)rosids);
>>> t.print_tree()
<BLANKLINE>
               /-Carica_papaya

          /---|

         |     \-Arabidopsis_thaliana

---- /---|

         |     /-Medicago_truncatula

          \---|

               \-Populus_trichocarpa
"""
import time
import logging

from urllib2 import urlopen
from ClientForm import ParseResponse
from BeautifulSoup import BeautifulSoup

from jcvi.utils.cbook import memoized


URL = "http://itol.embl.de/other_trees.shtml"


class TaxIDTree(object):

    def __init__(self, list_of_taxids):
        # If only one taxid provided, get full tree with nameExp
        # else, get default tree
        if isinstance(list_of_taxids, int):  # single taxon
            list_of_taxids = [list_of_taxids]
            form_element_id = "nameExp"
        else:
            form_element_id = "nameCol"

        # the data to send in
        form_data = "\n".join(str(x) for x in list_of_taxids)

        success = False
        while not success:
            try:
                response = urlopen(URL)
                success = True
            except (urllib2.URLError, urllib2.HTTPError,
                    RuntimeError, KeyError) as e:
                logging.error(e)
                logging.debug("wait 5 seconds to reconnect...")
                time.sleep(5)

        forms = ParseResponse(response, backwards_compat=False)
        form = forms[0]

        form["ncbiIDs"] = form_data
        page = urlopen(form.click()).read()
        soup = BeautifulSoup(page)

        for element in soup("textarea"):

            if element["id"] == form_element_id:
                self.newick = str(element.contents[0])

    def __str__(self):
        return self.newick

    def print_tree(self):
        from ete2 import Tree
        t = Tree(self.newick, format=8)
        print t


def get_names(list_of_taxids):
    """
    >>> mylist = [3702, 3649, 3694, 3880]
    >>> get_names(mylist)
    ['Arabidopsis thaliana', 'Carica papaya', 'Populus trichocarpa', 'Medicago truncatula']
    """

    list_of_taxids = [str(x) for x in list_of_taxids]
    from jcvi.apps.entrez import batch_taxonomy

    return list(batch_taxonomy(list_of_taxids))


@memoized
def MRCA(list_of_taxids):
    """
    This gets the most recent common ancester (MRCA) for a list of taxids

    >>> mylist = [3702, 3649, 3694, 3880]
    >>> MRCA(mylist)
    'rosids'
    """

    from ete2 import Tree

    t = TaxIDTree(list_of_taxids)
    t = Tree(str(t), format=8)

    leaves = [x.name for x in t.get_leaves()]
    ancestor = t.get_common_ancestor(*t.get_leaves())

    return ancestor.name


@memoized
def isPlantOrigin(taxid):
    """
    Given a taxid, this gets the expanded tree which can then be checked to
    see if the organism is a plant or not

    >>> isPlantOrigin(29760)
    True
    """

    assert isinstance(taxid, int)

    t = TaxIDTree(taxid)
    try:
        return "Viridiplantae" in str(t)
    except AttributeError:
        raise ValueError("{0} is not a valid ID".format(taxid))


if __name__ == '__main__':

    import doctest
    #doctest.testmod()
    print "Testing isPlantOrigin():"
    print 3702, isPlantOrigin(3702)    # Arabidopsis thaliana
    print 10090, isPlantOrigin(10090)  # Mus musculus

    print "\nTest cache by 10K calls:"
    for i in xrange(10000):
        isPlantOrigin(3702)
        isPlantOrigin(10090)
    print "done"

    print "\nTest invalid ID:"
    print 10099, isPlantOrigin(10099)  # Wrong ID
