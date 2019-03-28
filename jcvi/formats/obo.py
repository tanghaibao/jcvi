#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog obo_file

Parses obo_file and plot GO lineage
"""
from __future__ import print_function

import sys
import logging

from collections import deque
from jcvi.formats.base import read_until


typedef_tag, term_tag = "[Typedef]", "[Term]"


def after_colon(line):
    # macro for getting anything after the :
    return line.split(":", 1)[1].strip()


class OBOReader (object):
    """
    parse obo file, usually the most updated can be downloaded from
    http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo

    >>> reader = OBOReader()
    >>> for rec in reader:
            print rec
    """

    def __init__(self, obo_file="gene_ontology.1_2.obo"):

        try:
            self._handle = open(obo_file)
        except IOError as e:
            logging.error(e)
            logging.error("download obo first: " +
                          "http://geneontology.org/ontology/obo_format_1_2/")
            sys.exit(1)

    def __iter__(self):

        line = self._handle.readline()
        if not line.startswith(term_tag):
            read_until(self._handle, term_tag)
        while 1:
            yield next(self)

    def next(self):

        lines = []
        line = self._handle.readline()
        if not line or line.startswith(typedef_tag):
            raise StopIteration

        # read until the next tag and save everything in between
        while 1:
            pos = self._handle.tell()  # save current postion for roll-back
            line = self._handle.readline()
            if line.startswith(typedef_tag) or line.startswith(term_tag):
                self._handle.seek(pos)  # roll-back
                break
            lines.append(line)

        rec = GOTerm()
        for line in lines:
            if line.startswith("id:"):
                rec.id = after_colon(line)
            if line.startswith("alt_id:"):
                rec.alt_ids.append(after_colon(line))
            elif line.startswith("name:"):
                rec.name = after_colon(line)
            elif line.startswith("namespace:"):
                rec.namespace = after_colon(line)
            elif line.startswith("is_a:"):
                rec._parents.append(after_colon(line).split()[0])
            elif line.startswith("is_obsolete:") \
                    and after_colon(line) == "true":
                rec.is_obsolete = True

        return rec


class GOTerm (object):
    """
    GO term, actually contain a lot more properties than interfaced here
    """

    def __init__(self):
        self.id = ""              # GO:xxxxxx
        self.name = ""            # description
        self.namespace = ""       # BP, CC, MF
        self._parents = []        # is_a string of parents
        self.parents = []         # parent records
        self.children = []        # children records
        self.level = -1           # distance from root node
        self.is_obsolete = False  # is_obsolete
        self.alt_ids = []         # alternative identifiers

    def __str__(self):
        obsolete = "obsolete" if self.is_obsolete else ""
        return "%s\tlevel-%02d\t%s [%s] %s" % \
            (self.id, self.level, self.name,
             self.namespace, obsolete)

    def __repr__(self):
        return "GOTerm('%s')" % (self.id)

    def has_parent(self, term):
        for p in self.parents:
            if p.id == term or p.has_parent(term):
                return True
        return False

    def has_child(self, term):
        for p in self.children:
            if p.id == term or p.has_child(term):
                return True
        return False

    def get_all_parents(self):
        all_parents = set()
        for p in self.parents:
            all_parents.add(p.id)
            all_parents |= p.get_all_parents()
        return all_parents

    def get_all_children(self):
        all_children = set()
        for p in self.children:
            all_children.add(p.id)
            all_children |= p.get_all_children()
        return all_children

    def get_all_parent_edges(self):
        all_parent_edges = set()
        for p in self.parents:
            all_parent_edges.add((self.id, p.id))
            all_parent_edges |= p.get_all_parent_edges()
        return all_parent_edges

    def get_all_child_edges(self):
        all_child_edges = set()
        for p in self.children:
            all_child_edges.add((p.id, self.id))
            all_child_edges |= p.get_all_child_edges()
        return all_child_edges


class GODag(dict):

    def __init__(self, obo_file="gene_ontology.1_2.obo"):

        self.load_obo_file(obo_file)
        self.valid_names = set(x.name for x in self.values())

    def load_obo_file(self, obo_file):

        logging.debug("load obo file `%s`" % obo_file)
        obo_reader = OBOReader(obo_file)
        for rec in obo_reader:
            self[rec.id] = rec
            for alt in rec.alt_ids:
                self[alt] = rec

        self.populate_terms()
        logging.debug("%d nodes imported" % len(self))

    def populate_terms(self):

        def depth(rec):
            if rec.level < 0:
                if not rec.parents:
                    rec.level = 0
                else:
                    rec.level = min(depth(rec) for rec in rec.parents) + 1
            return rec.level

        # make the parents references to the GO terms
        for rec in self.itervalues():
            rec.parents = [self[x] for x in rec._parents]

        # populate children and levels
        for rec in self.itervalues():
            for p in rec.parents:
                p.children.append(rec)

            if rec.level < 0:
                depth(rec)

    def write_dag(self, out=sys.stdout):

        for rec_id, rec in sorted(self.items()):
            print(rec, file=out)

    def query_term(self, term, verbose=False):
        try:
            rec = self[term]
        except:
            print("Term %s not found!" % term, file=sys.stderr)
            return
        print(rec, file=sys.stderr)
        if verbose:
            print("all parents:", rec.get_all_parents(), file=sys.stderr)
            print("all children:", rec.get_all_children(), file=sys.stderr)

        return rec

    def _label_wrap(self, label):
        wrapped_label = r"%s\n%s" % (label,
                                     self[label].name.replace(",", r"\n"))
        return wrapped_label

    def draw_lineage(self, recs, nodecolor="mediumseagreen",
                     edgecolor="lightslateblue", dpi=96, verbose=False,
                     lineage_img="GO_lineage.png"):
        # draw AMIGO style network, lineage containing one query record
        import pygraphviz as pgv

        G = pgv.AGraph()
        edgeset = set()
        for rec in recs:
            edgeset.update(rec.get_all_parent_edges())
            edgeset.update(rec.get_all_child_edges())

        edgeset = [(self._label_wrap(a), self._label_wrap(b))
                   for (a, b) in edgeset]

        for src, target in edgeset:
            """
            Default layout in graphviz is top->bottom,
            so we invert the direction and plot using `dir="back"`
            """
            G.add_edge(target, src)

        G.graph_attr.update(dpi="%d" % dpi)
        G.node_attr.update(shape="box", style="rounded,filled",
                           fillcolor="beige", color=nodecolor)
        G.edge_attr.update(shape="normal", color=edgecolor,
                           dir="back", label="is_a")
        # highlight the query terms
        for rec in recs:
            try:
                q = G.get_node(self._label_wrap(rec.id))
                q.attr.update(fillcolor="plum")
            except:
                continue

        if verbose:
            print(G.to_string(), file=sys.stderr)

        print("lineage info for terms %s written to %s" %
              ([rec.id for rec in recs], lineage_img), file=sys.stderr)

        G.draw(lineage_img, prog="dot")

    def update_association(self, association):
        bad_terms = set()
        for key, terms in association.items():
            parents = set()
            for term in terms:
                try:
                    parents.update(self[term].get_all_parents())
                except:
                    bad_terms.add(term)
            terms.update(parents)
        if bad_terms:
            print("terms not found:", bad_terms, file=sys.stderr)


def load_GODag():
    """
    OBO file retrieved from http://obo.cvs.sourceforge.net/viewvc/obo/obo/ontology/genomic-proteomic/so.obo
    """
    from jcvi.apps.base import download

    so_file_url = "http://obo.cvs.sourceforge.net/viewvc/obo/obo/ontology/genomic-proteomic/so.obo"
    so_file = download(so_file_url, debug=False)

    return GODag(so_file)


def validate_term(term, so=None, method="verify"):
    """
    Validate an SO term against so.obo
    """
    if so is None:
        so = load_GODag()

    oterm = term
    if term not in so.valid_names:
        if "resolve" in method:
            if "_" in term:
                tparts = deque(term.split("_"))
                tparts.pop() if "prefix" in method else tparts.popleft()
                nterm = "_".join(tparts).strip()
                term = validate_term(nterm, so=so, method=method)
            if term is None:
                return None
        else:
            logging.error("Term `{0}` does not exist".format(term))
            sys.exit(1)

    if oterm != term:
        logging.debug("Resolved term `{0}` to `{1}`".format(oterm, term))
    return term


if __name__ == '__main__':

    import optparse
    p = optparse.OptionParser(__doc__)
    p.add_option("--term", dest="term", help="write the parents and children"
                 "of the query term", default=None)

    opts, args = p.parse_args()

    if len(args) != 1:
        sys.exit(p.print_help())

    obo_file = args[0]
    g = GODag(obo_file)
    g.write_dag()

    # run a test case
    if opts.term:
        rec = g.query_term(opts.term, verbose=True)
        g.draw_lineage([rec], verbose=True)
