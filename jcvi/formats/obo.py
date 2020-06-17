#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog obo_file

Parses obo_file and plot GO lineage
"""
from __future__ import print_function

import os.path as op
import sys
import logging

from collections import deque
from jcvi.apps.base import download
from jcvi.formats.base import read_until


Sample_GO_URL = "http://current.geneontology.org/ontology/go.obo"


class OBOReader(object):
    """
    parse obo file, usually the most updated can be downloaded from
    http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo

    >>> reader = OBOReader()
    >>> for rec in reader:
            print rec
    """

    def __init__(self, obo_file="go.obo"):
        if not op.exists(obo_file):
            logging.debug(
                "`{}` not found. Downloading from: {}".format(obo_file, Sample_GO_URL)
            )
            download(Sample_GO_URL)

        self.obo_file = obo_file
        self.format_version = None  # e.g., "1.2" of "format-version:" line
        self.data_version = (
            None  # e.g., "releases/2016-07-07" from "data-version:" line
        )
        self.typedefs = {}

    def __iter__(self):
        """Return one GO Term record at a time from an obo file."""
        # Wait to open file until needed. Automatically close file when done.
        with open(self.obo_file) as fstream:
            rec_curr = None  # Stores current GO Term
            typedef_curr = None  # Stores current typedef
            for line in fstream:
                # obo lines start with any of: [Term], [Typedef], /^\S+:/, or /^\s*/
                if self.data_version is None:
                    self._init_obo_version(line)
                if rec_curr is None and line[0:6].lower() == "[term]":
                    rec_curr = GOTerm()
                elif typedef_curr is None and line[0:9].lower() == "[typedef]":
                    typedef_curr = TypeDef()
                elif rec_curr is not None or typedef_curr is not None:
                    line = line.rstrip()  # chomp
                    if line:
                        self._add_to_obj(rec_curr, typedef_curr, line)
                    else:
                        if rec_curr is not None:
                            yield rec_curr
                            rec_curr = None
                        elif typedef_curr is not None:
                            # Save typedef.
                            self.typedefs[typedef_curr.item_id] = typedef_curr
                            typedef_curr = None
            # Return last record, if necessary
            if rec_curr is not None:
                yield rec_curr

    def _add_to_obj(self, rec_curr, typedef_curr, line):
        """Add information on line to GOTerm or Typedef."""
        if rec_curr is not None:
            self._add_to_ref(rec_curr, line)
        else:
            add_to_typedef(typedef_curr, line)

    def _init_obo_version(self, line):
        """Save obo version and release."""
        if line[0:14] == "format-version":
            self.format_version = line[16:-1]
        if line[0:12] == "data-version":
            self.data_version = line[14:-1]

    def _add_to_ref(self, rec_curr, line):
        """Add new fields to the current reference."""
        # Examples of record lines containing ':' include:
        #   id: GO:0000002
        #   name: mitochondrial genome maintenance
        #   namespace: biological_process
        #   def: "The maintenance of ...
        #   is_a: GO:0007005 ! mitochondrion organization
        if line[:4] == "id: ":
            assert not rec_curr.item_id
            item_id = line[4:]
            rec_curr.item_id = item_id
            rec_curr.id = item_id
        elif line[:8] == "alt_id: ":
            rec_curr.alt_ids.add(line[8:])
        elif line[:6] == "name: ":
            assert not rec_curr.name
            rec_curr.name = line[6:]
        elif line[:11] == "namespace: ":
            assert not rec_curr.namespace
            rec_curr.namespace = line[11:]
        elif line[:6] == "is_a: ":
            rec_curr._parents.add(line[6:].split()[0])
        elif line[:13] == "is_obsolete: " and line[13:] == "true":
            rec_curr.is_obsolete = True


# pylint: disable=too-few-public-methods
class TypeDef(object):
    """
        TypeDef term.
    """

    def __init__(self):
        self.id = ""  # GO:NNNNNNN  **DEPRECATED** RESERVED NAME IN PYTHON
        self.item_id = ""  # GO:NNNNNNN (will replace deprecated "id")
        self.name = ""  # description
        self.namespace = ""  # BP, CC, MF
        self._parents = set()  # is_a basestring of parents
        self.parents = set()  # parent records
        self.children = set()  # children records
        self.level = None  # shortest distance from root node
        self.depth = None  # longest distance from root node
        self.is_obsolete = False  # is_obsolete
        self.alt_ids = set()  # alternative identifiers

    def __str__(self):
        ret = []
        ret.append("Typedef - {} ({}):".format(self.item_id, self.name))
        ret.append(
            "  Inverse of: {}".format(self.inverse_of if self.inverse_of else "None")
        )
        if self.transitive_over:
            ret.append("  Transitive over:")
            for txo in self.transitive_over:
                ret.append("    - {}".format(txo))
        return "\n".join(ret)


def add_to_typedef(typedef_curr, obo_line):
    """Add new fields to the current typedef."""
    if obo_line[:4] == "id: ":
        assert not typedef_curr.item_id
        item_id = obo_line[4:]
        typedef_curr.item_id = item_id
    elif obo_line[:6] == "name: ":
        assert not typedef_curr.name
        typedef_curr.name = obo_line[6:]
    elif obo_line[:11] == "namespace: ":
        assert not typedef_curr.namespace
        typedef_curr.namespace = obo_line[11:]
    elif obo_line[17:] == "transitive_over: ":
        field_value = obo_line[17:].split("!")[0].rstrip()
        typedef_curr.transitive_over.append(field_value)
    elif obo_line[12:] == "inverse_of":
        assert not typedef_curr.inverse_of
        field_value = obo_line[12:].split("!")[0].rstrip()
        typedef_curr.inverse_of = field_value


class GOTerm(object):
    """
    GO term, actually contain a lot more properties than interfaced here
    """

    def __init__(self):
        self.id = ""  # GO:NNNNNNN  **DEPRECATED** RESERVED NAME IN PYTHON
        self.item_id = ""  # GO:NNNNNNN (will replace deprecated "id")
        self.name = ""  # description
        self.namespace = ""  # BP, CC, MF
        self._parents = set()  # is_a basestring of parents
        self.parents = set()  # parent records
        self.children = set()  # children records
        self.level = None  # shortest distance from root node
        self.depth = None  # longest distance from root node
        self.is_obsolete = False  # is_obsolete
        self.alt_ids = set()  # alternative identifiers

    header = "\t".join(("#id", "level", "name", "alt_ids"))

    def __str__(self):
        level = "level-{:02d}".format(self.level)
        description = "{} [{}]".format(self.name, self.namespace)
        if self.is_obsolete:
            description += " obsolete"
        alt_ids = ",".join(self.alt_ids)
        return "\t".join((self.id, level, description, alt_ids))

    def __repr__(self):
        return "GOTerm('%s')" % (self.id)

    def has_parent(self, term):
        """Return True if this GO object has a parent GO ID."""
        for parent in self.parents:
            if parent.item_id == term or parent.has_parent(term):
                return True
        return False

    def has_child(self, term):
        """Return True if this GO object has a child GO ID."""
        for parent in self.children:
            if parent.item_id == term or parent.has_child(term):
                return True
        return False

    def get_all_parents(self):
        """Return all parent GO IDs."""
        all_parents = set()
        for parent in self.parents:
            all_parents.add(parent.item_id)
            all_parents |= parent.get_all_parents()
        return all_parents

    def get_all_children(self):
        """Return all children GO IDs."""
        all_children = set()
        for parent in self.children:
            all_children.add(parent.item_id)
            all_children |= parent.get_all_children()
        return all_children

    def get_all_parent_edges(self):
        """Return tuples for all parent GO IDs, containing current GO ID and parent GO ID."""
        all_parent_edges = set()
        for parent in self.parents:
            all_parent_edges.add((self.item_id, parent.item_id))
            all_parent_edges |= parent.get_all_parent_edges()
        return all_parent_edges

    def get_all_child_edges(self):
        """Return tuples for all child GO IDs, containing current GO ID and child GO ID."""
        all_child_edges = set()
        for parent in self.children:
            all_child_edges.add((parent.item_id, self.item_id))
            all_child_edges |= parent.get_all_child_edges()
        return all_child_edges


class GODag(dict):
    def __init__(self, obo_file="go-basic.obo"):
        super(GODag, self).__init__()
        self.load_obo_file(obo_file)
        self.version, self.data_version = self.load_obo_file(obo_file)

    def load_obo_file(self, obo_file):
        """Read obo file. Store results."""
        reader = OBOReader(obo_file)

        # Save alt_ids and their corresponding main GO ID. Add to GODag after populating GO Terms
        alt2rec = {}
        for rec in reader:
            # Save record if:
            #   1) Argument load_obsolete is True OR
            #   2) Argument load_obsolete is False and the GO term is "live" (not obsolete)
            if not rec.is_obsolete:
                self[rec.item_id] = rec
                for alt in rec.alt_ids:
                    alt2rec[alt] = rec

        # Save the typedefs and parsed optional_attrs
        self.typedefs = reader.typedefs

        self._populate_terms()
        self._set_level_depth()

        # Add alt_ids to go2obj
        for goid_alt, rec in alt2rec.items():
            self[goid_alt] = rec
        desc = self._str_desc(reader)
        return desc, reader.data_version

    def _str_desc(self, reader):
        """String containing information about the current GO DAG."""
        data_version = reader.data_version
        if data_version is not None:
            data_version = data_version.replace("releases/", "")
        desc = "{OBO}: fmt({FMT}) rel({REL}) {N:,} GO Terms".format(
            OBO=reader.obo_file,
            FMT=reader.format_version,
            REL=data_version,
            N=len(self),
        )
        return desc

    def _populate_terms(self):
        """Convert GO IDs to GO Term record objects. Populate children."""
        # Make parents and relationships references to the actual GO terms.
        for rec in self.values():
            # Given parent GO IDs, set parent GO Term objects
            rec.parents = set([self[goid] for goid in rec._parents])

            # For each parent GO Term object, add it's child GO Term to the children data member
            for parent_rec in rec.parents:
                parent_rec.children.add(rec)

    def _set_level_depth(self):
        """Set level, depth and add inverted relationships."""

        def _init_level(rec):
            if rec.level is None:
                if rec.parents:
                    rec.level = min(_init_level(rec) for rec in rec.parents) + 1
                else:
                    rec.level = 0
            return rec.level

        def _init_depth(rec):
            if rec.depth is None:
                if rec.parents:
                    rec.depth = max(_init_depth(rec) for rec in rec.parents) + 1
                else:
                    rec.depth = 0
            return rec.depth

        def _init_reldepth(rec):
            if not hasattr(rec, "reldepth"):
                up_terms = rec.get_goterms_upper()
                if up_terms:
                    rec.reldepth = max(_init_reldepth(rec) for rec in up_terms) + 1
                else:
                    rec.reldepth = 0
            return rec.reldepth

        for rec in self.values():

            if rec.level is None:
                _init_level(rec)

            if rec.depth is None:
                _init_depth(rec)

    def populate_terms(self):
        def depth(rec):
            if rec.level < 0:
                if not rec.parents:
                    rec.level = 0
                else:
                    rec.level = min(depth(rec) for rec in rec.parents) + 1
            return rec.level

        # make the parents references to the GO terms
        for rec in self.values():
            rec.parents = [self[x] for x in rec._parents]

        # populate children and levels
        for rec in self.values():
            for p in rec.parents:
                p.children.append(rec)

            if rec.level < 0:
                depth(rec)

    def write_dag(self, out=sys.stdout):
        """Write info for all GO Terms in obo file, sorted numerically."""
        print(GOTerm.header, file=out)
        for _, rec in sorted(self.items()):
            print(rec, file=out)

    def query_term(self, term, verbose=False):
        """Given a GO ID, return GO object."""
        if term not in self:
            sys.stderr.write("Term %s not found!\n" % term)
            return

        rec = self[term]
        if verbose:
            print(rec)
            sys.stderr.write("all parents: {}\n".format(repr(rec.get_all_parents())))
            sys.stderr.write("all children: {}\n".format(repr(rec.get_all_children())))
        return rec

    def _label_wrap(self, label):
        wrapped_label = r"%s\n%s" % (label, self[label].name.replace(",", r"\n"))
        return wrapped_label

    def draw_lineage(
        self,
        recs,
        nodecolor="mediumseagreen",
        edgecolor="lightslateblue",
        dpi=96,
        verbose=False,
        lineage_img="GO_lineage.png",
    ):
        # draw AMIGO style network, lineage containing one query record
        import pygraphviz as pgv

        G = pgv.AGraph()
        edgeset = set()
        for rec in recs:
            edgeset.update(rec.get_all_parent_edges())
            edgeset.update(rec.get_all_child_edges())

        edgeset = [(self._label_wrap(a), self._label_wrap(b)) for (a, b) in edgeset]

        for src, target in edgeset:
            """
            Default layout in graphviz is top->bottom,
            so we invert the direction and plot using `dir="back"`
            """
            G.add_edge(target, src)

        G.graph_attr.update(dpi="%d" % dpi)
        G.node_attr.update(
            shape="box", style="rounded,filled", fillcolor="beige", color=nodecolor
        )
        G.edge_attr.update(shape="normal", color=edgecolor, dir="back", label="is_a")
        # highlight the query terms
        for rec in recs:
            try:
                q = G.get_node(self._label_wrap(rec.id))
                q.attr.update(fillcolor="plum")
            except:
                continue

        if verbose:
            print(G.to_string(), file=sys.stderr)

        print(
            "lineage info for terms %s written to %s"
            % ([rec.id for rec in recs], lineage_img),
            file=sys.stderr,
        )

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


if __name__ == "__main__":

    import optparse

    p = optparse.OptionParser(__doc__)
    p.add_option(
        "--term",
        dest="term",
        help="write the parents and children of the query term",
        default=None,
    )

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
