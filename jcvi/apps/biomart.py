#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Builds the queries for Globus and BioMart servie, usefu for extraction of
phytozome data sets.  Certain portion of the codes are ported from R package
`biomaRt` (thanks).
"""
from __future__ import print_function

import sys
import urllib
import logging

from urllib.parse import urljoin

from xml.etree.ElementTree import ElementTree, Element, SubElement, tostring

from jcvi.apps.base import OptionParser, ActionDispatcher, download


class GlobusXMLParser(ElementTree):
    def __init__(self, xml_file):
        """Parse an Globus directory listing XML file

        Args:
            xml_file (str): Path to the XML file
        """
        with open(xml_file) as fp:
            self.parse(fp)

    def get_genomes(self):
        """
        Only folders containing `assembly` and `annotation` are of interest.
        """
        root = PhytozomePath(next(self.iter(tag="organismDownloads")))
        genomes = {}
        for child in root.values():
            if child.has_genome_release:
                genomes[child.name] = child

        # early_release
        early_release = root.get("early_release")
        if early_release:
            for child in early_release.values():
                if child.has_genome_release:
                    genomes[child.name] = child

        return genomes


class PhytozomePath(dict):
    TAGS_OF_INTEREST = ("organismDownloads", "folder", "file")

    def __init__(self, element):
        """Deserialize XML => dict-like structure to ease navigation
        between folders. Keys are folder or file names.

        Args:
            element (ElementTree): XML parse tree
        """
        tag = element.tag
        assert tag in self.TAGS_OF_INTEREST
        self.url = None
        if tag == "file":
            self.name = element.attrib["filename"]
            self.url = element.attrib["url"]
        else:
            self.name = element.attrib["name"]
        self.tag = tag
        for child in element.getchildren():
            if child.tag not in self.TAGS_OF_INTEREST:
                continue
            child = PhytozomePath(child)
            self[child.name] = child

    @property
    def has_genome_release(self):
        """Only the folders that contain both `assembly` and `annotation` are of interest here.
        """
        return "assembly" in self and "annotation" in self

    def download(self, name, base_url, cookies, downloader=None):
        """Download the file if it has an URL. Otherwise, this will recursively search the children.

        See also:
        <https://genome.jgi.doe.gov/portal/help/download.jsf>

        Args:
            name (str, optional): Name of the file. Defaults to None.
            base_url (str): Link to the file on the internet.
            cookies (str, optional): cookies file. Defaults to None.
            downloader (str, optional): Use a given downloader. One of wget|curl|powershell|insecure.
            Defaults to None.
        """
        if self.name == name and base_url and self.url:
            url = urljoin(base_url, self.url)
            download(
                url, filename=name, debug=True, cookies=cookies, downloader=downloader
            )
        else:
            for child_name, child in self.items():
                if child_name == name:
                    child.download(name, base_url, cookies, downloader=downloader)
        return name

    def __repr__(self):
        return "{}: [{}]".format(self.name, ", ".join(repr(v) for v in self))


class MartXMLParser(ElementTree):
    def __init__(self, xml_data):
        self.parse(xml_data)

    def parse_marts(self):
        for t in self.getiterator("MartURLLocation"):
            if t.attrib["visible"] == "1":
                yield Mart(**t.attrib)

    def parse_configuration(self):
        # the attributes
        for t in self.getiterator("AttributeDescription"):
            yield Attribute(**t.attrib)

        # the filters
        for t in self.getiterator("FilterDescription"):
            f = Filter(**t.attrib)
            options = [Option(**x.attrib) for x in t.getiterator("Option")]
            f.add_options(options)
            yield f


class Mart(dict):
    def __init__(
        self,
        host="www.biomart.org",
        path="/biomart/martservice",
        port="80",
        name="ensembl",
        virtual_schema="default",
        **attrib
    ):

        self.__dict__ = attrib.copy()
        self.__dict__.update(
            x for x in locals().items() if x[0] not in ("self", "attrib")
        )

        self.registry = {}
        self.url = "http://{0}:{1}{2}".format(self.host, self.port, path)
        self.display_name = self.__dict__.get("displayName", "")
        self.virtual_schema = self.__dict__.get(
            "serverVirtualSchema", self.virtual_schema
        )

    def __str__(self):
        return "\t".join((self.name, self.display_name, self.virtual_schema))

    def get_registry(self, archive=False):
        type = "registry_archive" if archive else "registry"
        params = urllib.urlencode(dict(type=type))
        xml_data = urllib.urlopen(self.url, params)

        parser = MartXMLParser(xml_data)
        for t in parser.parse_marts():
            self.registry[t.name] = t

    def list_registry(self):
        if len(self.registry) == 0:
            self.get_registry()
        for m in sorted(self.registry.values()):
            print(m)

    def get_datasets(self):
        params = urllib.urlencode(dict(type="datasets", mart=self.name))
        web_data = urllib.urlopen(self.url, params)

        for row in web_data:
            atoms = row.strip().split("\t")
            if atoms[0] == "TableSet":
                name, description, last_updated = atoms[1], atoms[2], atoms[-1]
                self[name] = Dataset(name, description, last_updated, self)

    def list_datasets(self):
        if len(self) == 0:
            self.get_datasets()
        for m in sorted(self.values(), key=str):
            print(m)


class Dataset(object):
    """
    Connect to a specified dataset in the database
    """

    def __init__(self, name, description, last_updated, mart):
        self.name = name
        self.description = description
        self.last_updated = last_updated
        self.mart = mart

        self.attributes = {}
        self.filters = {}

    def __str__(self):
        return "\t".join((self.name, self.description, self.last_updated))

    def get_configuration(self):
        params = urllib.urlencode(dict(type="configuration", dataset=self.name))
        xml_data = urllib.urlopen(self.mart.url, params)

        parser = MartXMLParser(xml_data)
        for t in parser.parse_configuration():
            if isinstance(t, Attribute):
                self.attributes[t.internalName] = t
            elif isinstance(t, Filter):
                self.filters[t.internalName] = t

    def list_attributes(self):
        if len(self.attributes) == 0:
            self.get_configuration()
        for m in sorted(self.attributes.values()):
            print(m)

    def list_filters(self):
        if len(self.filters) == 0:
            self.get_configuration()
        for m in sorted(self.filters.values()):
            print(m)

    def query(self, filters={}, attributes=()):
        q = MartQuery(dataset=self)
        q.add_filters(**filters)
        q.add_attributes(attributes)
        return q.execute()


class MartQuery(object):
    def __init__(
        self, dataset=None, formatter="TSV", header="0", unique_rows="0", count="0"
    ):
        self.dataset = dataset
        self.url = dataset.mart.url
        self.virtual_schema = dataset.mart.virtual_schema
        self.formatter = formatter
        self.header = header
        self.unique_rows = unique_rows
        self.count = count
        self.name = dataset.name
        self.attributes = []
        self.filters = {}

    def add_filters(self, **filters):
        for key, val in filters.items():
            self.filters[key] = str(val)

    def add_attributes(self, attributes):
        for key in attributes:
            self.attributes.append(key)

    def set_header(self, flag):
        self.header = str(flag)

    def set_formatter(self, format="TSV"):
        self.formatter = format

    def build_query(self):
        query_t = Element(
            "Query",
            dict(
                virtualSchemaName=self.virtual_schema,
                formatter=self.formatter,
                header=self.header,
                uniqueRows=self.unique_rows,
                count=self.count,
                datasetConfigVersion="0.6",
            ),
        )
        dataset_t = SubElement(
            query_t, "Dataset", dict(name=self.name, interface="default")
        )
        for key, val in self.filters.items():
            SubElement(dataset_t, "Filter", dict(name=key, value=val))
        for attribute in self.attributes:
            SubElement(dataset_t, "Attribute", dict(name=attribute))

        return tostring(query_t)

    def execute(self, debug=False):
        xml_data = self.build_query()
        if debug:
            print(xml_data, file=sys.stderr)
        data = urllib.urlencode(dict(query=xml_data))
        return urllib.urlopen(self.url, data)


class MartArgument(object):
    def __init__(self, **attrib):
        self.__dict__ = attrib.copy()

    def __str__(self):
        return self.__class__.__name__ + str(self.__dict__)


class Attribute(MartArgument):
    """
    Attributes define the values that we are retrieving.

    For example, the gene start, stop, or chromosomes it belongs to
    """

    pass


class Filter(MartArgument):
    """
    Filters define a restriction on the query.

    For example, you can restrict output to all genes located on chr. 1
    then use the filter chromosome_name with value `1`
    """

    def add_options(self, options):
        self.options = dict((x.displayName, x) for x in options)


class Option(MartArgument):
    pass


class Sequence(object):
    def __init__(self, seq):
        self.seq = seq

    def export_fasta(self):
        pass


def test_biomart():
    bm = Mart()
    bm.list_registry()
    bm.list_datasets()
    return bm


def get_ensembl_dataset():
    bm = Mart()
    ensembl = bm.registry["ensembl"]
    ensembl.get_datasets()
    dataset = ensembl["mmusculus_gene_ensembl"]
    return dataset


def get_phytozome_dataset():
    # Either of the following method is okay
    # bm = Mart()
    # phytozome = bm.registry["phytozome_mart"]

    phytozome = Mart(
        host="www.phytozome.net",
        port="80",
        name="phytozome_mart",
        virtual_schema="zome_mart",
    )

    phytozome.get_datasets()
    dataset = phytozome["phytozome"]
    return dataset


def get_bed_from_phytozome(genelist):
    """
    >>> data = get_bed_from_phytozome(["AT5G54690", "AT1G01010"])
    >>> print data.read()  #doctest: +NORMALIZE_WHITESPACE
    Chr1	3631	5899	AT1G01010
    Chr5	22219224	22221840	AT5G54690
    <BLANKLINE>
    """
    genelist = ",".join(genelist)
    dataset = get_phytozome_dataset()
    filters = dict(gene_name_filter=genelist)
    attributes = "chr_name1,gene_chrom_start,gene_chrom_end,gene_name1".split(",")

    data = dataset.query(filters=filters, attributes=attributes)
    return data


def main():

    actions = (("bed", "get gene bed from phytozome"),)
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bed(args):
    """
    %prog bed genes.ids

    Get gene bed from phytozome. `genes.ids` contains the list of gene you want
    to pull from Phytozome. Write output to .bed file.
    """
    p = OptionParser(bed.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (idsfile,) = args
    ids = set(x.strip() for x in open(idsfile))
    data = get_bed_from_phytozome(list(ids))

    pf = idsfile.rsplit(".", 1)[0]
    bedfile = pf + ".bed"
    fw = open(bedfile, "w")
    for i, row in enumerate(data):
        row = row.strip()
        if row == "":
            continue

        print(row, file=fw)

    logging.debug("A total of {0} records written to `{1}`.".format(i + 1, bedfile))


if __name__ == "__main__":

    import doctest

    doctest.testmod()

    main()
