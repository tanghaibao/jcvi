JCVI utility libraries
======================

.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.31631.svg
    :target: http://dx.doi.org/10.5281/zenodo.31631
    :alt: DOI

.. image:: https://img.shields.io/pypi/v/jcvi.svg
    :target: https://pypi.python.org/pypi/jcvi
    :alt: Latest PyPI version

.. image:: https://img.shields.io/pypi/dm/jcvi.svg
    :target: https://pypi.python.org/pypi/jcvi
    :alt: Number of PyPI downloads

.. image:: https://scrutinizer-ci.com/g/tanghaibao/jcvi/badges/quality-score.png
    :target: https://scrutinizer-ci.com/g/tanghaibao/jcvi
    :alt: Scrutinizer

Collection of Python libraries to parse bioinformatics files, or perform
computation related to assembly, annotation, and comparative genomics.

:Author: Haibao Tang (`tanghaibao <http://github.com/tanghaibao>`_),
         Vivek Krishnakumar (`vivekkrish <https://github.com/vivekkrish>`_),
         Jingping Li (`Jingping <https://github.com/Jingping>`_)
:Email: tanghaibao@gmail.com
:License: `BSD <http://creativecommons.org/licenses/BSD/>`_

.. contents ::

Contents
---------
Following modules are available as generic Bioinformatics handling methods.

- ``algorithms``
    * Linear programming solver with SCIP and GLPK.
    * Supermap: find set of non-overlapping anchors in BLAST or NUCMER output.
    * Longest or heaviest increasing subsequence.
    * Matrix operations.

- ``apps``
    * GenBank entrez accession, phytozome, ensembl and SRA downloader.
    * Calculate (non)synonymous substitution rate between gene pairs.
    * Basic phylogenetic tree construction using PHYLIP, PhyML, or RAxML, and visualization.
    * Wrapper for BLAST+, LASTZ, LAST, BWA, BOWTIE2, CLC, CDHIT, CAP3, etc.

- ``formats``
    Currently supports ``.ace`` format (phrap, cap3, etc.), ``.agp`` (goldenpath),
    ``.bed`` format, ``.blast`` output, ``.btab`` format, ``.coords`` format (``nucmer`` output),
    ``.fasta`` format, ``.fastq`` format, ``.fpc`` format, ``.gff`` format, ``obo`` format (ontology),
    ``.psl`` format (UCSC blat, GMAP, etc.), ``.posmap`` format (Celera assembler output),
    ``.sam`` format (read mapping), ``.contig`` format (TIGR assembly format), etc.

- ``graphics``
    * BLAST or synteny dot plot.
    * Histogram using R and ASCII art.
    * Paint regions on set of chromosomes.
    * Macro-synteny and micro-synteny plots.

- ``utils``
    * Grouper can be used as disjoint set data structure.
    * range contains common range operations, like overlap and chaining.
    * Sybase connector to JCVI internal database.
    * Miscellaneous cookbook recipes, iterators decorators, table utilities.


Then there are modules that contain domain-specific methods.

- ``assembly``
    * K-mer histogram analysis.
    * Preparation and validation of tiling path for clone-based assemblies.
    * Scaffolding through BAMBUS, optical map and genetic map.
    * Pre-assembly and post-assembly QC procedures.

- ``annotation``
    * Training of *ab initio* gene predictors.
    * Calculate gene, exon and intron statistics.
    * Wrapper for PASA and EVM.
    * Launch multiple MAKER processes.

- ``compara``
    * C-score based BLAST filter.
    * Synteny scan (de-novo) and lift over (find nearby anchors).
    * Ancestral genome reconstruction using Sankoff's and PAR method.
    * Ortholog and tandem gene duplicates finder.


Applications
------------
Please visit `wiki <https://github.com/tanghaibao/jcvi/wiki>`_ for
full-fledged applications. Also visit our `Gallery
<https://github.com/tanghaibao/jcvi/wiki/Gallery>`_ to see our
graphics functionality for the production of publication-ready figures.


Dependencies
-------------
Following are a list of third-party python packages that are used by some
routines in the library. These dependencies are *not* mandatory since they are
only used by a few modules.

* `Biopython <http://www.biopython.org>`_
* `numpy <http://numpy.scipy.org>`_
* `matplotlib <http://matplotlib.org/>`_

There are other Python modules here and there in various scripts. The best way
is to install them via ``pip install`` when you see ``ImportError``.


Installation
------------
The easiest way is to install it via PyPI::

    easy_install jcvi

To install the development version::

    pip install git+git://github.com/tanghaibao/jcvi.git

Alternatively, if you want to install manually::

    cd ~/code  # or any directory of your choice
    git clone git://github.com/tanghaibao/jcvi.git
    export PYTHONPATH=~/code:$PYTHONPATH

Please replace ``~/code`` above with whatever you like, but it must contain ``jcvi``.
To avoid setting ``PYTHONPATH`` everytime, please insert the ``export`` command in your
``.bashrc`` or ``.bash_profile``.

In addition, a few module might ask for locations of external programs, if the extended
cannot be found in your ``PATH``. The external programs that are often used are:

* `Kent tools <http://hgdownload.cse.ucsc.edu/admin/jksrc.zip>`_
* `BEDTOOLS <http://code.google.com/p/bedtools/>`_
* `EMBOSS <http://emboss.sourceforge.net/>`_

Most of the scripts in this package contains multiple actions. To use the
``fasta`` example::

    Usage:
        python -m jcvi.formats.fasta ACTION


    Available ACTIONs:
              clean | Remove irregular chars in FASTA seqs
               diff | Check if two fasta records contain same information
            extract | Given fasta file and seq id, retrieve the sequence in fasta format
              fastq | Combine fasta and qual to create fastq file
             filter | Filter the records by size
             format | Trim accession id to the first space or switch id based on 2-column mapping file
            fromtab | Convert 2-column sequence file to FASTA format
               gaps | Print out a list of gap sizes within sequences
          identical | Given 2 fasta files, find all exactly identical records
                ids | Generate a list of headers
               info | Run `sequence_info` on fasta files
              ispcr | Reformat paired primers into isPcr query format
               join | Concatenate a list of seqs and add gaps in between
         longestorf | Find longest orf for CDS fasta
               pair | Sort paired reads to .pairs, rest to .fragments
        pairinplace | Starting from fragment.fasta, find if adjacent records can form pairs
               pool | Pool a bunch of fastafiles together and add prefix
               qual | Generate dummy .qual file based on FASTA file
             random | Randomly take some records
             sequin | Generate a gapped fasta file for sequin submission
               some | Include or exclude a list of records (also performs on .qual file if available)
               sort | Sort the records by IDs, sizes, etc.
            summary | Report the real no of bases and N's in fasta files
               tidy | Normalize gap sizes and remove small components in fasta
          translate | Translate CDS to proteins
               trim | Given a cross_match screened fasta, trim the sequence
          trimsplit | Split sequences at lower-cased letters
               uniq | Remove records that are the same

Then you need to use one action, you can just do::

    python -m jcvi.formats.fasta extract

This will tell you the options and arguments it expects.

**Feel free to check out other scripts in the package, it is not just for FASTA.**


Reference
---------
Haibao Tang et al. (2015). jcvi: JCVI utility libraries. Zenodo.
`10.5281/zenodo.31631 <http://dx.doi.org/10.5281/zenodo.31631>`_.
