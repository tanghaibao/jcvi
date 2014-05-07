
JCVI utility libraries
======================
My own collection of Python libraries to parse files, or perform
assembly-related calculations. Documentations will be lagging behind.

:Author: Haibao Tang (`tanghaibao <http://github.com/tanghaibao>`_),
         Vivek Krishnakumar (`vivekkrish <https://github.com/vivekkrish>`_),
         Jingping Li (`Jingping <https://github.com/Jingping>`_),
         Maria Kim (`msarmien <https://github.com/msarmien>`_),
         Xingtan Zhang (`tangerzhang <https://github.com/tangerzhang>`_)
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
    * GenBank entrez accession and Phytozome downloader.
    * Calculate (non)synonymous substitution rate between gene pairs.
    * Basic phylogenetic tree construction using PHYLIP, PhyML, or RAxML, and visualization.
    * Wrapper for BLAST+, LASTZ, LAST, BWA, BOWTIE2, CLC, CDHIT, CAP3, etc.

- ``formats``
    Currently supports ``.ace`` format (phrap, cap3, etc.), ``.agp`` (goldenpath),
    ``.bed`` format, ``.blast`` output, ``.btab`` format, ``.cas`` (CLC assembler output),
    ``.coords`` format (``nucmer`` output), ``.fasta`` format, ``.fastq`` format,
    ``.fpc`` format, ``.gff`` format, ``obo`` format (ontology),
    ``.psl`` format (UCSC blat, GMAP, etc.), ``.posmap`` format (Celera assembler output),
    ``.sam`` format (read mapping), ``.contig`` format (TIGR assembly format), etc.

- ``graphics``
    * BLAST or synteny dot plot.
    * Histogram using R and ASCII art.
    * Painting regions on set of chromosomes.
    * Heatmap from csv file.

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

- ``variation``
    * Convert between various flavors of SNP datasets.
    * Read deconvolution into taxa or samples.
    * Launch TASSEL pipeline.


Dependencies
-------------
Following are a list of third-party python packages that are used by some
routines in the library. These dependencies are *not* mandatory since they are
only used by a few modules.

* `Biopython <http://www.biopython.org>`_
* `numpy <http://numpy.scipy.org>`_
* `scipy <http://www.scipy.org>`_

There are other Python modules here and there in various scripts. The best way
is to install them via ``pip install`` when you see ``ImportError``.


Installation
------------
Resolve dependencies first. Then place the whole folder ``jcvi/`` on your
``PYTHONPATH``. Most scripts can both ``import`` or run as utility script. *This
is the preferred method*, as you can run regardless of the dir you are in::

    export PYTHONPATH=dir_contains_jcvi:$PYTHONPATH
    python -m jcvi.formats.fasta

Please replace ``dir_contains_jcvi`` above with whatever you like, but it must
contain ``jcvi``. To avoid setting ``PYTHONPATH`` everytime, please insert the last
command in your ``.bashrc`` or ``.bash_profile``.

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
        `extract`: given fasta file and seq id, retrieve the sequence in fasta format
        `longestorf`: find longest orf for CDS fasta
        `translate`: translate CDS to proteins
        `info`: run `sequence_info` on fasta files
        `summary`: report the real no of bases and N's in fasta files
        `uniq`: remove records that are the same
        `ids`: generate a list of headers
        `format`: trim accession id to the first space or switch id based on 2-column mapping file
        `pool`: pool a bunch of fastafiles together and add prefix
        `random`: randomly take some records
        `diff`: check if two fasta records contain same information
        `identical`: given 2 fasta files, find all exactly identical records
        `trim`: given a cross_match screened fasta, trim the sequence
        `sort`: sort the records by IDs, sizes, etc.
        `filter`: filter the records by size
        `pair`: sort paired reads to .pairs, rest to .fragments
        `pairinplace`: starting from fragment.fasta, find if adjacent records can form pairs
        `fastq`: combine fasta and qual to create fastq file
        `tidy`: normalize gap sizes and remove small components in fasta
        `sequin`: generate a gapped fasta file for sequin submission
        `gaps`: print out a list of gap sizes within sequences
        `join`: concatenate a list of seqs and add gaps in between
        `some`: include or exclude a list of records (also performs on .qual file if available)
        `clean`: remove irregular chars in FASTA seqs
        `ispcr`: reformat paired primers into isPcr query format
        `fromtab`: convert 2-column sequence file to FASTA format

Then you need to use one action, you can just do::

    python -m jcvi.formats.fasta extract

This will tell you the options and arguments it expects.

**Feel free to check out other scripts in the package, it is not just for FASTA.**
