
JCVI utility libraries
======================
My own collection of Python libraries to parse files, or perform
assembly-related calculations. Documentations will be lagging behind.

:Author: Haibao Tang (`tanghaibao <http://github.com/tanghaibao>`_)
:Email: tanghaibao@gmail.com
:License: `BSD <http://creativecommons.org/licenses/BSD/>`_

.. contents ::

Contents
---------
- ``algorithms``
    Algorithms for math intensive stuff, including:

    * Linear programming solver with SCIP and GLPK.
    * Synteny scan (de-novo) and lift over (find nearby anchors).
    * Supermap: find set of non-overlapping anchors in BLAST or NUCMER output.
    * Tandem gene duplicates finder.

- ``assembly``
    Scripts to prepare input data to assembler, and also post-assembly
    scaffolding, quality control, etc. In general, anything related to genome
    assembly and scaffolding:

    * K-mer histogram analysis
    * Prepare frg for Celera Assembler (CA)
    * Helper scripts for fixing unitig layout errors in CA
    * QC of potential scaffolding errors

- ``apps``
    Helper library to wrap command line programs and run jobs on JCVI grid
    engine (split jobs, check status, etc.). Driver scripts including:

    * BLAST filter that selects subset of anchors.
    * GenBank entrez accession downloader.
    * LASTZ wrapper.
    * Low complexity sequence masker with NCBI WindowMasker.

- ``formats``
    File parsers for various files used in genome assembly and comparisons. 
    Currents supports ``.agp`` (goldenpath), ``.bed`` format, 
    ``.blast`` output, ``.btab`` format, ``.cas`` (CLC assembler output),
    ``.coords`` format (``nucmer`` output), ``.fasta`` format, ``.fastq`` format, 
    ``.fpc`` format, ``.gff`` format, ``obo`` format (ontology),
    ``.posmap`` format (Celera assembler output), ``.sam`` format (read
    mapping).

- ``graphics``
    Graphics to visualize comparative genomics or assembly stuff. Including:

    * Genome assembly A50 plot.
    * BLAST or synteny dot plot.
    * Histogram using R.
    * Painting regions on set of chromosomes.

- ``utils``
    Data structures to simplify programming tasks. Most of the scripts are
    derived from ideas in the public domain, and are commonly used by other
    modules.  For example:

    * Grouper can be used as disjoint set data structure.
    * range contains common range operations, like overlap and chaining.
    * Sybase connector.
    * Miscellaneous cookbook recipes.


Dependencies
-------------
Following are a list of third-party python packages that are used by some
routines in the library. These dependencies are *not* mandatory since they are
only used by a few modules.

* `Biopython <http://www.biopython.org>`_
* `numpy <http://numpy.scipy.org>`_
* `scipy <http://www.scipy.org>`_


Installation
------------
Resolve dependencies first. Then place the whole folder ``jcvi/`` on your
``PYTHONPATH``. Most scripts can both ``import`` or run as utility script. 
For example::

    python -m jcvi.formats.fasta

Or just copy ``jcvi`` to the current folder (since Python searches current
folder by default)::

    python jcvi/formats/fasta.py   

Most of the scripts in this package contains multiple actions. To use the
``fasta`` example::

    available actions:
        `extract`: given fasta file and seq id, retrieve the sequence in fasta format
        `summary`: report the real no of bases and N's in fastafiles
        `uniq`: remove records that are the same
        `format`: trim accession id to the first space
        `random`: random take some records
        `diff`: check if two FASTA records contain same information
        `trim`: given a cross_match screened fasta, trim the sequence
        `pair`: sort paired reads to .pairs.fasta and remaining to .fragments.fasta
        `fastq`: combine fasta and qual to create fastq file
        `sequin`: generated a gapped fasta file suitable for sequin submission
        `gaps`: print out a list of gap sizes within sequences
        `join`: concatenate a list of seqs and add gaps in between
        `some`: include or exclude a list of records (also performs on .qual file if available)

Then you need to use one action, you can just do::

    python -m jcvi.formats.fasta extract

This will tell you the options and arguments it expects. 

**Feel free to check out other scripts in the package, it is not just for FASTA.**
 
