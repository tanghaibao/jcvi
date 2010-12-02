JCVI utility libraries
======================
My own collection of Python libraries to parse files, or perform
assembly-related calculations. Documentations will be lagging behind.

``formats/``
    File parsers for various files used in genome assembly and comparisons. 
    Currents supports ``.agp`` (goldenpath), ``.coords`` (``nucmer`` output), 
    ``.bed`` format, ``.fasta`` format, and tabular blast output. 

``utils/``
    Data structures to simplify programming tasks. ``Grouper`` object can be
    used as disjoint set, ``range`` contains common range operations.

``apps/``
    Helper library to wrap command line programs and run jobs on JCVI grid
    engine.
