# JCVI: A Versatile Toolkit for Comparative Genomics Analysis

[![Latest PyPI version](https://img.shields.io/pypi/v/jcvi.svg)](https://pypi.python.org/pypi/jcvi)
[![bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/jcvi/README.html?highlight=jcvi)
[![Github Actions](https://github.com/tanghaibao/jcvi/workflows/build/badge.svg)](https://github.com/tanghaibao/jcvi/actions)
[![Downloads](https://pepy.tech/badge/jcvi)](https://pepy.tech/project/jcvi)

Collection of Python libraries to parse bioinformatics files, or perform
computation related to assembly, annotation, and comparative genomics.

|         |                                                                  |
| ------- | ---------------------------------------------------------------- |
| Authors | Haibao Tang ([tanghaibao](http://github.com/tanghaibao))         |
|         | Vivek Krishnakumar ([vivekkrish](https://github.com/vivekkrish)) |
|         | Adam Taranto ([Adamtaranto](https://github.com/Adamtaranto))     |
|         | Xingtan Zhang ([tangerzhang](https://github.com/tangerzhang))    |
|         | Won Cheol Yim ([wyim-pgl](https://github.com/wyim-pgl))          |
| Email   | <tanghaibao@gmail.com>                                           |
| License | [BSD](http://creativecommons.org/licenses/BSD/)                  |

## How to cite

> [!TIP]
> JCVI is now published in iMeta!
>
> _Tang et al. (2024) JCVI: A Versatile Toolkit for Comparative Genomics
> Analysis. [iMeta](https://doi.org/10.1002/imt2.211)_

![MCSCAN example](https://www.dropbox.com/s/9vl3ys3ndvimg4c/grape-peach-cacao.png?raw=1)

![ALLMAPS animation](https://www.dropbox.com/s/jfs8xavcxix37se/ALLMAPS.gif?raw=1)

![GRABSEEDS example](https://www.dropbox.com/s/yu9ehsi6sqifuaa/bluredges.png?raw=1)

## Contents

Following modules are available as generic Bioinformatics handling
methods.

- <kbd>algorithms</kbd>

  - Linear programming solver with SCIP and GLPK.
  - Supermap: find set of non-overlapping anchors in BLAST or NUCMER output.
  - Longest or heaviest increasing subsequence.
  - Matrix operations.

- <kbd>apps</kbd>

  - GenBank entrez accession, Phytozome, Ensembl and SRA downloader.
  - Calculate (non)synonymous substitution rate between gene pairs.
  - Basic phylogenetic tree construction using PHYLIP, PhyML, or RAxML, and viualization.
  - Wrapper for BLAST+, LASTZ, LAST, BWA, BOWTIE2, CLC, CDHIT, CAP3, etc.

- <kbd>formats</kbd>

  Currently supports `.ace` format (phrap, cap3, etc.), `.agp`
  (goldenpath), `.bed` format, `.blast` output, `.btab` format,
  `.coords` format (`nucmer` output), `.fasta` format, `.fastq`
  format, `.fpc` format, `.gff` format, `obo` format (ontology),
  `.psl` format (UCSC blat, GMAP, etc.), `.posmap` format (Celera
  assembler output), `.sam` format (read mapping), `.contig`
  format (TIGR assembly format), etc.

- <kbd>graphics</kbd>

  - BLAST or synteny dot plot.
  - Histogram using R and ASCII art.
  - Paint regions on set of chromosomes.
  - Macro-synteny and micro-synteny plots.
  - Ribbon plots from whole genome alignments.

- <kbd>utils</kbd>
  - Grouper can be used as disjoint set data structure.
  - range contains common range operations, like overlap
    and chaining.
  - Miscellaneous cookbook recipes, iterators decorators,
    table utilities.

Then there are modules that contain domain-specific methods.

- <kbd>assembly</kbd>

  - K-mer histogram analysis.
  - Preparation and validation of tiling path for clone-based assemblies.
  - Scaffolding through ALLMAPS, optical map and genetic map.
  - Pre-assembly and post-assembly QC procedures.

- <kbd>annotation</kbd>

  - Training of _ab initio_ gene predictors.
  - Calculate gene, exon and intron statistics.
  - Wrapper for PASA and EVM.
  - Launch multiple MAKER processes.

- <kbd>compara</kbd>
  - C-score based BLAST filter.
  - Synteny scan (de-novo) and lift over (find nearby anchors).
  - Ancestral genome reconstruction using Sankoff's and PAR method.
  - Ortholog and tandem gene duplicates finder.

## Applications

Please visit [wiki](https://github.com/tanghaibao/jcvi/wiki) for
full-fledged applications.

## Dependencies

JCVI requires Python3 between v3.9 and v3.12.

Some graphics modules require the [ImageMagick](https://imagemagick.org/index.php) library.

On MacOS this can be installed using Conda (see next section). If you are using a linux system (i.e. Ubuntu) you can install ImageMagick using apt-get:

```bash
sudo apt-get update
sudo apt-get install libmagickwand-dev
```

See the [Wand](https://docs.wand-py.org/en/0.2.4/guide/install.html) docs for instructions on installing ImageMagick on other systems.

A few modules may ask for locations of external programs,
if the executable cannot be found in your `PATH`.

The external programs that are often used are:

- [Kent tools](http://hgdownload.cse.ucsc.edu/admin/jksrc.zip)
- [BEDTOOLS](http://code.google.com/p/bedtools/)
- [EMBOSS](http://emboss.sourceforge.net/)

### Managing dependencies with Conda

You can use the the YAML files in this repo to create an environment with basic JCVI dependencies.

If you are new to Conda, we recommend the [Miniforge](https://conda-forge.org/download/) distribution.

```bash
conda env create -f environment.yml

conda activate jcvi
```

Note: If you are using a Mac with an ARM64 (Apple Silicon) processor, some dependencies are not currently available from Bioconda for this architecture.

You can instead create a virtual OSX64 (intel) env like this:

```bash
conda env create -f env_osx64.yml

conda activate jcvi-osx64
```

After activating the Conda environment install JCVI using one of the following options.

## Installation

### Installation options

1) Use pip to install the latest development version directly from this repo.

```bash
pip install git+https://github.com/tanghaibao/jcvi.git
```

2) Install latest release from PyPi.

```bash
pip install jcvi
```

3) Alternatively, if you want to install in development mode.

```bash
git clone git://github.com/tanghaibao/jcvi.git && cd jcvi
pip install -e '.[tests]'
```

### Test Installation

If installed successfully, you can check the version with:

```bash
jcvi --version
```

## Usage

Use `python -m` to call any of the modules installed with JCVI.

Most of the modules in this package contains multiple actions. To use
the `fasta` example:

```console
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
             gc | Plot G+C content distribution
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
       simulate | Simulate random fasta file for testing
           some | Include or exclude a list of records (also performs on .qual file if available)
           sort | Sort the records by IDs, sizes, etc.
        summary | Report the real no of bases and N's in fasta files
           tidy | Normalize gap sizes and remove small components in fasta
      translate | Translate CDS to proteins
           trim | Given a cross_match screened fasta, trim the sequence
      trimsplit | Split sequences at lower-cased letters
           uniq | Remove records that are the same
```

Then you need to use one action, you can just do:

```console
python -m jcvi.formats.fasta extract
```

This will tell you the options and arguments it expects.

**Feel free to check out other scripts in the package, it is not just
for FASTA.**

## Star History

[![Star History
Chart](https://api.star-history.com/svg?repos=tanghaibao/jcvi&type=Date)](https://star-history.com/#tanghaibao/jcvi&Date)
