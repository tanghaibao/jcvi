#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog blast_file --qbed query.bed --sbed subject.bed

accepts .bed format: <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>
and a blast file.

local dup filter:
if the input is query.bed and subject.bed, the script files query.localdups and subject.localdups are created containing the parent|offspring dups, as inferred by subjects hitting the same query or queries hitting the same subject.

repeat filter:
adjust the evalues in a dagchainer/blast file by the number of times they occur.
query/subjects that appear often will have the evalues raise (made less significant).
adjusted_evalue(A, B) = evalue(A, B) ** ((counts_of_blast / counts_of_genes) / (counts(A) + counts(B)))

cscore filter:
see supplementary info for sea anemone genome paper <http://www.sciencemag.org/cgi/content/abstract/317/5834/86>, formula below
cscore(A,B) = score(A,B) / max(best score for A, best score for B)

Finally a .raw file (which is the input for the quota-align pipeline <http://github.com/tanghaibao/quota-alignment/>) is created
"""

import sys
import logging
import os.path as op
import collections
import itertools

from math import log10

from jcvi.formats.bed import Bed
from jcvi.formats.blast import BlastLine 
from jcvi.utils.grouper import Grouper
from jcvi.utils.cbook import gene_name
from jcvi.apps.base import debug
debug()


def main(blast_file, options):

    qbed_file, sbed_file = options.qbed, options.sbed

    # is this a self-self blast?
    is_self = (qbed_file == sbed_file)
    if is_self:
        logging.debug("... looks like a self-self BLAST to me")
    
    global_density_ratio = options.global_density_ratio
    tandem_Nmax = options.tandem_Nmax
    filter_repeats = options.filter_repeats
    cscore = options.cscore

    logging.debug("read bed files %s and %s" % (qbed_file, sbed_file))
    qbed = Bed(qbed_file)
    sbed = Bed(sbed_file)

    qorder = qbed.order
    sorder = sbed.order

    fp = file(blast_file)
    total_lines = sum(1 for line in fp)
    logging.debug("read BLAST file %s (total %d lines)" % \
            (blast_file, total_lines))
    fp.seek(0)
    blasts = sorted([BlastLine(line) for line in fp], \
            key=lambda b: b.score, reverse=True)

    filtered_blasts = []
    seen = set() 
    ostrip = options.strip_names
    for b in blasts:
        query, subject = b.query, b.subject
        if ostrip:
            query, subject = gene_name(query), gene_name(subject)
        if query not in qorder:
            logging.debug("WARNING: %s not in %s" % (query, qbed.filename))
            continue
        if subject not in sorder:
            logging.debug("WARNING: %s not in %s" % (subject, sbed.filename))
            continue

        qi, q = qorder[query]
        si, s = sorder[subject]
        
        if is_self and qi > si:
            # move all hits to same side when doing self-self BLAST
            query, subject = subject, query
            qi, si = si, qi
            q, s = s, q

        key = query, subject
        if key in seen: continue
        seen.add(key)
        b.query, b.subject = key

        b.qi, b.si = qi, si
        b.qseqid, b.sseqid = q['seqid'], s['seqid']
        
        filtered_blasts.append(b)


    if global_density_ratio:
        logging.debug("running the global_density filter" + \
                "(global_density_ratio=%d)..." % options.global_density_ratio)
        gene_count = len(qorder) + len(sorder)
        before_filter = len(filtered_blasts)
        filtered_blasts = filter_to_global_density(filtered_blasts, gene_count,
                                                   global_density_ratio)
        logging.debug("after filter (%d->%d)..." % (before_filter,
            len(filtered_blasts)))

    if tandem_Nmax:
        logging.debug("running the local dups filter (tandem_Nmax=%d)..." % \
                tandem_Nmax)

        qtandems = tandem_grouper(qbed, filtered_blasts,
                flip=True, tandem_Nmax=tandem_Nmax)
        standems = tandem_grouper(sbed, filtered_blasts, 
                flip=False, tandem_Nmax=tandem_Nmax)

        qdups_fh = open(op.splitext(qbed_file)[0] + ".localdups", "w")

        if is_self:
            for s in standems: qtandems.join(*s)
            qdups_to_mother = write_localdups(qdups_fh, qtandems, qbed)
            sdups_to_mother = qdups_to_mother
        else:
            qdups_to_mother = write_localdups(qdups_fh, qtandems, qbed)
            sdups_fh = open(op.splitext(sbed_file)[0] + ".localdups", "w")
            sdups_to_mother = write_localdups(sdups_fh, standems, sbed)

        if options.tandems_only:
            # just want to use this script as a tandem finder.
            sys.exit()

        # write out new .bed after tandem removal
        write_new_bed(qbed, qdups_to_mother)
        if not is_self:
            write_new_bed(sbed, sdups_to_mother)
        
        before_filter = len(filtered_blasts)
        filtered_blasts = list(filter_tandem(filtered_blasts, \
                qdups_to_mother, sdups_to_mother))
        logging.debug("after filter (%d->%d)..." % \
                (before_filter, len(filtered_blasts)))

        qnew_name = "%s.nolocaldups%s" % op.splitext(qbed.filename)
        snew_name = "%s.nolocaldups%s" % op.splitext(sbed.filename)

        qbed_new = Bed(qnew_name)
        sbed_new = Bed(snew_name)

        qorder = qbed_new.order
        sorder = sbed_new.order

    if filter_repeats:
        before_filter = len(filtered_blasts)
        logging.debug("running the repeat filter")
        filtered_blasts = list(filter_repeat(filtered_blasts))
        logging.debug("after filter (%d->%d)..." % (before_filter,
            len(filtered_blasts)))

    if cscore:
        before_filter = len(filtered_blasts)
        logging.debug("running the cscore filter (cscore>=%.2f)..." % cscore)
        filtered_blasts = list(filter_cscore(filtered_blasts, cscore=cscore))
        logging.debug("after filter (%d->%d)..." % (before_filter,
            len(filtered_blasts)))

    # this is the final output we will write to after BLAST filters
    raw_name = "%s.raw" % op.splitext(blast_file)[0]
    raw_fh = open(raw_name, "w")

    write_raw(qorder, sorder, filtered_blasts, raw_fh)

    if options.write_filtered_blast:
        write_new_blast(filtered_blasts)


def write_localdups(dups_fh, tandems, bed):

    logging.debug("write local dups to file %s" % dups_fh.name)

    tandem_groups = []
    for group in tandems:
        rows = [bed[i] for i in group]
        # within the tandem groups, genes are sorted with decreasing size
        rows.sort(key=lambda a: (-abs(a['end'] - a['start']), a['accn']))
        tandem_groups.append([row['accn'] for row in rows])

    dups_to_mother = {}
    for accns in sorted(tandem_groups):
        print >>dups_fh, "\t".join(accns)
        for dup in accns[1:]:
            dups_to_mother[dup] = accns[0]

    return dups_to_mother


def write_new_bed(bed, children):
    # generate local dup removed annotation files
    out_name = "%s.nolocaldups%s" % op.splitext(bed.filename)
    logging.debug("write tandem-filtered bed file %s" % out_name)
    fh = open(out_name, "w")
    for i, row in enumerate(bed):
        if row['accn'] in children: continue
        print >>fh, row
    fh.close()


def write_raw(qorder, sorder, filtered_blasts, raw_fh):

    logging.debug("write raw file %s" % raw_fh.name)
    for b in filtered_blasts: 
        qi, q = qorder[b.query]
        si, s = sorder[b.subject]
        qseqid, sseqid = q['seqid'], s['seqid']

        score = 50 if b.evalue == 0 else min(int(-log10(b.evalue)), 50)
        print >>raw_fh, "\t".join(map(str, (qseqid, qi, sseqid, si, score)))


def write_new_blast(filtered_blasts, fh=sys.stdout):
    for b in filtered_blasts:
        print >>fh, b

"""
All BLAST filters
"""

def filter_to_global_density(blast_list, gene_count, global_density_ratio):
    max_hits = int(gene_count * global_density_ratio)
    loigging.debug("cutting at: %d" % max_hits)
    return blast_list[:max_hits]

def filter_cscore(blast_list, cscore=.5):

    best_score = {}
    for b in blast_list:
        if b.query not in best_score or b.score > best_score[b.query]:
            best_score[b.query] = b.score
        if b.subject not in best_score or b.score > best_score[b.subject]:
            best_score[b.subject] = b.score

    for b in blast_list:
        cur_cscore = b.score / max(best_score[b.query], best_score[b.subject])
        if cur_cscore > cscore:
            yield b


def filter_repeat(blast_list, evalue_cutoff=.05):
    """
    adjust the evalues in a dagchainer/blast file by the number of times they occur.
    query/subjects that appear often will have the evalues raise (made less
    significant).
    """
    counts = collections.defaultdict(int)
    for b in blast_list:
        counts[b.query] += 1
        counts[b.subject] += 1

    expected_count = len(blast_list) * 1. / len(counts)
    logging.debug("(expected_count=%d)..." % expected_count)

    for b in blast_list:
        count = counts[b.query] + counts[b.subject]
        adjusted_evalue = b.evalue ** (expected_count / count)

        if adjusted_evalue < evalue_cutoff: yield b


def filter_tandem(blast_list, qdups_to_mother, sdups_to_mother):
    
    mother_blast = []
    for b in blast_list:
        if b.query in qdups_to_mother: b.query = qdups_to_mother[b.query]
        if b.subject in sdups_to_mother: b.subject = sdups_to_mother[b.subject]
        mother_blast.append(b)
    
    mother_blast.sort(key=lambda b: b.score, reverse=True)
    seen = {}
    for b in mother_blast:
        if b.query==b.subject: continue
        key = b.query, b.subject
        if key in seen: continue
        seen[key] = None
        yield b


def tandem_grouper(bed, blast_list, tandem_Nmax=10, flip=True):
    if not flip:
        simple_blast = [(b.query, (b.sseqid, b.si)) for b in blast_list if b.evalue < 1e-10] 
    else:
        simple_blast = [(b.subject, (b.qseqid, b.qi)) for b in blast_list if b.evalue < 1e-10] 

    simple_blast.sort()

    standems = Grouper()
    for name, hits in itertools.groupby(simple_blast, key=lambda x:x[0]):
        # these are already sorted.
        hits = [x[1] for x in hits]
        for ia, a in enumerate(hits[:-1]):
            b = hits[ia + 1]
            # on the same chromosome and rank difference no larger than tandem_Nmax
            if b[1] - a[1] <= tandem_Nmax and b[0] == a[0]: 
                standems.join(a[1], b[1])

    return standems


if __name__ == "__main__":
    import optparse

    parser = optparse.OptionParser(__doc__)
    parser.add_option("--qbed", dest="qbed", 
            help="path to qbed")
    parser.add_option("--sbed", dest="sbed", 
            help="path to sbed")
    parser.add_option("--no_strip_names", dest="strip_names", action="store_false", default=True,
            help="do not strip alternative splicing (e.g. At5g06540.1 -> At5g06540)")
    parser.add_option("--tandems-only", dest="tandems_only", action="store_true", default=False,
            help="only calculate tandems, write .localdup file and exit.")

    filter_group = optparse.OptionGroup(parser, "BLAST filters")
    filter_group.add_option("--tandem_Nmax", dest="tandem_Nmax", type="int", default=None, 
            help="merge tandem genes within distance [default: %default]")
    filter_group.add_option("--filter_repeats", dest="filter_repeats", action="store_true", default=False,
            help="require higher e-value for repetitive matches BLAST.")
    filter_group.add_option("--cscore", type="float", default=None,
            help="retain hits that have good bitscore. a value of 0.5 means "
                 "keep all values that are 50% or greater of the best hit. "
                 "higher is more stringent [default: %default]")
    filter_group.add_option("--global_density_ratio", type="float", default=None,
            help="maximum ratio of blast hits to genes a good value is 10. "
                 "if there are more blasts, only the those with the lowest "
                 "are kept. [default: %default]")
    filter_group.add_option("--write-filtered-blast", dest="write_filtered_blast",
                            action="store_true", default=False,
                            help="write a new blast with localdups and reps filtered")

    parser.add_option_group(filter_group)

    (options, blast_files) = parser.parse_args()

    if not (len(blast_files) == 1 and options.qbed and options.sbed):
        sys.exit(parser.print_help())

    main(blast_files[0], options)

