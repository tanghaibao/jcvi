# cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

"""
Cythonized (fast) version of BlastLine

Stolen from brentp's biostuff (thanks):
<https://github.com/brentp/bpbio/blob/master/biostuff/biostuff/cblastline.pyx>
"""
import sys
from libc.stdio cimport FILE, EOF, fopen, fscanf, rewind, fclose, sscanf, \
            fgets, sprintf
from libc.string cimport strcpy


cdef const char *blast_format = "%s\t%s\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%f"
cdef const char *blast_format_line = "%s\t%s\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%f\n"
cdef const char *blast_output = "%s\t%s\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2g\t%.3g"
cdef const char *bed_output = "%s\t%d\t%d\t%s:%d-%d\t%.2g\t%c"


cdef class Blast:
    cdef:
        FILE* fh
        object filename

    def __cinit__(self, char* filename):
        self.fh = fopen(filename, 'r')
        self.filename = filename

    def __iter__(self):
        rewind(self.fh)
        return self

    def __next__(self):
        cdef:
            float pct = 0.0, evalue = 0.0, bit = 0.0
            char qname[128]
            char sname[128]
            int hlen, nmiss, ngap, qstart, qstop, sstart, sstop
            char *tmp
            int success

        success = fscanf(self.fh, blast_format_line, qname, sname, \
                         &pct, &hlen, &nmiss, &ngap, &qstart, &qstop,\
                         &sstart, &sstop, &evalue, &bit )
        if success == EOF:
            raise StopIteration
        return create_blast_line(qname, sname, pct, hlen, nmiss, ngap,
                        qstart, qstop, sstart, sstop, evalue, bit)

    def __dealloc__(self):
        fclose(self.fh)

    def __repr__(self):
        return "Blast('%s')" % (self.filename, )

# Python 2 and 3 differ in str and unicode handling
# https://github.com/PySlurm/pyslurm/wiki/Strings-and-bytes-in-Cython
cdef bytes c_str(str s):
    return s.encode("UTF-8")

cdef str py_str(bytes s):
    return s.decode("UTF-8", "replace")


cdef class BlastLine:
    """
    Given a string of tab-delimited (-m 8) blast output, parse it and create
    an object with the usual attrs:

    >>> b = BlastLine("Os09g11510	Os08g13650	92.31	39	3	0	2273	2311	3237	3199	0.001	54.0")
    >>> b.query
    'Os09g11510'
    >>> attrs = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
    ...  'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score')
    >>> [getattr(b, attr) for attr in attrs]  # doctest: +ELLIPSIS
    ['Os09g11510', 'Os08g13650', 92.3..., 39, 3, 0, 2273, 2311, 3237, 3199, 0.001..., 54.0]
    """

    __slots__ = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
                 'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score', \
                 'qseqid', 'sseqid', 'qi', 'si', 'orientation')

    cdef public:
        char _query[128]
        char _subject[128]
        int hitlen, nmismatch, ngaps, qstart, qstop, sstart, sstop
        float pctid, score
        double evalue
        object qseqid, sseqid
        int qi, si
        char orientation

    property query:
        def __get__(self):
            return py_str(self._query)
        def __set__(self, val: str):
            strcpy(self._query, c_str(val))

    property subject:
        def __get__(self):
            return py_str(self._subject)
        def __set__(self, val: str):
            strcpy(self._subject, c_str(val))

    def __init__(self, s):
        sline = c_str(s)
        sscanf(sline, blast_format, self._query, self._subject,
                &self.pctid, &self.hitlen, &self.nmismatch, &self.ngaps,
                &self.qstart, &self.qstop,
                &self.sstart, &self.sstop,
                &self.evalue, &self.score)

        self.orientation = '+'
        if self.qstart > self.qstop:
            self.qstart, self.qstop = self.qstop, self.qstart
            self.orientation = '-'
        if self.sstart > self.sstop:
            self.sstart, self.sstop = self.sstop, self.sstart
            self.orientation = '-'

    def __richcmp__(BlastLine self, BlastLine other, size_t op):
        if op == 2: # ==
            if self.query != other.query and self.qstart != other.qstart:
                return False
            return self.subject == other.subject and \
                    self.qstop == other.qstop and \
                    self.sstop == other.sstop and \
                    self.evalue == other.evalue and \
                    self.hitlen == other.hitlen

        elif op == 3: # !=
            return not self.__richcmp__(other, 2)
        else:
            raise Exception("that comparison not implemented")

    def __hash__(self):
        return id(self)

    def __repr__(self):
        return "BlastLine('%s' to '%s', eval=%.3f, score=%.1f)" % \
                (self.query, self.subject, self.evalue, self.score)

    def __str__(self):
        args = [getattr(self, attr) for attr in BlastLine.__slots__[:12]]
        if self.orientation == '-':
            args[8], args[9] = args[9], args[8]

        cdef char result[512]
        sprintf(result, blast_output, self._query, self._subject,
            self.pctid, self.hitlen, self.nmismatch, self.ngaps,
            self.qstart, self.qstop,
            self.sstart, self.sstop,
            self.evalue, self.score)

        return py_str(result)

    @property
    def has_score(self):
        return hasattr(self, "score")

    @property
    def swapped(self):
        """
        Swap query and subject.
        """
        args = [getattr(self, attr) for attr in BlastLine.__slots__[:12]]
        args[0:2] = [self.subject, self.query]
        args[6:10] = [self.sstart, self.sstop, self.qstart, self.qstop]
        if self.orientation == '-':
            args[8], args[9] = args[9], args[8]
        b = "\t".join(str(x) for x in args)
        return BlastLine(b)

    @property
    def bedline(self):
        cdef char result[512]
        sprintf(result, bed_output,
                self._subject, self.sstart - 1, self.sstop,
                self._query, self.qstart, self.qstop,
                self.score, self.orientation)

        return py_str(result)

    def __reduce__(self):
        return create_blast_line, (
            self.query, self.subject, self.pctid, self.hitlen, self.nmismatch,
            self.ngaps, self.qstart, self.qstop, self.sstart, self.sstop,
            self.evalue, self.score)


cdef BlastLine create_blast_line(char *query, char *subject, float pctid, int hitlen,
                       int nmismatch, int ngaps, int qstart, int qstop,
                       int sstart, int sstop, float evalue, float score):
    """ Factory method.
    """
    cdef BlastLine b = BlastLine.__new__(BlastLine)
    b.query = query
    b.subject = subject
    b.pctid = pctid
    b.hitlen = hitlen
    b.nmismatch = nmismatch
    b.ngaps = ngaps
    b.qstart = qstart
    b.qstop = qstop
    b.sstart = sstart
    b.sstop = sstop
    b.evalue = evalue
    b.score = score
    return b
