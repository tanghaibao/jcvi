"""
Cythonized (fast) version of BlastLine

Stolen from brentp's biostuff (thanks):
<https://github.com/brentp/bpbio/blob/master/biostuff/biostuff/cblastline.pyx>
"""
cdef extern from *:
    ctypedef char* const_char_star "const char*"

import sys
cimport libc.stdlib
from libc.stdio cimport FILE, EOF, fopen, fscanf, rewind, fclose, sscanf, \
            fgets, sprintf


cdef extern from "Python.h":
    char *PyString_AsString(object)


cdef const_char_star blast_format = "%s\t%s\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%f"
cdef const_char_star blast_format_line = "%s\t%s\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%f\n"
cdef const_char_star blast_output = "%s\t%s\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2g\t%.3g"


cdef class Blast:
    cdef FILE* fh
    cdef object filename
    def __cinit__(self, char* filename):
        self.fh = fopen(filename, 'r')
        self.filename = filename

    def __iter__(self):
        rewind(self.fh)
        return self

    def __next__(self):
        cdef float pct = 0.0, evalue = 0.0, bit = 0.0
        cdef char qname[128]
        cdef char sname[128]
        cdef int hlen, nmiss, ngap, qstart, qstop, sstart, sstop
        cdef char *tmp
        cdef int success
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
        return "BlastFile('%s')" % (self.filename, )


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
    cdef public int hitlen, nmismatch, ngaps, qstart, qstop, sstart, sstop
    cdef public float pctid, score
    cdef public double evalue
    cdef public char orientation
    cdef public int qi, si

    cdef char _cqseqid[32]
    cdef char _csseqid[32]
    cdef object _pyqseqid, _pysseqid
    cdef char _cquery[128]
    cdef char _csubject[128]
    cdef object _pysubject, _pyquery
    __slots__ = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
                 'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score', \
                 'qseqid', 'sseqid', 'qi', 'si', 'orientation')

    property qseqid:
        def __get__(self):
            if self._pyqseqid is None:
                return self._cqseqid
            return self._pyqseqid
        def __set__(self, val):
            self._pyqseqid = val

    property sseqid:
        def __get__(self):
            if self._pysseqid is None:
                return self._csseqid
            return self._pysseqid
        def __set__(self, val):
            self._pysseqid = val

    property query:
        def __get__(self):
            if self._pyquery is None:
                return self._cquery
            return self._pyquery
        def __set__(self, val):
            self._pyquery = val

    property subject:
        def __get__(self):
            if self._pysubject is None:
                return self._csubject
            return self._pysubject
        def __set__(self, val):
            self._pysubject = val

    def __init__(self, char *sline):
        if sline != NULL:
            sscanf(sline, blast_format, self._cquery, self._csubject,
                &self.pctid, &self.hitlen, &self.nmismatch, &self.ngaps,
                &self.qstart, &self.qstop,
                &self.sstart, &self.sstop,
                &self.evalue, &self.score)

        if self.sstart > self.sstop:
            self.sstart, self.sstop = self.sstop, self.sstart
            self.orientation = '-'
        else:
            self.orientation = '+'

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

    def __repr__(self):
        return "BlastLine('%s' to '%s', eval=%.3f, score=%.1f)" % \
                (self.query, self.subject, self.evalue, self.score)

    def __str__(self):
        args = [getattr(self, attr) for attr in BlastLine.__slots__[:12]]
        if self.orientation == '-':
            args[8], args[9] = args[9], args[8]

        cdef char result[512]
        sprintf(result, blast_output, self._cquery, self._csubject,
            self.pctid, self.hitlen, self.nmismatch, self.ngaps,
            self.qstart, self.qstop,
            self.sstart, self.sstop,
            self.evalue, self.score)

        return str(result)

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
        return "\t".join(str(x) for x in \
                (self.subject, self.sstart - 1, self.sstop, self.query,
                 self.score, self.orientation))

    def __reduce__(self):
        return create_blast_line, (
            self.query, self.subject, self.pctid, self.hitlen, self.nmismatch,
            self.ngaps, self.qstart, self.qstop, self.sstart, self.sstop,
            self.evalue, self.score)


cdef extern from "pnew.h":
        cdef BlastLine NEW_BLASTLINE "PY_NEW" (object t)


cpdef BlastLine create_blast_line(char *query, char*subject, float pctid, int hitlen,
                       int nmismatch, int ngaps, int qstart, int qstop,
                       int sstart, int sstop, float evalue, float score):
    """ Factory method.
    """
    cdef BlastLine b = NEW_BLASTLINE(BlastLine)
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
