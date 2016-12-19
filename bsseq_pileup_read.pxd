from pysam.libcalignedsegment cimport AlignedSegment, PileupRead, makeAlignedSegment
from pysam.libcalignmentfile cimport AlignmentFile, PileupColumn
from pysam.libchtslib cimport *

cdef class BSSeqPileupRead(PileupRead):
    cdef int32_t _bsseq_strand_ind
    cdef uint32_t _meth_status_flag
    cdef str _observed_watson_base
    cpdef get_meth_status_at_pileup_pos(self, str watson_ref_base)
    # TODO: type correct?
    cdef uint32_t _baseq_at_pos
    cdef uint8_t _overlap_flag
    cdef uint8_t _trimm_flag
    cdef uint8_t _qc_fail_flag


cdef inline make_bsseq_pileup_read(bam_pileup1_t * src, AlignmentFile alignment_file)

cdef inline get_bsseq_strand_index(uint32_t flag)
