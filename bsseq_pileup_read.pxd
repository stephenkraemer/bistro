from pysam.libcalignedsegment cimport AlignedSegment, PileupRead, makeAlignedSegment
from pysam.libcalignmentfile cimport AlignmentFile, PileupColumn
from pysam.libchtslib cimport *

cdef class BSSeqPileupRead(PileupRead):
    cdef readonly uint32_t bs_seq_strand_flag
    cdef readonly uint32_t meth_status_flag
    cdef readonly str observed_watson_base
    cpdef get_meth_status_at_pileup_pos(self, str watson_ref_base)
    # TODO: type correct?
    cdef public uint32_t baseq_at_pos
    cdef public uint8_t overlap_flag
    cdef public uint8_t trimm_flag


cdef inline make_bsseq_pileup_read(bam_pileup1_t * src, AlignmentFile alignment_file)

