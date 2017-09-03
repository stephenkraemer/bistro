from pysam.libcalignedsegment cimport AlignedSegment, PileupRead, makeAlignedSegment
from pysam.libcalignmentfile cimport AlignmentFile, PileupColumn
from pysam.libchtslib cimport *
from cpython cimport array

cdef extern from "htslib_util.h":
    uint32_t * pysam_bam_get_cigar(bam1_t * b)
    uint16_t pysam_get_n_cigar(bam1_t * b)
    char * pysam_bam_get_qname(bam1_t * b)
    uint8_t * pysam_bam_get_seq(bam1_t * b)
    uint8_t * pysam_bam_get_qual(bam1_t * b)
    uint8_t * pysam_bam_get_aux(bam1_t * b)
    int pysam_bam_get_l_aux(bam1_t * b)
    char pysam_bam_seqi(uint8_t * s, int i)
    uint16_t pysam_get_flag(bam1_t * b)
    uint8_t pysam_get_l_qname(bam1_t * b)
    uint8_t pysam_get_qual(bam1_t * b)


cdef class BSSeqPileupRead(PileupRead):
    cdef int32_t _bsseq_strand_ind
    # correct type?
    cdef str _observed_watson_base
    cdef str _expected_watson_base
    # TODO: type correct?
    cdef uint32_t _baseq_at_pos
    cdef uint8_t _overlap_flag
    cdef uint8_t _trimm_flag
    cdef uint8_t _qc_fail_flag
    cdef uint8_t _bad_phred_flag
    cdef int32_t _pos_in_read
    cdef uint8_t _bad_mapq_flag
    cdef uint16_t _meth_status_flag

cdef inline make_bsseq_pileup_read(bam_pileup1_t * src, AlignmentFile alignment_file,
                                   str expected_watson_base)

cdef inline get_bsseq_strand_index(int32_t flag)
