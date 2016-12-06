from pysam.libcalignedsegment cimport AlignedSegment, PileupRead, makeAlignedSegment, PileupColumn
from pysam.libcalignmentfile cimport AlignmentFile

DEF W_BC = 96
DEF C_BC = 80
DEF W_BC_RV = 144
DEF C_BC_RV = 160
DEF MATE_AND_DIR_BITS = 240

DEF IS_NA = 16
DEF IS_METHYLATED = 8
DEF IS_UNMETHYLATED = 4
DEF IS_SNP = 2
DEF IS_REF = 1

cdef dict meth_dict = {W_BC: {'C': {'C': IS_METHYLATED,
                                          'T': IS_UNMETHYLATED,
                                          'G': IS_SNP,
                                          'A': IS_SNP,
                                          'N': IS_NA},
                                    'G': {'C': IS_SNP,
                                          'T': IS_SNP,
                                          'G': IS_REF,
                                          'A': IS_SNP,
                                          'N': IS_NA}},
                             W_BC_RV: {'C': {'C': IS_METHYLATED,
                                             'T': IS_UNMETHYLATED,
                                             'G': IS_SNP,
                                             'A': IS_SNP,
                                             'N': IS_NA},
                                       'G': {'C': IS_SNP,
                                             'T': IS_SNP,
                                             'G': IS_REF,
                                             'A': IS_SNP,
                                             'N': IS_NA}},
                             C_BC: {'G': {'C': IS_SNP,
                                          'T': IS_SNP,
                                          'G': IS_METHYLATED,
                                          'A': IS_UNMETHYLATED,
                                          'N': IS_NA},
                                    'C': {'C': IS_REF,
                                          'T': IS_SNP,
                                          'G': IS_SNP,
                                          'A': IS_SNP,
                                          'N': IS_NA}},
                             C_BC_RV: {'G': {'C': IS_SNP,
                                             'T': IS_SNP,
                                             'G': IS_METHYLATED,
                                             'A': IS_UNMETHYLATED,
                                             'N': IS_NA},
                                       'C': {'C': IS_REF,
                                             'T': IS_SNP,
                                             'G': IS_SNP,
                                             'A': IS_SNP,
                                             'N': IS_NA}}}

"""
switch bs_seq_strand {
    case W_BC {
         if watson_ref_base == 'C':
             switch observed_base {
                   case 'C':


}
}
}

switch is the same implementatio as hash table on the assembly code level
"""


cdef class BSSeqPileupRead(PileupRead):
    cpdef get_meth_status_at_pileup_pos(self, str watson_ref_base):
        if self._is_del or self._is_refskip:
            return IS_NA
        if self.observed_watson_base == 'N':
            return IS_NA
        # E.g. if only first in pair is known, but no alignment, therefore no direct., would be 64
        if self.bs_seq_strand_flag not in [W_BC, W_BC_RV, C_BC, C_BC_RV]:
            return IS_NA
        return meth_dict[self.bs_seq_strand_flag][watson_ref_base][self.observed_watson_base]



cdef inline make_bsseq_pileup_read(bam_pileup1_t * src,
                                   AlignmentFile alignment_file):
    """return a PileupRead object construted from a bam_pileup1_t * object."""
    cdef BSSeqPileupRead dest = BSSeqPileupRead.__new__(BSSeqPileupRead)

    dest._alignment = makeAlignedSegment(src.b, alignment_file)
    dest._qpos = src.qpos
    dest._indel = src.indel
    dest._level = src.level
    dest._is_del = src.is_del
    dest._is_head = src.is_head
    dest._is_tail = src.is_tail
    dest._is_refskip = src.is_refskip

    dest.meth_status_flag = IS_NA
    dest.bs_seq_strand_flag = dest._alignment._delegate.core.flag & MATE_AND_DIR_BITS
    #TODO: _qpos may have incorrect value
    dest.observed_watson_base = dest._alignment.query_sequence[dest._qpos]
    dest.baseq_at_pos = dest.alignment.query_qualities[dest._qpos]
    dest.overlap_flag = 0
    dest.trimm_flag = 0
    return dest

def pileups(PileupColumn pileup_column):
    pileups = []

    if pileup_column.plp == NULL or pileup_column.plp[0] == NULL:
        raise ValueError("PileupColumn accessed after iterator finished")

    # TODO: warning: there could be problems if pileup_column.n and self.buf are
    # out of sync.
    cdef int x
    for x from 0 <= x < pileup_column.n_pu:
        pileups.append(make_bsseq_pileup_read(&(pileup_column.plp[0][x]),
                                              pileup_column._alignment_file))
    return pileups
