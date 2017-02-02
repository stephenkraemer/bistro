from pysam.libcalignedsegment cimport PileupColumn
from pysam.libcalignmentfile cimport AlignmentFile
from pysam.libchtslib cimport *

DEF W_BC_FLAG = 96
DEF C_BC_FLAG = 80
DEF W_BC_RV_FLAG = 144
DEF C_BC_RV_FLAG = 160
DEF MATE_AND_DIR_BITS = 240

DEF C_BC_IND = 0
DEF C_BC_RV_IND = 1
DEF W_BC_IND = 2
DEF W_BC_RV_IND = 3
DEF NA_STRAND_IND = -1

DEF IS_NA = 16
DEF IS_METHYLATED = 8
DEF IS_UNMETHYLATED = 4
DEF IS_SNP = 2
DEF IS_REF = 1

cdef dict meth_dict = {W_BC_IND: {'C': {'C': IS_METHYLATED,
                                          'T': IS_UNMETHYLATED,
                                          'G': IS_SNP,
                                          'A': IS_SNP,
                                          'N': IS_NA},
                                    'G': {'C': IS_SNP,
                                          'T': IS_SNP,
                                          'G': IS_REF,
                                          'A': IS_SNP,
                                          'N': IS_NA}},
                             W_BC_RV_IND: {'C': {'C': IS_METHYLATED,
                                             'T': IS_UNMETHYLATED,
                                             'G': IS_SNP,
                                             'A': IS_SNP,
                                             'N': IS_NA},
                                       'G': {'C': IS_SNP,
                                             'T': IS_SNP,
                                             'G': IS_REF,
                                             'A': IS_SNP,
                                             'N': IS_NA}},
                             C_BC_IND: {'G': {'C': IS_SNP,
                                          'T': IS_SNP,
                                          'G': IS_METHYLATED,
                                          'A': IS_UNMETHYLATED,
                                          'N': IS_NA},
                                    'C': {'C': IS_REF,
                                          'T': IS_SNP,
                                          'G': IS_SNP,
                                          'A': IS_SNP,
                                          'N': IS_NA}},
                             C_BC_RV_IND: {'G': {'C': IS_SNP,
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
    # TODO: calculate for every read right away? or at least make property with default caching! (probably better alternative)
    # Note that calculating the methylation status requires knowledge of the watson base!
    cpdef get_meth_status_at_pileup_pos(self, str watson_ref_base):
        if self._is_del or self._is_refskip:
            return IS_NA
        if self.observed_watson_base == 'N':
            return IS_NA
        # E.g. if only first in pair is known, but no alignment, therefore no direct., would be 64
        if self._bsseq_strand_ind == NA_STRAND_IND:
            return IS_NA
        return meth_dict[self.bsseq_strand_ind][watson_ref_base][self.observed_watson_base]

    @property
    def baseq_at_pos(self):
        return self._baseq_at_pos

    @property
    def overlap_flag(self):
        return self._overlap_flag
    @overlap_flag.setter
    def overlap_flag(self, value):
        self._overlap_flag = value

    @property
    def trimm_flag(self):
        return self._trimm_flag
    @trimm_flag.setter
    def trimm_flag(self, value):
        self._trimm_flag = value

    @property
    def bsseq_strand_ind(self):
        return  self._bsseq_strand_ind

    @property
    def meth_status_flag(self):
        return self._meth_status_flag

    @property
    def observed_watson_base(self):
        return self._observed_watson_base

    @property
    def qc_fail_flag(self):
        return self._qc_fail_flag

    @qc_fail_flag.setter
    def qc_fail_flag(self, value):
        self._qc_fail_flag = value

    @property
    def pos_in_read(self):
        """Zero-based position in read

        Currently this value is determined assuming that an aligner such as
        bwa-mem is used, which will write the Watson sequence for all alignments.
        This means that queries corresponding to the Crick strand (C-BC, W-BC-RV)
        are represented by the reverse-complement of their actual sequence.

        Then, the query position returned by htslib is the position in the
        sequence given in the alignment, not the position in the read.
        """

        return self._pos_in_read


cdef inline make_bsseq_pileup_read(bam_pileup1_t * src,
                                   AlignmentFile alignment_file):
    """return a PileupRead object construted from a bam_pileup1_t * object."""
    cdef BSSeqPileupRead dest = BSSeqPileupRead.__new__(BSSeqPileupRead)

    dest._alignment = makeAlignedSegment(src.b, alignment_file)
    dest._indel = src.indel
    dest._level = src.level
    dest._is_del = src.is_del
    dest._is_head = src.is_head
    dest._is_tail = src.is_tail
    dest._is_refskip = src.is_refskip
    dest._overlap_flag = 0
    dest._trimm_flag = 0

    dest._qpos = src.qpos
    # TODO-fixme: what happens if no sequence is given? Then query_length == 0, what about
    #             query sequence?
    dest._observed_watson_base = dest._alignment.query_sequence[dest._qpos]
    dest._baseq_at_pos = dest.alignment.query_qualities[dest._qpos]

    dest._bsseq_strand_ind = get_bsseq_strand_index(dest._alignment._delegate.core.flag)

    # TODO-doublecheck: _pos_in_read computation
    # TODO-document: aligner must write all query sequences as Watson sequences
    """ Set position in read
    This assumes that the aligner always prints the query sequence along the watson strand
    In the consequence, the strands which actually map to the crick strand are written out
    reverse complemented.
    Then, the actual position in the read is not given by the position in the query sequence
    for the C-BC and W-BC-RV strands. Instead, one must substract the query position from the
    read length to get the actual read position in these cases
    """
    # baml1_t.core.l_qseq is 0 if no sequence was given
    query_length = dest._alignment._delegate.core.l_qseq
    if dest._bsseq_strand_ind == C_BC_IND or dest._bsseq_strand_ind == W_BC_RV_IND:
        dest._pos_in_read = query_length - (dest._qpos + 1)
    else:
        dest._pos_in_read = dest._qpos

    return dest

cdef inline get_bsseq_strand_index(uint32_t flag):
    # TODO: are these flag fields always defined, or can they have arbitrary values in some situations?
    # TODO: smarter way of doing these tests
    if flag & C_BC_FLAG == C_BC_FLAG:
        return C_BC_IND
    elif flag & C_BC_RV_FLAG == C_BC_RV_FLAG:
        return C_BC_RV_IND
    elif flag & W_BC_FLAG == W_BC_FLAG:
        return W_BC_IND
    elif flag & W_BC_RV_FLAG == W_BC_RV_FLAG:
        return W_BC_RV_IND
    else:  # BSSeq-strand not defined
        return -1

def pileups(PileupColumn pileup_column):
    pileups = []

    if pileup_column.plp == NULL or pileup_column.plp[0] == NULL:
        raise ValueError("PileupColumn accessed after iterator finished")

    # TODO: warning: there could be problems if pileup_column.n and self.buf are
    # out of sync.
    cdef int x
    for x in range(pileup_column.n_pu):
        pileups.append(make_bsseq_pileup_read(&(pileup_column.plp[0][x]),
                                              pileup_column._alignment_file))
    return pileups
