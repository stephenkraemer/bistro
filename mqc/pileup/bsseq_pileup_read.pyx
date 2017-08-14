from pysam.libcalignedsegment cimport PileupColumn
from pysam.libcalignmentfile cimport AlignmentFile
from pysam.libchtslib cimport *

from pysam.libcutils cimport force_bytes, force_str, \
    charptr_to_str, charptr_to_bytes

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

DEF QC_FAIL_SAM_FLAG = 1
DEF QC_FAIL_PHRED_SCORE = 2
DEF QC_FAIL_MAPQ = 4
DEF QC_FAIL_MISSING_INFO = 8
DEF QC_FAIL_SOFTCLIPPED = 16
DEF QC_FAIL_POS_IN_READ_EXCEEDS_TLEN = 32
DEF QC_FAIL_OVERLAP = 64

cdef uint32_t flag_proper_pair_mask = 3
cdef uint32_t flag_qc_fail_mask = 3852

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


cdef class BSSeqPileupRead(PileupRead):

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
                                   AlignmentFile alignment_file,
                                   str expected_watson_base):
    """return a PileupRead object construted from a bam_pileup1_t * object.

    Implementation notes:
    - I got a segmentation fault when I made filling up the dest properties
    conditional on passing the qc tests. Ideally, these properties would only
    be accessed when the qc_fail_flag is 0, but I haven't checked where
    the segfault is coming from yet. For now, I just fill up all essential
    properties regardless of the quality checks
    - The calculation of the position in the read assumes that the aligner always prints the query sequence along the watson strand
    In the consequence, the strands which actually map to the crick strand are written out
    reverse complemented.
    Then, the actual position in the read is not given by the position in the query sequence
    for the C-BC and W-BC-RV strands. Instead, one must substract the query position from the
    read length to get the actual read position in these cases
    - The calculation of the position in the read assumes that softclips at the beginning of the read
      are not part of the insert - to be concordant with the SAM TLEN definition. Other implementations
      will be added in the future
    """

    # TODO: what happens if no sequence is found?
    #       Then query_length == 0, what about query sequence?
    # TODO-doublecheck: _pos_in_read computation
    # TODO-document: aligner must write all query sequences as Watson sequences
    # TODO-fixme: baml1_t.core.l_qseq is 0 if no sequence was given

    cdef BSSeqPileupRead dest = BSSeqPileupRead.__new__(BSSeqPileupRead)
    cdef bam1_t * bam_line = src.b

    # Fill basic properties
    dest._alignment = makeAlignedSegment(bam_line, alignment_file)
    dest._indel = src.indel
    dest._level = src.level
    dest._is_del = src.is_del
    dest._is_head = src.is_head
    dest._is_tail = src.is_tail
    dest._is_refskip = src.is_refskip
    dest._overlap_flag = 0
    dest._trimm_flag = 0
    dest._qc_fail_flag = 0
    dest._qpos = src.qpos
    dest._expected_watson_base = expected_watson_base

    dest._observed_watson_base = dest._alignment.query_sequence[dest._qpos]
    dest._baseq_at_pos = dest.alignment.query_qualities[dest._qpos]
    dest._bsseq_strand_ind = get_bsseq_strand_index(dest._alignment._delegate.core.flag)


    cdef uint8_t * qualities = pysam_bam_get_qual(bam_line)
    cdef uint32_t * cigar_p = pysam_bam_get_cigar(src.b)
    # Discard reads without cigar info - we can't check for softclips, which
    # are necessary for calculation of the read position and to check that
    # we are not in a softclipped region
    if (bam_line.core.l_qseq == 0
        or qualities[0] == 0xff
        or cigar_p == NULL):
        dest._qc_fail_flag |= QC_FAIL_MISSING_INFO
        return dest


    cdef uint16_t flag = dest._alignment._delegate.core.flag
    if (flag & flag_proper_pair_mask != 3 or flag & flag_qc_fail_mask != 0):
        dest._qc_fail_flag |= QC_FAIL_SAM_FLAG
        return dest

    # add meth status
    if (dest._is_del
        or dest._is_refskip
        or dest._observed_watson_base == 'N'
        or dest._bsseq_strand_ind == NA_STRAND_IND):
        dest._meth_status_flag = IS_NA
    else:
        dest._meth_status_flag = (meth_dict[dest._bsseq_strand_ind]
                                  [dest._expected_watson_base][dest._observed_watson_base])


    # Check whether qpos is in a softclip region, set qc_fail_flag if so
    cdef uint32_t n_cig_ops = pysam_get_n_cigar(src.b)
    cdef int32_t query_length = dest._alignment._delegate.core.l_qseq
    cdef int op, k


    cdef softclips_start = 0
    k = 0
    op = cigar_p[k] & BAM_CIGAR_MASK
    if op == BAM_CHARD_CLIP:
        op = cigar_p[1] & BAM_CIGAR_MASK
        k = 1
    if op == BAM_CSOFT_CLIP:
        # there is a softclip at the beginning of the read!
        softclips_start = cigar_p[k] >> BAM_CIGAR_SHIFT
        # print('Found softclip at start: ', softclips_start)
        # print('qpos: ', dest._qpos)

    cdef softclips_end = 0
    k = n_cig_ops - 1
    op = cigar_p[k] & BAM_CIGAR_MASK
    if op == BAM_CHARD_CLIP:
        k -= 1
        op = cigar_p[k] & BAM_CIGAR_MASK
    if op == BAM_CSOFT_CLIP:
        # there is a softclip at the beginning of the read!
        softclips_end = cigar_p[k] >> BAM_CIGAR_SHIFT
        # print('Found softclip at end: ', softclips_end)
        # print('qpos: ', dest._qpos)

    if not (softclips_start < dest._qpos + 1 < query_length - softclips_end + 1):
        dest._qc_fail_flag |= QC_FAIL_SOFTCLIPPED
        return dest

    # Determine read position under the assumption that softclipped bases
    # at the beginning of the read are not part of the original fragment
    # (which is the assumption underlying the SAM TLEN info, see docstring)
    """
    # print('qpos', dest._qpos)
    # print('softclips end', softclips_end)
    # print('softclips start', softclips_end)
    # print('qlen', query_length)
    # print('strand', dest._bsseq_strand_ind)
    """
    if dest._bsseq_strand_ind == C_BC_IND or dest._bsseq_strand_ind == W_BC_RV_IND:
        # print('adjusting read pos')
        dest._pos_in_read = query_length - (dest._qpos + 1) - softclips_end
    else:  # W_BC or C_BC_RV
        # print('adjusting read pos')
        dest._pos_in_read = dest._qpos - softclips_start

    cdef int i
    if dest._pos_in_read + 1 > abs(src.b.core.isize):
        """
        print('position in read greater than tlen')
        print('tlen ', src.b.core.isize)
        print('pos: ', dest._pos_in_read)
        print('softclips end', softclips_end)
        print('softclips start', softclips_start)
        print('strand', dest._bsseq_strand_ind)
        print('query length: ', query_length)
        print('qpos: ', dest._qpos)

        if pysam_get_l_qname(bam_line) == 0:
            qname = ''
        else:
            qname = charptr_to_str(<char *>pysam_bam_get_qname(bam_line))
        print('qname: ', qname)
        print('pos: ', bam_line.core.pos)
        print('cigar ops:')
        for i in range(0, n_cig_ops):
            print(cigar_p[i] & BAM_CIGAR_MASK)
            print(cigar_p[i] >> BAM_CIGAR_SHIFT)
        """
        dest._qc_fail_flag |= QC_FAIL_POS_IN_READ_EXCEEDS_TLEN

    # print('softclips start', softclips_start)
    return dest


cdef inline get_bsseq_strand_index(int32_t flag):
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
        return NA_STRAND_IND

def pileups(PileupColumn pileup_column, str expected_watson_base):
    pileups = []

    if pileup_column.plp == NULL or pileup_column.plp[0] == NULL:
        raise ValueError("PileupColumn accessed after iterator finished")

    # TODO: warning: there could be problems if pileup_column.n and self.buf are
    # out of sync.
    cdef int x
    for x in range(pileup_column.n_pu):
        pileups.append(make_bsseq_pileup_read(&(pileup_column.plp[0][x]),
                                              pileup_column._alignment_file,
                                              expected_watson_base))
    return pileups
