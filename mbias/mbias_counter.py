import numpy as np

# Methylation status flags
IS_NA = 16
IS_METHYLATED = 8
IS_UNMETHYLATED = 4
IS_SNP = 2
IS_REF = 1

# BS-Seq strand flags
W_BC = 96
C_BC = 80
W_BC_RV = 144
C_BC_RV = 160
MATE_AND_DIR_BITS = 240

# Indices
C_BC_IND = 0
C_BC_RV_IND = 2
W_BC_IND = 4
W_BC_RV_IND = 6


class MbiasCounter:
    def __init__(self, max_read_length, min_phred_score, max_flen_considered_for_trimming):
        self.counter = np.zeros([8, max_flen_considered_for_trimming + 1, max_read_length + 1], dtype='i4')
        self.phred_score_threshold = min_phred_score
        self.max_flen_considered_for_trimming = max_flen_considered_for_trimming

    def update(self, motif_pileups, index_position):
        watson_motif_seq = index_position.watson_motif
        for motif_base, pileup_reads in zip(watson_motif_seq, motif_pileups):
            if motif_base in ['C', 'G']:
                for pileup_read in pileup_reads:
                    # pileup_read: mqc.bsseq_pileup_read.BSSeqPileupRead
                    if (pileup_read.qc_fail_flag
                            or pileup_read.overlap_flag
                            or pileup_read.trimm_flag):
                        continue

                    # TODO: tlen should not return lower number than number of bases in read
                    # TODO: note thoughts on using the tlen field
                    tlen = abs(pileup_read.alignment.template_length)
                    if tlen > self.max_flen_considered_for_trimming:
                        tlen = self.max_flen_considered_for_trimming

                    pos_in_read = pileup_read.query_position
                    meth_status_flag = pileup_read.get_meth_status_at_pileup_pos(motif_base)

                    if pileup_read.bs_seq_strand_flag == C_BC:
                        if meth_status_flag & IS_METHYLATED:
                            strand_and_meth_status_based_index = C_BC_IND
                        else:
                            strand_and_meth_status_based_index = C_BC_IND + 1

                    elif pileup_read.bs_seq_strand_flag == C_BC_RV:
                        if meth_status_flag & IS_METHYLATED:
                            strand_and_meth_status_based_index = C_BC_RV_IND
                        else:
                            strand_and_meth_status_based_index = C_BC_RV_IND + 1

                    elif pileup_read.bs_seq_strand_flag == W_BC:
                        if meth_status_flag & IS_METHYLATED:
                            strand_and_meth_status_based_index = W_BC_IND
                        else:
                            strand_and_meth_status_based_index = W_BC_IND + 1

                    elif pileup_read.bs_seq_strand_flag == W_BC_RV:
                        if meth_status_flag & IS_METHYLATED:
                            strand_and_meth_status_based_index = W_BC_RV_IND
                        else:
                            strand_and_meth_status_based_index = W_BC_RV_IND + 1

                    else:
                        continue

                    self.counter[strand_and_meth_status_based_index][tlen][pos_in_read] += 1
