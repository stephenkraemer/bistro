from mqc.visitors import Visitor
from mqc.pileup.pileup import MotifPileup
import mqc.flag_and_index_values as mfl
import numpy as np

b_inds = mfl.bsseq_strand_indices
b_na_ind = mfl.bsseq_strand_na_index
m_flags = mfl.methylation_status_flags
bstrat = mfl.beta_value_stratum_indices
mcall_ind = mfl.meth_calling_indices

class MethCaller(Visitor):
    def __init__(self):
        pass
    def process(self, motif_pileup: MotifPileup):

        n_meth = 0
        n_total = 0
        for curr_read in motif_pileup.reads:
            if (curr_read.qc_fail_flag
                or curr_read.trimm_flag
                or curr_read.overlap_flag):
                continue

            meth_status_flag = curr_read.meth_status_flag
            if meth_status_flag == m_flags.is_methylated:
                n_meth += 1
                n_total += 1
            elif meth_status_flag == m_flags.is_unmethylated:
                n_total += 1
            else:  # SNP, Ref, NA
                continue

        motif_pileup.n_meth = n_meth
        motif_pileup.n_total = n_total

        try:
            motif_pileup.beta_value = n_meth / n_total
        except ZeroDivisionError:
            motif_pileup.beta_value = np.nan


class StratifiedMethCaller(Visitor):
    def __init__(self):
        pass
    def process(self, motif_pileup: MotifPileup):

        meth_counts_arr = np.zeros((7, 2), dtype='uint16')

        for curr_read in motif_pileup.reads:
            if (curr_read.qc_fail_flag
                or curr_read.trimm_flag):
                continue

            #TODO: Write acceptance test for StratifiedMethCaller + StratifiedBetaCounter

            mate_idx = bstrat.mate1 if curr_read.alignment.is_read1 else bstrat.mate2
            # mate indices are 4,5
            # all reads stratum has index 6

            meth_status_flag = curr_read.meth_status_flag
            if meth_status_flag == m_flags.is_methylated:
                meth_counts_arr[curr_read.bsseq_strand_ind, :] += 1
                meth_counts_arr[mate_idx, :] += 1
                if not curr_read.overlap_flag:
                    meth_counts_arr[bstrat.all, :] += 1
            elif meth_status_flag == m_flags.is_unmethylated:
                meth_counts_arr[curr_read.bsseq_strand_ind, mcall_ind.n_total] += 1
                meth_counts_arr[mate_idx, mcall_ind.n_total] += 1
                if not curr_read.overlap_flag:
                    meth_counts_arr[bstrat.all, mcall_ind.n_total] += 1
            else:  # SNP, Ref, NA
                continue

        with np.errstate(divide='ignore', invalid='ignore'):
            motif_pileup.strat_beta_arr = np.divide(meth_counts_arr[:, mcall_ind.n_meth], meth_counts_arr[:, mcall_ind.n_total])
            motif_pileup.strat_beta_arr[~ np.isfinite(motif_pileup.strat_beta_arr)] = np.nan

        motif_pileup.n_meth = meth_counts_arr[bstrat.all, mcall_ind.n_meth]
        motif_pileup.n_total = meth_counts_arr[bstrat.all, mcall_ind.n_total]
        motif_pileup.beta_value = motif_pileup.strat_beta_arr[bstrat.all]
        motif_pileup.meth_counts_arr = meth_counts_arr
