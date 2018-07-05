"""Caller classes: add methylation stats to MotifPileup

Global methylation stats (beta_value, n_meth, n_total): MethCaller

Global methylation stats plus stratified methylation
stats (by BS-Seq strand and mate): StratifiedMethCaller
"""

import numpy as np

import mqc.flag_and_index_values as mfl
from mqc.pileup.pileup import MotifPileup
from mqc.visitors import Visitor

b_inds = mfl.bsseq_strand_indices
b_na_ind = mfl.bsseq_strand_na_index
m_flags = mfl.methylation_status_flags
strat_call_ids = mfl.strat_call_indices
mstat_ind = mfl.meth_status_indices


class MethCaller(Visitor):
    """Add global methylation stats to MotifPileup

    Currently, the attributes beta_value, n_meth and n_unmeth
    are *not* initialized to None during MotifPileup init.

    The MethCaller must be called before any access to these
    attributes.
    """
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
    """ Add stratified meth. event counts and beta values to MotifPileup

    Will also add global attributes n_meth, n_total, beta_value.

    Events from overlapping reads are counted once for the global event counts,
    but both events are considered for the mate- and strand-stratified statistics.

    Added attributes for stratified calls
    -------------------------------------
    meth_counts_arr: np.array(i4)
    strat_beta_arr: np.array(float)
        Will not contain inf values, but may contain np.nan
    """

    def __init__(self):
        pass

    def process(self, motif_pileup: MotifPileup):

        meth_counts_arr = np.zeros((7, 2), dtype='i4')

        for curr_read in motif_pileup.reads:
            if (curr_read.qc_fail_flag
                    or curr_read.trimm_flag):
                continue

            if curr_read.alignment.is_read1:
                mate_idx = strat_call_ids.mate1
            else:
                mate_idx = strat_call_ids.mate2

            meth_status_flag = curr_read.meth_status_flag
            if meth_status_flag == m_flags.is_methylated:
                meth_counts_arr[curr_read.bsseq_strand_ind, :] += 1
                meth_counts_arr[mate_idx, :] += 1
                if not curr_read.overlap_flag:
                    meth_counts_arr[strat_call_ids.all, :] += 1
            elif meth_status_flag == m_flags.is_unmethylated:
                meth_counts_arr[curr_read.bsseq_strand_ind, mstat_ind.n_total] += 1
                meth_counts_arr[mate_idx, mstat_ind.n_total] += 1
                if not curr_read.overlap_flag:
                    meth_counts_arr[strat_call_ids.all, mstat_ind.n_total] += 1
            else:  # SNP, Ref, NA
                continue

        # Division n_meth / n_total is either invalid (0/0) -> np.nan
        # or valid (0 / t or m / t with m, t > 0). It can never be a
        # zero division error m / 0. Therefore, the array will not contain
        # inf values
        with np.errstate(divide='raise', invalid='ignore'):
            motif_pileup.strat_beta_arr = np.divide(meth_counts_arr[:, mstat_ind.n_meth],
                                                    meth_counts_arr[:, mstat_ind.n_total])

        motif_pileup.n_meth = meth_counts_arr[strat_call_ids.all, mstat_ind.n_meth]
        motif_pileup.n_total = meth_counts_arr[strat_call_ids.all, mstat_ind.n_total]
        motif_pileup.beta_value = motif_pileup.strat_beta_arr[strat_call_ids.all]

        motif_pileup.meth_counts_arr = meth_counts_arr
