from mqc.visitors import Visitor
from mqc.pileup.pileup import MotifPileup
import mqc.flag_and_index_values as mfl
import numpy as np

b_inds = mfl.bsseq_strand_indices
b_na_ind = mfl.bsseq_strand_na_index
m_flags = mfl.methylation_status_flags

class MethCaller(Visitor):
    def __init__(self):
        pass
    def process(self, motif_pileup: MotifPileup):

        n_meth = 0
        n_unmeth = 0
        for curr_read in motif_pileup.reads:
            if (curr_read.overlap_flag
                or curr_read.qc_fail_flag
                # TODO: not necessary to check the bsseq strand (?)
                #       this should be done implicitely when I read out the
                #       meth status flag. Robust also for future to leave this out?
                or curr_read.bsseq_strand_ind == b_na_ind):
                continue

            meth_status_flag = curr_read.meth_status_flag
            if meth_status_flag == m_flags.is_methylated:
                n_meth += 1
            elif meth_status_flag == m_flags.is_unmethylated:
                n_unmeth += 1
            else:  # SNP, Ref, NA
                continue

        motif_pileup.n_meth = n_meth
        motif_pileup.n_unmeth = n_unmeth

        try:
            motif_pileup.beta_value = n_meth / (n_meth + n_unmeth)
        except ZeroDivisionError:
            motif_pileup.beta_value = np.nan

