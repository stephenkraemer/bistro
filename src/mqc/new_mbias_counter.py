"""
own module for new counter and new tests to make it clearer which commits
change exp counter and which fix production counter

put in mqc labbook: add handling of sequence contexts i) contain Ns ii)
smaller than 5bp, because we are at end of chrom sequence We want to i)
accept such positions as not 5bp classifiable ii) still count these
positions in a way that allows us to add them to CG, CHG or CHH context
analyses, which will be performed downstream We do not want global CG,
CHG and CHH counters in the same dimension as the sequence contexts, because
these are different levels of aggregation The solution could be to add
additional "5bp motifs" - NCG - NCHG - NCHH these collect events which
cannot be classified in 5bp motifs, but which can be classified into classic
motifs

put in mqc labbook: counting of CN in mbias stats
"""

from itertools import product
from mqc.pileup.bsseq_pileup_read import BSSeqPileupRead
from mqc.mbias import MotifPileup
from mqc.visitors import Counter
import mqc.flag_and_index_values as mfl
from typing import Dict
from math import floor, ceil

import numpy as np
import pandas as pd

b_inds = mfl.bsseq_strand_indices
b_na_ind = mfl.bsseq_strand_na_index
m_flags = mfl.methylation_status_flags

config = {'trimming': {
    'max_flen_with_single_flen_resolution': 200,
    'max_flen_bin': 500,
    'flen_bin_size': 15,
    'max_phred': 40,
    'phred_bin_size': 5,
}}


class MbiasCounter(Counter):
    """Count stratified M-bias stats

    *Implementation notes:*
    The Fragment length dimension includes the length 0, so that it can be
    indexed by 1-based values. The read position indexes on the other hand
    are zero-based, for better interaction with the C/cython parts of the
    program.
    """

    def __init__(self, config: Dict):

        self.save_stem = config["paths"]["mbias_counts"]
        self.max_read_length = config["data_properties"]["max_read_length_bp"]
        sts = config["stats"]
        self.max_flen = sts["max_flen"]
        self.max_single_flen = sts["max_flen_with_single_flen_resolution"]
        self.flen_bin_size = sts["flen_bin_size"]
        self.max_phred = sts["max_phred"]
        self.phred_bin_size = sts["phred_bin_size"]
        self.seq_context_size = config["stats"]["seq_context_size"]

        dim_names = ["seq_context", "bs_strand", "flen",
                     "phred", "pos", "meth_status"]

        flen_breaks = (
            list(range(0, self.max_single_flen + 1)) +
            list(range(self.max_single_flen + 1,
                       self.max_flen + 2,
                       self.flen_bin_size)))
        if flen_breaks[-1] != self.max_flen + 1:
            flen_breaks.append(self.max_flen + 1)

        flen_intervals = (pd.IntervalIndex
                          .from_breaks(flen_breaks, closed="left")
                          .values.tolist())
        self.last_flen_bin_idx = len(flen_intervals) - 1

        phred_breaks = list(range(0, self.max_phred + 2, self.phred_bin_size))
        if phred_breaks[-1] != self.max_phred + 1:
            phred_breaks.append(self.max_phred + 1)
        phred_intervals = (pd.IntervalIndex
                           .from_breaks(phred_breaks, closed="left")
                           .values.tolist())
        self.last_phred_bin_idx = len(phred_intervals) - 1

        self.seq_ctx_idx_dict = get_sequence_context_to_array_index_table(
            self.seq_context_size)

        # Note: 1-based position labels in dataframe, 0-based position indices
        # in array
        ordered_seq_contexts = [seq_context
                                for seq_context, idx in sorted(
                self.seq_ctx_idx_dict.items(),
                key=lambda x: x[1])]
        dim_levels = [ordered_seq_contexts,
                      ['c_bc', 'c_bc_rv', 'w_bc', 'w_bc_rv'],
                      flen_intervals,
                      phred_intervals,
                      range(1, self.max_read_length + 1),
                      ['n_meth', 'n_unmeth']]

        array_shape = [len(self.seq_ctx_idx_dict),
                       4,  # BSSeq-strands
                       len(flen_intervals),
                       len(phred_intervals),
                       self.max_read_length,
                       2]  # meth status
        counter_array = np.zeros(array_shape, dtype='i4')

        super().__init__(dim_names=dim_names,
                         dim_levels=dim_levels,
                         counter_array=counter_array,
                         save_stem=self.save_stem)

    def process(self, motif_pileup: MotifPileup):
        """Extract M-bias stats from MotifPileup

        Reads are discarded if

        - they have a qc_fail_flag
        - their bsseq strand could not be determined
        - they have methylation calling status: NA, SNP or Ref

        Stratified by:

        - motif
        - BSSeq-strand
        - fragment length
        - position in read
        - methylation status
        """

        # Can't handle Ns and too short motifs
        try:
            seq_ctx_idx = self.seq_ctx_idx_dict[
                motif_pileup.idx_pos.seq_context]
        except KeyError:
            return  # can't count this MotifPileup

        curr_read: BSSeqPileupRead
        for curr_read in motif_pileup.reads:
            # TODO: currently this sorts out any qc_fail, including phred
            # score fails, phred score fails should be kept here
            if (curr_read.qc_fail_flag
                or curr_read.bsseq_strand_ind == b_na_ind):
                continue

            meth_status_flag = curr_read.meth_status_flag
            if meth_status_flag == m_flags.is_methylated:
                meth_status_index = 0
            elif meth_status_flag == m_flags.is_unmethylated:
                meth_status_index = 1
            else:  # SNP, Ref, NA
                continue

            tlen = abs(curr_read.alignment.template_length)
            if tlen < self.max_single_flen:
                tlen_idx = tlen
            elif tlen > self.max_flen:
                tlen_idx = self.last_flen_bin_idx
            else:
                tlen_idx = (self.max_single_flen
                            + ceil((tlen - self.max_single_flen)
                                   / self.flen_bin_size))

            phred = curr_read.baseq_at_pos
            if phred > self.max_phred:
                phred_idx = self.last_phred_bin_idx
            else:
                phred_idx = floor(phred / self.phred_bin_size)

            event_class = (seq_ctx_idx,
                           curr_read.bsseq_strand_ind,
                           tlen_idx,
                           phred_idx,
                           curr_read.pos_in_read,
                           meth_status_index)

            self.counter_array[event_class] += 1


def get_sequence_context_to_array_index_table(motif_size: int):
    if motif_size % 2 != 1:
        raise ValueError("Motif size must be an uneven number")

    all_bases = ['C', 'G', 'T', 'A']
    three_letter_bases = ['C', 'G', 'W']

    n_bp_per_side = (motif_size - 1) // 2
    binned_bases_set = ([three_letter_bases] * n_bp_per_side
                        + [['C']] + [three_letter_bases] * n_bp_per_side)

    # note that indicies are given in alphabetical sorting order
    all_binned_motifs = sorted([''.join(motif)
                                for motif in product(*binned_bases_set)])

    binned_motif_to_idx_mapping = {motif: i
                                   for i, motif in
                                   enumerate(all_binned_motifs)}

    l2 = [all_bases] * n_bp_per_side + [['C']] + [all_bases] * n_bp_per_side
    all_5bp_motifs = [''.join(motif) for motif in product(*l2)]

    _5bp_to_three_letter_motif_index_mapping = {
        motif: binned_motif_to_idx_mapping[
            motif.translate(str.maketrans('CGTA', 'CGWW'))]
        for motif in all_5bp_motifs}

    return _5bp_to_three_letter_motif_index_mapping


def map_seq_ctx_to_motif(seq_ctx, use_classical=True):
    """Map sequence context strings containing [ACGTW] to motifs

    Motifs may be classical: [CG, CHG, CHH] or extended (composed of C,G,W)
    # TODO: Ns? at the end of chroms?
    """

    middle_idx = len(seq_ctx) // 2
    if seq_ctx[middle_idx:(middle_idx + 2)] == 'CG':
        return 'CG'

    if use_classical:
        base_mapping = str.maketrans('ACTW', 'HHHH')  # G unchanged
    else:
        base_mapping = str.maketrans('AT', 'WW')  # CGW unchanged

    seq_suffix = seq_ctx[(middle_idx + 1):(middle_idx + 3)]
    motif_suffix = seq_suffix.translate(base_mapping)

    motif = 'C' + motif_suffix
    return motif
