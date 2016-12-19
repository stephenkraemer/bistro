import numpy as np
import pandas as pd
import mqc

from typing import List

"""
BetaValueCounter
- stratified by mate
- stratified by region

BetaValueDistribution
- beta value stats in dataframe format

BetaValuePlotter

It must be possible to add BetaValueCounters
"""

MATE1_IDX = 4
MATE2_IDX = 5
TOTAL_IDX = 6

MATE1_STRAND_IDXS = [0, 2]
MATE2_STRAND_IDXS = [1, 3]

METH_IDX = 0
UNMETH_IDX = 1

m_flags = mqc.flag_and_index_values.methylation_status_flags


class BetaValueCounter:
    """Extract stratified beta values from MotifPileup AND update counter"""

    def __init__(self, min_cov):
        # Minimal coverage must be >= 1, because we rely on the minimal coverage condition
        # to avoid ZeroDivisionErrors
        self.min_cov = min_cov if min_cov >= 1 else 1
        self.beta_counter = np.zeros([7, 1001])

    def update_and_return_total_beta_value(self, motif_pileups, index_position: 'mqc.IndexPosition'):
        total_counts, event_counter = self._get_stratified_meth_event_counts_for_motif_position(
            motif_pileups, index_position)

        total_beta = self._update_total_beta_value_counts_and_return_beta_value(total_counts)

        self._update_stratified_beta_value_counts(event_counter)

        # return total beta values to avoid duplicate computations for methylation calling
        return total_counts[METH_IDX], total_counts[UNMETH_IDX], total_beta

    @staticmethod
    def _get_stratified_meth_event_counts_for_motif_position(motif_pileups, index_position: 'mqc.IndexPosition'):
        total_counts = np.zeros([2])
        events_counter = np.zeros([4, 2])
        # fill BSSeq-strand stratified events counter
        for watson_base, pileup_reads in zip(index_position.watson_motif, motif_pileups):
            # pileup_reads: List[mqc.BSSeqPileupRead]
            if watson_base in ['C', 'G']:
                for pr in pileup_reads:
                    # pr: mqc.BSSeqPileupRead
                    # TODO: flag bad reads once in own annotation function
                    if pr.qc_fail_flag or pr.trimm_flag:
                        continue
                    meth_status_flag = pr.get_meth_status_at_pileup_pos(watson_base)
                    if meth_status_flag & m_flags.is_methylated:
                        meth_status_idx = METH_IDX
                    elif meth_status_flag & m_flags.is_unmethylated:
                        meth_status_idx = UNMETH_IDX
                    else:
                        continue
                    events_counter[pr.bsseq_strand_ind, meth_status_idx] += 1
                    if not pr.overlap_flag:
                        total_counts[meth_status_idx] += 1

        return total_counts, events_counter

    def _update_total_beta_value_counts_and_return_beta_value(self, total_counts):
        # Global beta value is returned for methylation calling
        # Therefore, it is always computed, even if the coverage is below the threshold for
        # consideration in the beta value distribution
        num_events_total = total_counts.sum()
        try:
            total_beta = total_counts[METH_IDX] / num_events_total
        except ZeroDivisionError:
            total_beta = np.nan

        if num_events_total > self.min_cov:
            self._add_beta_value_to_counter(total_beta, TOTAL_IDX)

        return total_beta

    def _update_stratified_beta_value_counts(self, events_counter):
        for curr_strand_idx in range(4):
            num_events_curr_strand = events_counter[curr_strand_idx, :].sum()
            if num_events_curr_strand > self.min_cov:
                beta_value_curr_strand = events_counter[curr_strand_idx, METH_IDX] / num_events_curr_strand
                self._add_beta_value_to_counter(beta_value_curr_strand, curr_strand_idx)

        # calculate mate stratified beta values
        for mate_idx, mate_strand_idxs in [(MATE1_IDX, MATE1_STRAND_IDXS), (MATE2_IDX, MATE2_STRAND_IDXS)]:
            num_events_curr_mate = events_counter[mate_strand_idxs, :].sum().sum()
            if num_events_curr_mate > self.min_cov:
                beta_value_curr_mate = events_counter[mate_strand_idxs, METH_IDX].sum() / num_events_curr_mate
                self._add_beta_value_to_counter(beta_value_curr_mate, mate_idx)

    def _add_beta_value_to_counter(self, beta_value, counter_row_idx):
        beta_value_idx = np.int16(1000 * np.round(beta_value, 3))
        self.beta_counter[counter_row_idx, beta_value_idx] += 1


class BetaValueData:
    """Provide beta value statistics in series format"""

    def __init__(self):
        self.ser = pd.DataFrame()

    def add_data_from_counter(self, beta_value_counter: 'BetaValueCounter',
                              region_str, trimming_status_str):
        df = pd.DataFrame(beta_value_counter.beta_counter)

        # TODO: avoid hard coding of strand order and number of strand categories (7 at the moment)
        multi_idx = pd.MultiIndex.from_arrays([[region_str] * 7, [trimming_status_str] * 7,
                                               'c_bc c_bc_rv w_bc w_bc_rv mate1 mate2 total'.split()],
                                              names=['Region', 'Trimming status', 'BS-Seq strand'])
        df.index = multi_idx

        df.columns.name = 'Position'
        ser = df.stack()
        ser.name = 'Beta Value'

        self.ser = ser

    def __str__(self):
        return self.ser.head().__str__()


class BetaValuePlotter:
    pass
