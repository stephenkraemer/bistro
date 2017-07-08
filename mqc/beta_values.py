import numpy as np
import pandas as pd
import mqc
from typing import List
import scipy.ndimage.filters as sfilters

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import itertools as it

MATE1_IDX = 4
MATE2_IDX = 5
TOTAL_IDX = 6

MATE1_STRAND_IDXS = [0, 2]
MATE2_STRAND_IDXS = [1, 3]

METH_IDX = 0
UNMETH_IDX = 1

m_flags = mqc.flag_and_index_values.methylation_status_flags


class StratifiedBetaValueCounter:
    """Extract stratified beta values from MotifPileup AND update counter"""

    def __init__(self, config):
        # Minimal coverage must be >= 1, because we rely on the minimal coverage condition
        # to avoid ZeroDivisionErrors
        config_min_cov = config['beta_value_dist_stats']['min_cov']
        self.min_cov = config_min_cov if config_min_cov >= 1 else 1
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
        # TODO-unidir: overlapping reads are counted in stratified meth events counter
        #              this is compatiable with mate-stratified beta value plots for the directional protocol
        #              is this also true for the undirectional protocol?
        """

        Parameters
        ----------
        motif_pileups
        index_position

        Returns
        -------
        total_counts: np.ndarray(int32)
            Total number of methylated and unmethylated reads, overlapping reads are discarded
        events_counter: np.ndarray(int32)
            Total number of meth. and unmeth. reads, stratified by BS-Seq strand. Overlapping reads are kept(!).


        """
        total_counts = np.zeros([2], dtype='int32')
        events_counter = np.zeros([4, 2], dtype='int32')
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
        """ Calculate the total beta value and update counter if coverage sufficient

        The global beta value is returned for use in methylation calling. Therefore, it is always computed, even if the
         coverage is below the threshold for consideration in the beta value distribution.

        Parameters
        ----------
        total_counts: np.ndarray
            [meth_counts, unmeth_counts]

        Returns
        -------
        total_beta: np.float

        """

        num_events_total = total_counts.sum()
        if num_events_total == 0:
            # Checking num_events_total==0 is cleaner than catching division by zero error
            # due to the way numpy handles division by infinitesimal numbers
            total_beta = np.nan
        else:
            total_beta = total_counts[METH_IDX] / num_events_total

        if num_events_total > self.min_cov:
            self._add_beta_value_to_counter(total_beta, TOTAL_IDX)

        return total_beta

    def _update_stratified_beta_value_counts(self, events_counter):
        """ Calculate beta values stratified by BS-Seq strand and by mate

        Parameters
        ----------
        events_counter: np.ndarray
        """

        # Division by zero cannot occur, because beta values are only computed if self.min_cov is met
        # by definition, self.min_cov > 0, as enforced during class initialisation
        for curr_strand_idx in range(4):
            num_events_curr_strand = events_counter[curr_strand_idx, :].sum()
            if num_events_curr_strand >= self.min_cov:
                beta_value_curr_strand = events_counter[curr_strand_idx, METH_IDX] / num_events_curr_strand
                self._add_beta_value_to_counter(beta_value_curr_strand, curr_strand_idx)

        # calculate mate stratified beta values
        for mate_idx, mate_strand_idxs in [(MATE1_IDX, MATE1_STRAND_IDXS), (MATE2_IDX, MATE2_STRAND_IDXS)]:
            num_events_curr_mate = events_counter[mate_strand_idxs, :].sum().sum()
            if num_events_curr_mate >= self.min_cov:
                beta_value_curr_mate = events_counter[mate_strand_idxs, METH_IDX].sum() / num_events_curr_mate
                self._add_beta_value_to_counter(beta_value_curr_mate, mate_idx)

    def _add_beta_value_to_counter(self, beta_value, counter_row_idx):
        beta_value_idx = np.int16(1000 * np.round(beta_value, 3))
        self.beta_counter[counter_row_idx, beta_value_idx] += 1


class BetaValueData:
    """Compute beta value statistics and provide in series/dataframe format"""

    def __init__(self, config):
        """
        Attributes
        ----------
        self.df: pd.DataFrame()
            index levels: 'region' 'trimming_status' 'bsseq_strand' 'beta_value'
            column labels: 'count'    'Frequency'    'Smoothed_frequency'
        """
        self.sigma_gaussian_kde = config['beta_value_dist_stats']['kde_sigma']
        self.df = pd.DataFrame()
        self.counters: List['mqc.beta_values.StratifiedBetaValueCounter'] = []

    def add_data_from_counter(self, beta_value_counter: 'StratifiedBetaValueCounter',
                              region_str, trimming_status_str):

        self.counters.append(beta_value_counter)

        df = pd.DataFrame(beta_value_counter.beta_counter)

        # TODO: avoid hard coding of strand order and number of strand categories (7 at the moment)
        multi_idx = pd.MultiIndex.from_arrays([[region_str] * 7, [trimming_status_str] * 7,
                                               'c_bc c_bc_rv w_bc w_bc_rv mate1 mate2 total'.split()],
                                              names=['region', 'trimming_status', 'bsseq_strand'])
        df.index = multi_idx

        df.columns = df.columns.get_values() / 1000
        df.columns.name = 'beta_value'
        df = df.stack().to_frame('count')
        self.df = pd.concat([self.df, df], axis='index')
        self.df = self.df.sort_index(level=1, sort_remaining=True)

    # TODO: do this automatically
    def add_smoothed_beta_value_frequencies(self) -> None:
        self._compute_frequencies()
        def gaussian_kde(beta_values_group_ser: pd.Series, sigma):
            beta_values_group_arr = sfilters.gaussian_filter1d(
                    input=beta_values_group_ser, sigma=sigma)
            beta_values_group_ser[:] = beta_values_group_arr
            return beta_values_group_ser

        grouped = self.df.groupby(
                level=['region', 'trimming_status', 'bsseq_strand'])
        self.df['Smoothed_frequency'] = grouped['Frequency'].transform(
                gaussian_kde, sigma=self.sigma_gaussian_kde)

    def _compute_frequencies(self) -> None:
        def get_frequencies_in_stratum(group_df):
            freqs: pd.Series
            freqs = group_df['count'] / group_df['count'].sum()
            df = freqs.to_frame('Frequency')
            return df

        freq_df = self.df.groupby(level=['region', 'trimming_status', 'bsseq_strand']).apply(
                get_frequencies_in_stratum
        )

        self.df['Frequency'] = freq_df['Frequency']

    def __str__(self):
        return self.df.head().__str__()


class BetaValuePlotter:
    def __init__(self, beta_value_data: 'BetaValueData', config):
        self.beta_value_data = beta_value_data
        self.sample_name = config['sample']['name']

    def beta_value_dist_plot(self, out_basepath_abs):
        """Plot total, mate- and strand-stratified beta value distribution plots in a given region"""
        df = self.beta_value_data.df
        idx = pd.IndexSlice

        # TODO-doc: 'total' not 'global'
        bsseq_strand_groups_with_name = [
            (['total'], 'total_beta'),
            (['mate1', 'mate2'], 'mate_beta'),
            ('c_bc c_bc_rv w_bc w_bc_rv'.split(), 'strand_beta')
        ]
        regions = df.index.get_level_values('region').unique()
        plotting_strata = it.product(regions, bsseq_strand_groups_with_name)

        for curr_region, (curr_str_group, curr_str_name) in plotting_strata:
            rows = idx[curr_region, :, curr_str_group, :]
            # both trimming modes and all beta values
            plotting_data = df.loc[rows].reset_index()
            # TODO-doc: legend kwarg only in factorplot, not FacetGrid
            # TODO-snippet
            g = sns.FacetGrid(data=plotting_data,
                              col='bsseq_strand', col_wrap=2,
                              hue='trimming_status',
                              size=cm_to_inch(8), aspect=1,
                              legend_out=True,
                              )

            g = g.map(plt.plot, 'beta_value', 'Smoothed_frequency')

            g.set_axis_labels('Beta value', 'Frequency')

            g.set_titles(
                    col_template = '{col_name}', row_template = '{row_name}')

            g.fig.subplots_adjust(top=0.85)
            g.fig.suptitle(self.sample_name, ha='center')

            g = g.add_legend(title='Trimming')

            out_path = f'{out_basepath_abs}.{curr_region}.{curr_str_name}.png'
            g.savefig(out_path)


def cm_to_inch(cm):
    return cm / 2.54
