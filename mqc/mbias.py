import os
import numpy as np
import pandas as pd
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from abc import ABCMeta, abstractmethod
from typing import Dict
from itertools import product
import pickle

import pandas as pd
idx = pd.IndexSlice

import matplotlib
from matplotlib.axes import Axes # for autocompletion in pycharm
from matplotlib.figure import Figure  # for autocompletion in pycharm
matplotlib.use('Agg') # import before pyplot import!
import matplotlib.pyplot as plt
import seaborn as sns

from mqc.pileup.bsseq_pileup_read import BSSeqPileupRead
from mqc.pileup.pileup import MotifPileup
from mqc.visitors import Counter
from mqc.utils import convert_array_to_df
import mqc.flag_and_index_values as mfl

b_inds = mfl.bsseq_strand_indices
b_na_ind = mfl.bsseq_strand_na_index
m_flags = mfl.methylation_status_flags


class MbiasCounter(Counter):
    """Count stratified M-bias stats

    *Implementation notes:*
    The Fragment length dimension includes the length 0, so that it can be
    indexed by 1-based values. The read position indexes on the other hand
    are zero-based, for better interaction with the C/cython parts of the
    program.
    """

    def __init__(self, config: Dict):

        max_read_length = config['data_properties']['max_read_length_bp']
        motifs = config['run']['motifs']
        self.max_flen_considered_for_trimming = (
            config['trimming']['max_flen_considered_for_trimming'])
        idx_motif_tuples = enumerate(config['run']['motifs'], start=0)
        self.motif_idx_dict = {motif: i for i, motif in idx_motif_tuples}

        dim_names = ['motif', 'bs_strand', 'flen', 'pos', 'meth_status']
        # Note: 1-based position labels in dataframe, 0-based position indices
        # in array
        dim_levels = [motifs,
                      ['c_bc', 'c_bc_rv', 'w_bc', 'w_bc_rv'],
                      range(0, self.max_flen_considered_for_trimming + 1),
                      range(1, max_read_length+1),
                      ['n_meth', 'n_unmeth']]

        array_shape = [len(motifs),
                       4,  # BSSeq-strands
                       self.max_flen_considered_for_trimming + 1,
                       max_read_length,
                       2] # meth status
        counter_array = np.zeros(array_shape, dtype='i4')

        super().__init__(dim_names=dim_names,
                         dim_levels=dim_levels,
                         counter_array=counter_array)

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

        # import pdb; pdb.set_trace()

        curr_motif = motif_pileup.idx_pos.motif
        curr_motif_idx = self.motif_idx_dict[curr_motif]

        curr_read: BSSeqPileupRead
        for curr_read in motif_pileup.reads:
                # TODO: currently this sorts out any qc_fail, including phred score fails, phred score fails should be kept here
                if (curr_read.qc_fail_flag or
                    curr_read.bsseq_strand_ind == b_na_ind):
                    continue

                meth_status_flag = curr_read.meth_status_flag
                if meth_status_flag == m_flags.is_methylated:
                    meth_status_index = 0
                elif meth_status_flag == m_flags.is_unmethylated:
                    meth_status_index = 1
                else:  # SNP, Ref, NA
                    continue

                tlen = abs(curr_read.alignment.template_length)
                if tlen > self.max_flen_considered_for_trimming:
                    tlen = self.max_flen_considered_for_trimming

                event_class = (curr_motif_idx,
                               curr_read.bsseq_strand_ind,
                               tlen,
                               curr_read.pos_in_read,
                               meth_status_index)

                self.counter_array[event_class] += 1


class MbiasData:
    """Provide basic M-bias stats in dataframe format

    Attributes
    ----------
    mbias_stats_df: pd.DataFrame
        index levels: bsseq_strand flen pos
        column labels: n_meth n_unmeth beta_values
                       n_meth_smoothed, n_unmeth_smoothed, beta_values_smoothed
    """

    def __init__(self, mbias_counter, config):

        self.min_flen_considered_for_methylation_calling = (config["data_properties"]
                                                            ["min_flen_considered_for_methylation_calling"])
        self.mbias_counter = mbias_counter
        self.mbias_stats_array = mbias_counter.counter

        self._adjusted_masked_mbias_stats_df = pd.DataFrame()
        self._minimal_masked_mbias_stats_df = pd.DataFrame()
        self._adjusted_cutting_sites_df = pd.DataFrame()
        self._minimal_masked_df = pd.DataFrame()
        self._adjusted_masked_df = pd.DataFrame()
        self.mbias_stats_df = pd.DataFrame()

        self.max_window_size_for_smoothing = config["trimming"]["max_window_size_for_smoothing"]
        self.required_n_events_for_cutting_site_determination = (config["trimming"]
                                                                 ["required_n_events_for_cutting_site_determination"])
        self.max_read_length_bp = config["data_properties"]["max_read_length_bp"]
        self.max_flen_considered_for_trimming = config["trimming"]["max_flen_considered_for_trimming"]

    def convert_mbias_arr_info_to_df_format(self):
        """Convert M-bias counts in ndarray format into dataframe

        Output:

                               n_meth  n_unmeth   \
        bsseq_strand flen pos
        c_bc         1    1                    1.0                    1.0
                          2                    1.0                    1.0

                          n_meth_smoothed n_unmeth_smoothed
                          NA                           NA
                          NA                           NA

        Flen and pos are 1-based labels
        """
        rows = []
        for bsseq_strand_name, bsseq_strand_ind in b_inds._asdict().items():
            for flen in range(0, self.max_flen_considered_for_trimming):
                for pos in range(0, self.max_read_length_bp):
                    row = dict()
                    row['bsseq_strand'] = bsseq_strand_name
                    row['flen'] = flen + 1
                    row['pos'] = pos + 1
                    row['n_meth'] = self.mbias_stats_array[bsseq_strand_ind, flen, pos, 0]
                    row['n_unmeth'] = self.mbias_stats_array[bsseq_strand_ind, flen, pos, 1]
                    rows.append(row)

        df = (pd.DataFrame(rows)
              .set_index(['bsseq_strand', 'flen', 'pos'])
              .sort_index(axis='index', level=0, sort_remaining=True))

        df['beta_values'] = df['n_meth'] / ( df['n_meth'] + df['n_unmeth'] )
        df['n_total'] = df['n_meth'] + df['n_unmeth']

        df['n_meth_smoothed'] = np.nan
        df['n_unmeth_smoothed'] = np.nan
        df['beta_values_smoothed'] = np.nan

        self.mbias_stats_df = (df.astype(dtype=np.float64, copy=True)
                               .sort_index(axis='index', level=0, sort_remaining=True))

    def add_smoothed_mbias_stats(self):
        # TODO-document
        """Smooth M-bias curves for (bsseq_strand, flen) strata with low coverage

        The curves are smoothed by adding data from adjacent, lower flens for the same bsseq_strand,
        based on the observation that lower flens generally have the same or worse M-bias problems, but
        they are generally not less affected. Therefore, smoothing may lead to an over-estimation,
        but not to an underestimation of M-bias.

        Algorithm:
        the ['n_meth_smoothed', 'n_unmeth_smoothed'] columns initially all have
        np.nan values. where possible (enough coverage before/after smoothing), the NA values are replaced with
        either the original beta value counts (where coverage is high enough) or with smoothed beta values (where
        smoothing is possible according to the config parameters)


        Returns
        ---------------
        The M-bias stats dataframe is modified inplace
        """

        df = self.mbias_stats_df
        cols_to_fill = ['n_meth_smoothed', 'n_unmeth_smoothed']
        # both columns are initialized to np.nan upstream
        meth_event_cols = ['n_meth', 'n_unmeth']
        for curr_bsseq_strand in b_inds._fields:
            # State variable set to true if all following flens will have too low coverage
            minimal_flen_with_enough_coverage_reached = False
            # TODO-learn
            for curr_flen in range(self.max_flen_considered_for_trimming,
                                   self.min_flen_considered_for_methylation_calling - 1, -1):
                if minimal_flen_with_enough_coverage_reached:
                    break
                # get index pointing to all positions in the current (bsseq_strand, flen) stratum
                idx_rows_curr_stratum = pd.IndexSlice[curr_bsseq_strand, curr_flen, :]
                meth_events_curr_stratum_arr = df.loc[idx_rows_curr_stratum, meth_event_cols].values
                total_events = meth_events_curr_stratum_arr.sum().sum()
                if total_events >= self.required_n_events_for_cutting_site_determination:
                    df.loc[idx_rows_curr_stratum, cols_to_fill] = meth_events_curr_stratum_arr
                else:
                    # Counts are insufficient for cutting site determination
                    # -> Try to get 'smoothed' curve with more counts by joining adjacent flens
                    for ind in range(1, self.max_window_size_for_smoothing + 1):
                        next_flen = curr_flen - ind
                        if next_flen == 0:
                            minimal_flen_with_enough_coverage_reached = True
                            break
                        # noinspection PyTypeChecker
                        meth_events_curr_stratum_arr += df.loc[(curr_bsseq_strand, next_flen, slice(None)), meth_event_cols].values
                        total_events = meth_events_curr_stratum_arr.sum().sum()
                        if total_events >= self.required_n_events_for_cutting_site_determination:
                            df.loc[idx_rows_curr_stratum, cols_to_fill] = meth_events_curr_stratum_arr
                            break


        df['beta_values_smoothed'] = (
            df['n_meth_smoothed'] / (df['n_meth_smoothed'] + df['n_unmeth_smoothed']))

        df['n_total_smoothed'] = df['n_meth_smoothed'] + df['n_unmeth_smoothed']


    def get_masked_mbias_df(self, trimming_mode: str, cutting_sites):
        """Get M-bias counts in dataframe format, with positions in M-bias trimming zones set to NA"""

        # TODO: do I want to cache?
        # # Look for cached result first
        # if trimming_mode == 'minimal' and not self._minimal_masked_mbias_stats_df.empty:
        #     return self._minimal_masked_mbias_stats_df
        # elif trimming_mode == 'adjusted' and not self._adjusted_masked_mbias_stats_df.empty:
        #     return self._adjusted_masked_mbias_stats_df

        # We have to compute the masked dataframe
        cutting_sites_df = cutting_sites.get_df()

        if cutting_sites_df.empty:
            warnings.warn('Cutting sites dataframe for mode {} unknown. Return unmasked dataframe'.format(
                trimming_mode
            ))
            return self.mbias_stats_df

        def mask_read_positions_in_trimming_zone(df):
            # Hierarchical index with three levels: bsseq_strand, flen, pos; get strand and flen in current group df
            # The cutting site df start and end positions are the start and end of the good plateau -> they are still ok
            # only the adjacent positions have to be merged
            labels = df.index.unique()[0][0:2]

            start_slice_idx, end_slice_idx = cutting_sites_df.loc[labels, ['start', 'end']]
            plateau_start_label = start_slice_idx + 1
            plateau_end_label = end_slice_idx

            idx = pd.IndexSlice
            df.loc[idx[:, :, 1:(plateau_start_label-1)], :] = np.nan
            df.loc[idx[:, :, (plateau_end_label + 1):], :] = np.nan
            return df

        masked_df = (self.mbias_stats_df
                     .groupby(level=['bsseq_strand', 'flen'], axis='index')
                     .apply(mask_read_positions_in_trimming_zone))

        # TODO: do I want to cache this?
        # TODO: this code is weird!
        # if trimming_mode == 'minimal':
        #     setattr(self, '_minimal_masked_mbias_stats_df', masked_df)
        # elif trimming_mode == 'adjusted':
        #     setattr(self, '_adjusted_masked_mbias_stats_df', masked_df)
        # else:
        #     raise ValueError('Unknown trimming mode')

        return masked_df


class MbiasDataPlotter:
    def __init__(self, mbias_data: MbiasData, config):
        self.mbias_data = mbias_data
        # TODO-refactor: bad way of getting config vars!
        self.max_read_length = config['data_properties']['max_read_length_bp']
        self.distance_between_displayed_flen = config["mbias_plots"]["distance_between_displayed_flen"]
        self.sample_name = config['sample']['name']
        self.max_read_length = config['data_properties']['max_read_length_bp']
        self.flens_to_show = config['trimming']['flens_to_show_in_plots']
        # TODO-hardcoded
        self.ylim = [0.5, 1]


        self.min_mean_raw_data = 20

    def total_plot(self, output_path, adjusted_cutting_sites=None):

        if adjusted_cutting_sites:
            df = self.mbias_data.get_masked_mbias_df(
                    'adjusted', cutting_sites=adjusted_cutting_sites)
        else:
            df = self.mbias_data.mbias_stats_df

        df_total_counts = (
            df.loc[:, ['n_meth', 'n_unmeth']]
            .groupby(level=['bsseq_strand', 'pos'])
            .agg(np.sum)
            .reset_index()
            .assign(beta_value=lambda x: x.n_meth /
                         (x.n_meth + x.n_unmeth)
                    )
        )


        g = (sns.FacetGrid(data=df_total_counts,
                           col='bsseq_strand', col_wrap=2,
                           sharex=True, sharey=True,
                           size=cm_to_inch(8), aspect=2,
                           margin_titles=True)
                .map(plt.plot, 'pos', 'beta_value', linewidth=1, alpha=0.4)
                .set_axis_labels('Position in read', 'Average methylation')
                .set(xlim=[1, self.max_read_length],
                     ylim=self.ylim)
        )

        g.set_titles(col_template = '{col_name}',
                     row_template = '{row_name}')
        g.fig.subplots_adjust(top=0.9)
        g.fig.suptitle(self.sample_name, fontsize=12)

        g.savefig(output_path)

    def flen_strat_plot(self, output_path, cutting_sites=None,
                        plot_smoothed_values=True):


        if plot_smoothed_values:
            beta_value_column_name = 'beta_values_smoothed'
        else:
            beta_value_column_name = 'beta_values'

        # TODO-reformat: trimming mode argument is unnecessary for masked df
        #                computation
        if cutting_sites is None:
            mbias_stats_df = self.mbias_data.mbias_stats_df
        else:
            mbias_stats_df = self.mbias_data.get_masked_mbias_df(
                    trimming_mode='adjusted', cutting_sites=cutting_sites)

        df = (mbias_stats_df.loc[pd.IndexSlice[:, self.flens_to_show, :], :])

        if not plot_smoothed_values:
            # TODO-check: mean of NA?
            # discard (bsseq_strand, flen) strata with too low coverage
            # since smoothing was not applied
            df = (df.groupby(level=['bsseq_strand', 'flen'])
                    .filter(lambda df: df['n_total'].mean() >=
                                       self.min_mean_raw_data))

        df = df.loc[:, beta_value_column_name].dropna()

        g = (sns.FacetGrid(data=df.reset_index(), col='bsseq_strand',
                           col_order='w_bc w_bc_rv c_bc c_bc_rv'.split(),
                           col_wrap=2,
                           hue='flen',
                           xlim=[0, self.max_read_length],
                           ylim=self.ylim)
             .map(plt.plot, 'pos', beta_value_column_name)
             .set_axis_labels('Position in read', 'Average methylation')
             .add_legend())

        g.set_titles(col_template = '{col_name}',
                     row_template = '{row_name}')

        g.fig.subplots_adjust(top=0.9)
        g.fig.suptitle(self.sample_name, fontsize=12)

        g.fig.savefig(output_path)


class CuttingSites(metaclass=ABCMeta):
    # TODO: Make sure that plateau intervals are 1-based labels, not zero-based indices
    @abstractmethod
    def get_array(self):
        pass

    @abstractmethod
    def get_df(self):
        pass


class FixedRelativeCuttingSites(CuttingSites):
    """Calculate and provide M-bias cutting sites in dataframe and array format"""

    def __init__(self, config):

        self._minimal_cutting_sites_array = np.array([])
        self._minimal_cutting_sites_df = pd.DataFrame()

        self.max_read_length_bp = config["data_properties"]["max_read_length_bp"]
        self.max_flen = config["trimming"]["max_flen_considered_for_trimming"]
        self.n_bp_to_discard_at_frag_ends = config["trimming"]["relative_to_fragment_ends_dict"]

        self.dim_levels = [b_inds._fields,
                           range(0, self.max_flen + 1),
                           ['left', 'right']]
        self.dim_names = 'bs_strand flen end'.split()
        self.value_column_name = 'inclusive_bound'
    def get_array(self):
        """
        Given the number of bp to be removed from the 5' and 3' fragment ends,
        define cutting sites expressed as zero-based positions in the read.

        The cutting sites are given as plateau start and end, i.e. as the first
        and last acceptable position in the read.

        Cutting sites are given for every (bsseq_strand, flen) stratum with
        flen \in [0, max_flen]

        The (input) 'relative cutting sites' may look like this:
        {
            'w_bc': [0, 9]
            'c_bc': [0, 9]
            'w_bc_rv': [9, 0]
            'c_bc_rv': [9, 0]
        }

        The, the result should be a multidimensional array, which in table format would look like this
        bsseq_strand      where       flen         pos_in_read
        'w_bc'            'start'     1            10
                                      2            10
                                      ...
                                      max_flen
                          end'        1            90
                                      ...
        'w_bc_rv'         ...         ...
        """

        if not self._minimal_cutting_sites_array.size:
            # Flen is 1-based
            min_trimmsite_arr = np.zeros([
                4, (self.max_flen + 1), 2], dtype=np.int32)

            for bsseq_strand_name, bsseq_strand_index in b_inds._asdict().items():

                # Set start position
                # Start position is the zero-based plateau start in the read
                # the relative_cutting_site_dict gives the number of bases
                # to be removed from the start of the read
                min_trimmsite_arr[bsseq_strand_index, :, 0] = (
                    self.n_bp_to_discard_at_frag_ends[bsseq_strand_name][0])

                # Set last included position (zero-based)
                for curr_flen in range(0, self.max_flen):
                    max_allowed_pos_in_fragment = (
                        curr_flen - self.n_bp_to_discard_at_frag_ends[bsseq_strand_name][1] - 1)
                    max_allow_pos_in_read = (max_allowed_pos_in_fragment
                                             if max_allowed_pos_in_fragment <= (self.max_read_length_bp - 1)
                                             else (self.max_read_length_bp - 1) )
                    min_trimmsite_arr[bsseq_strand_index, curr_flen, 1] = max_allow_pos_in_read

            # cache
            self._minimal_cutting_sites_array = min_trimmsite_arr
            return min_trimmsite_arr

        else:
            # return cached value
            return self._minimal_cutting_sites_array

    def get_df(self):
        if self._minimal_cutting_sites_df.empty:
            self._minimal_cutting_sites_df = convert_array_to_df(
                arr=self.get_array(),
                dim_levels=self.dim_levels,
                dim_names=self.dim_names,
                value_column_name=self.value_column_name
            )
        return self._minimal_cutting_sites_df






def cm_to_inch(cm):
    return cm / 2.54




class AdjustedCuttingSites(CuttingSites):

    def __init__(self, mbias_df, config):
        self.max_flen = config['trimming']['max_flen_considered_for_trimming']
        self._df = self._compute_df(mbias_df, config)
        self._arr = np.array([])

    def get_df(self) -> pd.DataFrame:
        return self._df

    def get_array(self):
        if not self._arr.size:
            self._arr = convert_cutting_sites_df_to_array(self._df,
                                                          self.max_flen)
        return self._arr

    @staticmethod
    def _compute_df(mbias_df, config):
        df =  (mbias_df
               .loc['CG', :]
               .groupby(axis='index',
                        level=['bs_strand', 'flen'])
               .apply(fit_normalvariate_plateau, config))
        df.columns.name = 'cut_end'
        df = df.stack().to_frame('cut_pos')
        df = df.sort_index()
        return df


def convert_cutting_sites_df_to_array(df, max_flen):
    # May be shared between AdjustedCuttingSites and
    # UserCuttingSites in the future
    df = df.reset_index()
    level_name_to_index_mappings = {
        'bs_strand': b_inds._asdict(),
        'cut_end': {'left_cut_end': 0, 'right_cut_end': 1},
    }
    df_int_idx = df.replace(level_name_to_index_mappings)
    arr = np.zeros([4, max_flen + 1, 2])
    arr[df_int_idx['bs_strand'].tolist(),
        df_int_idx['flen'].tolist(),
        df_int_idx['cut_end'].tolist()] = df_int_idx['cut_pos']
    return arr

def fit_normalvariate_plateau(group_df: pd.DataFrame, config) -> pd.Series:
    """Find the longest possible plateau of good quality

    Good quality: std along the plateau below a threshold, ends of the plateau
    not different from other read positions (not implemented yet)

    Algorithm:
    1. Start with longest possible plateau length (full read)
    2. Check if the plateau has good quality
    3. If yes: use the plateau and break
    4. Decrease the plateau length by one
    5. Check all possibilities of finding a plateau of the given length within the read
    6. If there are one or more way of fitting the plateau with good quality: break and return the best solution
    7. If not: go to step 4
    """

    min_perc = config["trimming"]["min_plateau_perc"]
    max_std = config["trimming"]["max_std_within_plateau"]
    min_flen = config['trimming']['min_flen_considered_for_trimming']

    beta_values = group_df['beta_value']

    effective_read_length = len(beta_values)
    if effective_read_length < min_flen:
        return pd.Series([0, 0],
                         index=['left_cut_end', 'right_cut_end'])

    if beta_values.isnull().all():
        return pd.Series([0, 0],
                         index=['left_cut_end', 'right_cut_end'])

    min_plateau_length = int(effective_read_length * min_perc)

    for plateau_length in range(effective_read_length, min_plateau_length - 1, -1):
        max_start_pos = effective_read_length - plateau_length
        std_to_beat = max_std
        best_start = None
        for start_pos in range(0, max_start_pos):
            end_pos = start_pos + plateau_length
            curr_beta_values = beta_values.iloc[start_pos:end_pos]
            curr_std = curr_beta_values.std()
            if curr_std < std_to_beat:
                plateau_height = curr_beta_values.mean()
                left_end_bad = ( abs(curr_beta_values[0:4] - plateau_height) > 2 * curr_std).any()
                right_end_bad = ( abs(curr_beta_values[-4:] - plateau_height) > 2 * curr_std).any()
                if not (left_end_bad or right_end_bad):
                    std_to_beat = curr_std
                    best_start = start_pos
                    best_end = end_pos
        if best_start is not None:
            break
    else:
        best_start = 0
        best_end = 0

    return pd.Series([best_start, best_end],
                     index=['left_cut_end', 'right_cut_end'])


def cutting_sites_plot(cutting_sites_df, config):
    g = sns.FacetGrid(data=cutting_sites_df.reset_index(),
                      col='bs_strand', col_wrap=2,
                      col_order=['w_bc', 'w_bc_rv', 'c_bc', 'c_bc_rv'],
                      hue='cut_end')
    g.set(xticks=range(0,500, 50))
    g = g.map(plt.plot, 'flen', 'cut_pos', alpha=0.7)
    g.savefig(config['paths']['adj_cutting_sites_plot'])


def compute_mbias_stats_df(mbias_counts_df):
    # discard positions beyond current flen
    # unstack n_meth, n_unmeth
    # calculate beta_value
    df = (mbias_counts_df['counts']
          .unstack()
          .groupby(level=['motif', 'bs_strand', 'flen'])
          .apply(lambda group_df: group_df.loc[(group_df.name[0], group_df.name[1], group_df.name[2]), :].query(
            "{pos_col_name} <= {curr_flen}".format(
                pos_col_name='pos', curr_flen=group_df.name[2])))
          )
    df['beta_value'] = df['n_meth'] / (df['n_meth'] + df['n_unmeth'])
    df = df[['beta_value', 'n_meth', 'n_unmeth']]
    return df


# # old version
# def mask_mbias_stats_df(df: pd.DataFrame, cutting_sites_df: pd.DataFrame):
#
#     def true_if_needs_masking(bs_strand, flen, pos):
#         left = cutting_sites_df.loc[(bs_strand, flen, 'left_cut_end'), 'cut_pos']
#         right = cutting_sites_df.loc[(bs_strand, flen, 'right_cut_end'), 'cut_pos']
#         if not left <= pos <= right:
#             return True
#         return False
#
#     needs_masking = [true_if_needs_masking(bs_strand, flen, pos)
#                      for bs_strand, flen, pos in zip(
#                         df.index.get_level_values('bs_strand'),
#                         df.index.get_level_values('flen'),
#                         df.index.get_level_values('pos'))]
#
#     res = df.copy(deep=True)
#     res.loc[needs_masking, :] = np.nan
#
#     return res

def mask_mbias_stats_df(df: pd.DataFrame, cutting_sites_df: pd.DataFrame):
    def mask_positions(group_df):
        bs_strand, flen = group_df.name
        left, right = cutting_sites_df.loc[(bs_strand, flen, ['left_cut_end', 'right_cut_end']), 'cut_pos']
        group_df.loc[idx[:, bs_strand, flen, 1:(left-1)], :] = np.nan
        group_df.loc[idx[:, bs_strand, flen, (right+1):], :] = np.nan
        return group_df
    return (df.groupby(level=['bs_strand', 'flen'])
            .apply(mask_positions))


def pos_vs_beta_plots(mbias_stats_dfs_dict, config):
    trunk_path = config['paths']['mbias_plots_trunk']
    aes_mappings = [
        {'row': 'motif', 'col': 'bs_strand', 'hue': None},
        {'row': 'motif', 'col': 'bs_strand', 'hue': 'flen'},
        # {'row': 'motif', 'col': 'bs_strand', 'hue': 'Phred'},
        # {'row': 'flen', 'col': 'bs_strand', 'hue': 'Phred'},
        # {'row': 'Phred', 'col': 'bs_strand', 'hue': 'flen'},
    ]
    for (curr_name, curr_df), curr_aes_mapping in product(
            mbias_stats_dfs_dict.items(), aes_mappings):

        groupby_vars = ([val for val in curr_aes_mapping.values() if val is not None]
                        + ['pos'])
        agg_df = curr_df.groupby(level=groupby_vars).sum()
        agg_df['beta_value'] = agg_df['n_meth']/(agg_df['n_meth'] + agg_df['n_unmeth'])

        plot_df = agg_df.reset_index()

        if 'flen' in groupby_vars:
            row_is_displayed = plot_df['flen'].isin(config['plots']['mbias_flens_to_display'])
            plot_df = plot_df.loc[row_is_displayed, :]
        if 'phred' in groupby_vars:
            row_is_displayed = plot_df['phred'].isin(config['plots']['mbias_phreds_to_display'])
            plot_df = plot_df.loc[row_is_displayed, :]

        p = (sns.FacetGrid(plot_df, **curr_aes_mapping)
             .map(plt.plot, 'pos', 'beta_value')
             .add_legend())
        strat_name = '_'.join([f"{aes}-{var}"
                               for aes, var in curr_aes_mapping.items()])
        p.savefig(f"{trunk_path}_{curr_name}_{strat_name}.png")



def freq_plot_per_motif(mbias_stats_dfs_dict, config):

    trunk_path = config['paths']['cg_occurence_plot_trunk']

    aes_mappings = [
        {'row': 'bs_strand', 'col': 'meth_status', 'hue': None},
        {'row': 'bs_strand', 'col': 'meth_status', 'hue': 'flen'},
        # {'row': 'bs_strand', 'col': 'meth_status', 'hue': 'Phred'},
    ]

    for (curr_name, curr_df), curr_aes_mapping in product(
            mbias_stats_dfs_dict.items(), aes_mappings):

        curr_df = (curr_df[['n_meth', 'n_unmeth']]
                   .stack()
                   .to_frame('counts'))

        for curr_motif, curr_motif_df in curr_df.groupby(level='motif'):

            groupby_vars = (['motif']
                            + [val for val in curr_aes_mapping.values() if val is not None]
                            + ['pos'])
            agg_df = curr_motif_df.groupby(level = groupby_vars).sum()

            plot_df = agg_df.reset_index()

            if 'flen' in groupby_vars:
                row_is_displayed = plot_df['flen'].isin(config['plots']['mbias_flens_to_display'])
                plot_df = plot_df.loc[row_is_displayed, :]
            if 'phred' in groupby_vars:
                row_is_displayed = plot_df['phred'].isin(config['plots']['mbias_phreds_to_display'])
                plot_df = plot_df.loc[row_is_displayed, :]

            p = sns.FacetGrid(plot_df, sharey=False, margin_titles=True,
                              legend_out=True, **curr_aes_mapping)
            p.map(plt.plot, 'pos', 'counts')
            p.add_legend()
            p.fig.tight_layout()

            strat_name = '_'.join([f"{aes}-{var}" for aes, var in curr_aes_mapping.items()])
            p.savefig(f"{trunk_path}_{curr_name}_{strat_name}_{curr_motif}.png")


def create_mbias_stats_plots(mbias_stats_dfs_dict, config):
    pos_vs_beta_plots(mbias_stats_dfs_dict, config)
    freq_plot_per_motif(mbias_stats_dfs_dict, config)


def analyze_mbias_counts(config):

    mbias_evaluate_paths = config['paths']
    # TODO: create paths subdirs in separate logical unit
    os.makedirs(mbias_evaluate_paths['qc_stats_dir'], exist_ok=True, mode=0o770)

    mbias_counts_df = pd.read_pickle(config['paths']['mbias_counts_p'])

    mbias_stats_df = compute_mbias_stats_df(mbias_counts_df)
    mbias_stats_df.to_pickle(mbias_evaluate_paths['mbias_stats_p'])
    mbias_stats_df.reset_index().to_csv(mbias_evaluate_paths['mbias_stats_tsv'],
                                        header=True, index=False, sep='\t')

    adjusted_cutting_sites = AdjustedCuttingSites(mbias_stats_df, config)
    with open(mbias_evaluate_paths['adjusted_cutting_sites_obj_p'], 'wb') as fobj:
        pickle.dump(adjusted_cutting_sites, fobj)
    adjusted_cutting_sites.get_df().to_pickle(
        mbias_evaluate_paths['adjusted_cutting_sites_df_p'])
    adjusted_cutting_sites.get_df().reset_index().to_csv(
        mbias_evaluate_paths['adjusted_cutting_sites_df_tsv'],
        sep = "\t", header = True, index = False)

    masked_mbias_stats_df = mask_mbias_stats_df(mbias_stats_df,
                                                adjusted_cutting_sites.get_df())
    masked_mbias_stats_df.to_pickle(mbias_evaluate_paths['mbias_stats_masked_p'])
    masked_mbias_stats_df.reset_index().to_csv(mbias_evaluate_paths['mbias_stats_masked_tsv'],
                                               sep='\t', header=True, index=False)

    mbias_stats_dfs_dict = {'full': mbias_stats_df,
                            'trimmed': masked_mbias_stats_df}
    create_mbias_stats_plots(mbias_stats_dfs_dict, config)
    cutting_sites_plot(adjusted_cutting_sites.get_df(), config)
    plt.close('all')
