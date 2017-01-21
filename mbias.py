# TODO: make conversion of counter array to pandas dataframe (for adjusted cutting sites)
#  easier by using Pandas built-in capabilities?
# TODO: add conversion of minimal cutting array to dataframe method
# TODO: test plot with minimal cutting site masking

import warnings
from abc import ABCMeta, abstractmethod
from typing import Tuple

import pandas as pd
import numpy as np
import mqc

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

b_inds = mqc.flag_and_index_values.bsseq_strand_indices
b_na_ind = mqc.flag_and_index_values.bsseq_strand_na_index
m_flags = mqc.flag_and_index_values.methylation_status_flags


class MbiasCounter:
    """Count M-bias stats from a stream of motif pileups and provide as np.ndarray"""

    def __init__(self, config):
        self.phred_score_threshold = config["basic_quality_filtering"]["min_phred_score"]
        self.max_flen_considered_for_trimming = config["trimming"]["max_flen_considered_for_trimming"]
        self.max_read_length = config["data_properties"]["max_read_length_bp"]
        self.counter = np.zeros([4, self.max_flen_considered_for_trimming + 1, self.max_read_length + 1, 2], dtype='i4')

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

                    if meth_status_flag == 8:
                        meth_status_index = 0
                    elif meth_status_flag == 4:
                        meth_status_index = 1
                    else:  # SNP, Ref base
                        continue

                    if pileup_read.bsseq_strand_ind == b_na_ind:
                        continue

                    self.counter[pileup_read.bsseq_strand_ind][tlen][pos_in_read][meth_status_index] += 1


class MbiasData:
    """Provide basic M-bias stats in dataframe format"""

    def __init__(self, mbias_counter, config):

        self.min_flen_considered_for_methylation_calling = (config["data_properties"]
                                                            ["min_flen_considered_for_methylation_calling"])
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

                               meth_events_per_pos  unmeth_events_per_pos   \
        bsseq_strand flen pos
        c_bc         1    1                    1.0                    1.0
                          2                    1.0                    1.0

                          smoothed_meth_events_per_pos smoothed_unmeth_events_per_pos
                          NA                           NA
                          NA                           NA
        """
        rows = []
        for bsseq_strand_name, bsseq_strand_ind in b_inds._asdict().items():
            for flen in range(1, self.max_flen_considered_for_trimming + 1):
                for pos in range(1, self.max_read_length_bp):
                    row = dict()
                    row['bsseq_strand'] = bsseq_strand_name
                    row['flen'] = flen
                    row['pos'] = pos
                    row['meth_events_per_pos'] = self.mbias_stats_array[bsseq_strand_ind, flen, pos, 0]
                    row['unmeth_events_per_pos'] = self.mbias_stats_array[bsseq_strand_ind, flen, pos, 1]
                    rows.append(row)

        df = (pd.DataFrame(rows)
              .set_index(['bsseq_strand', 'flen', 'pos'])
              .sort_index(axis='index', level=0, sort_remaining=True))
        df.loc[:, 'smoothed_meth_events_per_pos'] = np.nan
        df.loc[:, 'smoothed_unmeth_events_per_pos'] = np.nan
        self.mbias_stats_df = (df.astype(dtype=np.float64, copy=True)
                               .sort_index(axis='index', level=0, sort_remaining=True))

    def add_smoothed_mbias_stats(self):
        df = self.mbias_stats_df
        win_cols = ['smoothed_meth_events_per_pos', 'smoothed_unmeth_events_per_pos']
        flen_data_cols = ['meth_events_per_pos', 'unmeth_events_per_pos']
        for bsseq_strand in b_inds._fields:
            minimal_flen_with_enough_coverage_reached = False
            for flen in range(self.max_flen_considered_for_trimming,
                              self.min_flen_considered_for_methylation_calling, -1):
                if minimal_flen_with_enough_coverage_reached:
                    break
                curr_flen_rows = pd.IndexSlice[bsseq_strand, flen, :]
                curr_win_data = df.loc[curr_flen_rows, flen_data_cols].values
                total_events = curr_win_data.sum().sum()
                if total_events < self.required_n_events_for_cutting_site_determination:
                    for ind in range(1, self.max_window_size_for_smoothing):
                        next_flen = flen - ind
                        if next_flen == 0:
                            minimal_flen_with_enough_coverage_reached = True
                            break
                        # noinspection PyTypeChecker
                        curr_win_data += df.loc[(bsseq_strand, next_flen, slice(None)), flen_data_cols].values
                        total_events = curr_win_data.sum().sum()
                        if total_events >= self.required_n_events_for_cutting_site_determination:
                            df.loc[curr_flen_rows, win_cols] = curr_win_data
                            break


    def add_beta_values_from_smoothed_data(self):
        self.mbias_stats_df['smoothed_beta_values'] = self.mbias_stats_df['smoothed_meth_events_per_pos'] / (
            self.mbias_stats_df['smoothed_meth_events_per_pos'] + self.mbias_stats_df['smoothed_unmeth_events_per_pos'])

    def get_masked_mbias_df(self, trimming_mode: str, cutting_sites: 'mqc.mbias.CuttingSites'):
        """Get M-bias counts in dataframe format, with positions in M-bias trimming zones set to NA"""

        # Look for cached result first
        if trimming_mode == 'minimal' and not self._minimal_masked_mbias_stats_df.empty:
            return self._minimal_masked_mbias_stats_df
        elif trimming_mode == 'adjusted' and not self._adjusted_masked_mbias_stats_df.empty:
            return self._adjusted_masked_mbias_stats_df

        # We have to compute the masked dataframe
        cutting_sites_df = cutting_sites.get_df()

        if cutting_sites_df.empty:
            warnings.warn('Cutting sites dataframe for mode {} unknown. Return unmasked dataframe'.format(
                trimming_mode
            ))
            return self.mbias_stats_df

        def mask_read_positions_in_trimming_zone(df):
            # Hierarchical index with three levels: bsseq_strand, flen, pos; get strand and flen in current group df
            labels = df.index.unique()[0][0:2]
            start, end = cutting_sites_df.loc[labels, ['start', 'end']]
            idx = pd.IndexSlice
            df.loc[idx[:, :, 1:(start - 1)], :] = np.nan
            df.loc[idx[:, :, (end + 1):], :] = np.nan
            return df

        masked_df = (self.mbias_stats_df
                     .groupby(level=['bsseq_strand', 'flen'], axis='index')
                     .apply(mask_read_positions_in_trimming_zone))

        if trimming_mode == 'minimal':
            setattr(self, '_minimal_masked_mbias_stats_df', masked_df)
        elif trimming_mode == 'adjusted':
            setattr(self, '_adjusted_masked_mbias_stats_df', masked_df)
        else:
            raise ValueError('Unknown trimming mode')

        return masked_df


class MbiasDataPlotter:
    def __init__(self, mbias_data: MbiasData, config):
        self.mbias_data = mbias_data
        self.distance_between_displayed_flen = config["mbias_plots"]["distance_between_displayed_flen"]

    def flen_strat_plot(self, output_path, cutting_sites=None, trimming_mode=None):

        if cutting_sites is None:
            mbias_stats_df = self.mbias_data.mbias_stats_df
        elif trimming_mode in ['minimal', 'adjusted']:  #
            mbias_stats_df = self.mbias_data.get_masked_mbias_df(trimming_mode,
                                                                 cutting_sites=cutting_sites)
        else:
            raise NotImplementedError('Trimming mode unknown')

        # TODO: parameter for number of fragments (or: all fragments?)
        plotting_data = (mbias_stats_df
                         .loc[pd.IndexSlice[:, range(self.mbias_data.min_flen_considered_for_methylation_calling,
                                                     self.mbias_data.max_flen_considered_for_trimming + 1,
                                                     self.distance_between_displayed_flen), :], 'smoothed_beta_values']
                         .reset_index())

        plotting_data = plotting_data.dropna(axis='index', how='any')
        g = (sns.FacetGrid(data=plotting_data, col='bsseq_strand',
                           col_order='w_bc w_bc_rv c_bc c_bc_rv'.split(),
                           col_wrap=2,
                           hue='flen')
             .map(plt.plot, 'pos', 'smoothed_beta_values')
             .add_legend())
        g.fig.savefig(output_path)


class CuttingSites(metaclass=ABCMeta):
    # TODO: Make sure that plateau intervals are 1-based labels, not zero-based indices
    @abstractmethod
    def get_array(self):
        pass

    @abstractmethod
    def get_df(self):
        pass


class MinimalCuttingSites:
    """Calculate and provide M-bias cutting sites in dataframe and array format"""

    def __init__(self, config):
        self.max_read_length_bp = config["data_properties"]["max_read_length_bp"]
        self.max_flen_considered_for_trimming = config["trimming"]["max_flen_considered_for_trimming"]
        self.relative_cutting_site_dict = config["trimming"]["relative_to_fragment_ends_dict"]
        self._minimal_cutting_sites_df = pd.DataFrame()
        self._minimal_cutting_sites_array = np.array([])

    def get_array(self):
        """
        Given minimal cutting sites (defined relative to the fragment ends),
        define an multidimensional array of absolute cutting sites (at the beginning and end
        of the fragment) for every allowed fragment length on every BS-Seq sequencing strand

        The (input) relative cutting sites may look like this:
        {
            'w_bc': [0, 9]
            'c_bc': [0, 9]
            'w_bc_rv': [9, 0]
            'c_bc_rv': [9, 0]
        }

        The, the result should be a multidimensional array, which in table format would look like this
        bsseq_strand      where       flen         absolute position
        'w_bc'            'start'     1            10
                                      2            10
                                      ...
                                      max_flen
                          end'        1            90
                                      ...
        'w_bc_rv'         ...         ...


        The array contains values for all fragment lengths between 0 and the max. flen, and
        flen below min flen are not set to discard via (0,0)
        """
        if not self._minimal_cutting_sites_array.size:
            min_cutting_sites_array = np.zeros([4, 2, self.max_flen_considered_for_trimming + 1], dtype=np.int32)

            for bsseq_strand_name, bsseq_strand_index in b_inds._asdict().items():
                # Set start position
                min_cutting_sites_array[bsseq_strand_index, 0, :] = self.relative_cutting_site_dict[bsseq_strand_name][
                    0]

                # Set end position
                for curr_flen in range(self.max_flen_considered_for_trimming + 1):
                    max_allowed_pos_in_fragment = curr_flen - self.relative_cutting_site_dict[bsseq_strand_name][1]
                    max_allow_pos_in_read = (max_allowed_pos_in_fragment
                                             if max_allowed_pos_in_fragment <= self.max_read_length_bp
                                             else self.max_read_length_bp)
                    min_cutting_sites_array[bsseq_strand_index, 1, curr_flen] = max_allow_pos_in_read

            self._minimal_cutting_sites_array = min_cutting_sites_array
            return min_cutting_sites_array
        else:
            return self._minimal_cutting_sites_array

    def get_df(self):
        if self._minimal_cutting_sites_df.empty:
            raise NotImplementedError
        else:
            return self._minimal_cutting_sites_df


class AdjustedMbiasCuttingSites:
    """Calculate and provide M-bias cutting sites in dataframe and array format"""

    def __init__(self, mbias_data: MbiasData, calling_mode, config):
        self.mbias_data = mbias_data
        self.calling_mode = calling_mode
        self.config = config
        self.max_read_length_bp = config["data_properties"]["max_read_length_bp"]
        self.max_flen_considered_for_trimming = config["trimming"]["max_flen_considered_for_trimming"]
        self.cutting_site_calling_mode_dict = {
            'standard': self.fit_normalvariate_plateau
        }
        self._adjusted_cutting_sites_array = np.array([])
        self._adjusted_cutting_sites_df = pd.DataFrame()

    def get_array(self):
        """Convert the M-bias cutting sites dataframe to multidim. array

        The result should be a multidimensional array, which in table format would look like this

        bsseq_strand      where       flen         absolute position
        'w_bc'            'start'     1            10
                                      2            10
                                      ...
                                      max_flen
                          end'        1            90
                                      ...
        'w_bc_rv'         ...         ...

        The array also contains fragment lengths which are below the threshold for consideration
        in methylation calling. For these fragment lengths, the cutting sites for start and end
        will always be (0,0), which is equivalent to fully discarding these fragments.

        For performance reasons, this part of the array is not used directly, it is merely included
        so that the indices can correspond directly to the fragment lengths
        """
        if self._adjusted_cutting_sites_array.size == 0:
            if self._adjusted_cutting_sites_df.empty:
                df = self.get_df()
            else:
                df = self._adjusted_cutting_sites_df
            df.columns.name = 'start_or_end'
            df_formatted = (df.stack()
                            .to_frame('absolute_cutting_position')
                            .reset_index()
                            .loc[:, ['bsseq_strand', 'start_or_end', 'flen', 'absolute_cutting_position']]
                            )
            df_formatted.bsseq_strand.replace(b_inds._fields, b_inds, inplace=True)
            df_formatted.start_or_end.replace(['start', 'end'], [0,1], inplace=True)
            arr = np.zeros([4, 2, self.max_flen_considered_for_trimming + 1])
            # +1 so that indices correspond to 1-based fragment lengths
            df_formatted: pd.DataFrame
            for row in df_formatted.itertuples():
                # row is a named tuple, first element is the index
                bsseqstrand_startend_flen_tuple = row[1:4]
                cutting_pos = row[4]
                arr[bsseqstrand_startend_flen_tuple] = cutting_pos
            self._adjusted_cutting_sites_array = arr
            return arr
        else:
            return self._adjusted_cutting_sites_array

    def get_df(self) -> pd.DataFrame:
        """ Compute M-bias cutting sites from smoothed M-bias stats (where available)

        Iterates over all (bsseq_strand, flen) strata in the M-bias stats dataframe (must have smoothed beta values).
        If smoothed beta values are available for a given stratum (if the coverage is too low, the stratum will have np.nan instead of beta values),
        cutting sites will be determined using the calling function specified in the config file.

        All calling functions take a dataframe for a given stratum and return a Series with start and end of the plateau, the average height of the plateau as well as arbitrary other fields

        The returned dataframe looks like this:

                              start    end   average_methylation
        bsseq_strand  flen
        'c_bc'        1       20        81   0.712312

        If a fragment has to be discarded completely (no smoothed beta values available or bad quality), start and end
        are given as (0,0)

        """
        # TODO: is this true:  Note that the start and end of the plateau interval are 1-based labels, not 0-based indices
        if self._adjusted_cutting_sites_df.empty:
            calling_function = self.cutting_site_calling_mode_dict[self.calling_mode]
            cutting_site_df = (self.mbias_data.mbias_stats_df.groupby(axis='index', level=['bsseq_strand', 'flen'])
                               .apply(calling_function, self.config))
            self._adjusted_cutting_sites_df = cutting_site_df
            return cutting_site_df
        else:
            return self._adjusted_cutting_sites_df

    @staticmethod
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

        def get_plateau_start_and_end(beta_values: pd.Series) -> Tuple[int, int]:
            """Get plateau position for a single (bsseq_strand, flen) stratum"""

            # Either all values for a given flen stratum are NAN, or none are NAN
            if np.isnan(beta_values.iloc[0]):
                return 0,0

            min_len = config["trimming"]["min_plateau_length"]
            max_std = config["trimming"]["max_std_within_plateau"]
            read_length = len(beta_values)

            for max_start_pos, plateau_length in enumerate(range(read_length, min_len - 1, -1)):
                std_to_beat = max_std
                best_start = None
                for start_pos in range(0, max_start_pos + 1):
                    end_pos = start_pos + plateau_length
                    curr_std = beta_values[start_pos:end_pos].std()
                    if curr_std < std_to_beat:
                        std_to_beat = curr_std
                        best_start = start_pos
                if best_start is not None:
                    return best_start, best_start + plateau_length

            return 0,0


        plateau_start, plateau_end = get_plateau_start_and_end(group_df['smoothed_beta_values'])

        # TODO: is this integer indexing or label indexing?
        plateau_height = group_df['smoothed_beta_values'][plateau_start:plateau_end].mean()
        return pd.Series([plateau_start, plateau_end, plateau_height], index=['start', 'end', 'average_methylation'])

    @staticmethod
    def max_plateau_length_threshold_std(group_df, config):
        """This method is a backup and will be discarded later"""
        def get_plateau_start_and_end(beta_values: pd.Series):

            min_len = config["trimming"]["min_plateau_length"]
            max_std = config["trimming"]["max_std_within_plateau"]

            read_length = len(beta_values)

            best_interval = (0, 0)
            best_interval_length = 0
            best_interval_std = 100
            for start_pos in range(0, read_length - min_len):
                for end_pos in range(start_pos + min_len, read_length):
                    curr_std = beta_values[start_pos:end_pos].std()
                    if curr_std < max_std:
                        curr_interval_length = end_pos - start_pos
                        if curr_interval_length > best_interval_length:
                            best_interval = (start_pos, end_pos)
                            best_interval_std = curr_std
                            best_interval_length = curr_interval_length
                        elif (curr_interval_length == best_interval_length
                              and curr_std < best_interval_std):
                            best_interval = (start_pos, end_pos)
                            best_interval_std = curr_std
                            best_interval_length = curr_interval_length
            return best_interval

        plateau_start, plateau_end = get_plateau_start_and_end(group_df['smoothed_beta_values'])

        plateau_height = group_df['smoothed_beta_values'][plateau_start:plateau_end].mean()
        return pd.Series([plateau_start, plateau_end, plateau_height], index=['start', 'end', 'average_methylation'])

def cutting_sites_area_plot(mbias_cutting_sites: 'AdjustedMbiasCuttingSites', output_path):
    plot_df = mbias_cutting_sites.get_df().reset_index()
    g = sns.FacetGrid(data=plot_df, col='bsseq_strand', col_wrap=2, col_order='w_bc w_bc_rv c_bc c_bc_rv'.split())
    g = g.map(plt.fill_between, 'flen', 'start', 'end')
    g.fig.savefig(output_path)
