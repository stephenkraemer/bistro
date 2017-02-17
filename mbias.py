# TODO: make conversion of counter array to pandas dataframe (for adjusted cutting sites)
#  easier by using Pandas built-in capabilities?
# TODO: add conversion of minimal cutting array to dataframe method
# TODO: test plot with minimal cutting site masking

import warnings
from abc import ABCMeta, abstractmethod

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
        """

        Parameters
        ----------
        config: dict


        Attributes
        ----------
        self.counter: np.ndarray
            The 4 dimensions in order: BS-Seq strands(4), flen (max_flen param), read_pos (max_read_length param), methylation status (2)
            The Fragment length dimension includes the length 0, so that it can be indexed by 1-based values

            the read position indexes are zero-based
        """
        self.phred_score_threshold = config["basic_quality_filtering"]["min_phred_score"]
        self.max_flen_considered_for_trimming = config["trimming"]["max_flen_considered_for_trimming"]
        self.max_read_length = config["data_properties"]["max_read_length_bp"]
        self.counter = np.zeros([4, self.max_flen_considered_for_trimming + 1, self.max_read_length, 2], dtype='i4')

    def update(self, motif_pileups, index_position):
        watson_motif_seq = index_position.watson_motif
        for motif_base, pileup_reads in zip(watson_motif_seq, motif_pileups):
            if motif_base in ['C', 'G']:
                for pileup_read in pileup_reads:
                    pileup_read: mqc.bsseq_pileup_read.BSSeqPileupRead
                    if (pileup_read.qc_fail_flag or
                        pileup_read.bsseq_strand_ind == b_na_ind):
                        continue

                    meth_status_flag = pileup_read.get_meth_status_at_pileup_pos(motif_base)
                    if meth_status_flag == 8:
                        meth_status_index = 0
                    elif meth_status_flag == 4:
                        meth_status_index = 1
                    else:  # SNP, Ref base
                        continue

                    # TODO: tlen should not return lower number than number of bases in read
                    # TODO: note thoughts on using the tlen field
                    tlen = abs(pileup_read.alignment.template_length)
                    if tlen > self.max_flen_considered_for_trimming:
                        tlen = self.max_flen_considered_for_trimming

                    # pileup_read: mqc.bsseq_pileup_read.BSSeqPileupRead
                    pos_in_read = pileup_read.pos_in_read

                    self.counter[pileup_read.bsseq_strand_ind][tlen - 1][pos_in_read][meth_status_index] += 1


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


    def get_masked_mbias_df(self, trimming_mode: str, cutting_sites: 'mqc.mbias.CuttingSites'):
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


class MinimalCuttingSites:
    """Calculate and provide M-bias cutting sites in dataframe and array format"""

    def __init__(self, config):
        self.max_read_length_bp = config["data_properties"]["max_read_length_bp"]
        self.max_flen = config["trimming"]["max_flen_considered_for_trimming"]
        self.n_bp_to_discard_at_frag_ends = config["trimming"]["relative_to_fragment_ends_dict"]
        self._minimal_cutting_sites_df = pd.DataFrame()
        self._minimal_cutting_sites_array = np.array([])

    def get_array(self):
        # TODO-algorithm: make sure that the computed absolute cutting sites provide the start and end of the plateau
        # TODO-doublecheck
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
            min_trimmsite_arr = np.zeros([
                4, 2, self.max_flen], dtype=np.int32)

            for bsseq_strand_name, bsseq_strand_index in b_inds._asdict().items():

                # Set start position
                # Start position is the zero-based plateau start in the read
                # the relative_cutting_site_dict gives the number of bases
                # to be removed from the start of the read
                min_trimmsite_arr[bsseq_strand_index, 0, :] = (
                    self.n_bp_to_discard_at_frag_ends[bsseq_strand_name][0])

                # Set end position
                # all positions are zero-based
                # TODO-format: use min_flen parameter?
                for curr_flen in range(0, self.max_flen):
                    max_allowed_pos_in_fragment = (
                        curr_flen - self.n_bp_to_discard_at_frag_ends[bsseq_strand_name][1] + 1)
                    max_allow_pos_in_read = (max_allowed_pos_in_fragment
                                             if max_allowed_pos_in_fragment <= self.max_read_length_bp
                                             else self.max_read_length_bp)
                    min_trimmsite_arr[bsseq_strand_index, 1, curr_flen] = max_allow_pos_in_read

            # cache
            self._minimal_cutting_sites_array = min_trimmsite_arr
            return min_trimmsite_arr

        else:
            # return cached value
            return self._minimal_cutting_sites_array

    def get_df(self):
        if not self._minimal_cutting_sites_df.empty:
            return self._minimal_cutting_sites_df
        else:
            df = counter_array_to_dataframe(
                    self._minimal_cutting_sites_array,
                    labels=[b_inds._fields, ['start', 'end'], [], []],
                    column_names=['bsseq_strand', 'start_or_end', 'flen', 'cutting_site'],
                    increment_integer_labels=True)
            # TODO: unstack dataframe
            # This code is an unfinished draft
            raise NotImplementedError


class AdjustedMbiasCuttingSites:
    """Calculate and provide M-bias cutting sites in dataframe and array format

    Attributes
    ----------


    """

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

        bsseq_strand      where       flen         absolute_cutting_position
        'w_bc'            'start'     0            0
                                      1            10
                                      ...
                                      max_flen
                          'end'       0            0
                          'end'       1            90
                                      ...
        'w_bc_rv'         ...         ...

        The index along the flen dimension corresponds to the flen, therefore, a flen of zero
        is included, but will not be filled.

        The array also contains fragment lengths which are below the threshold for consideration
        in methylation calling. For these fragment lengths, the cutting sites for start and end
        will always be (0,0), which is equivalent to fully discarding these fragments. For performance reasons, this part of the array is not used directly, it is merely included
        so that the indices can correspond directly to the fragment lengths
        """
        if not self._adjusted_cutting_sites_array.size == 0:
            # we have already computed and cached the array
            return self._adjusted_cutting_sites_array
        else:
            df = self.get_df()  # cutting sites dataframe

            """Note that the cutting sites dataframe will always have a 'start' and 'end' column
            by specification, while additional columns are optional, but may be present.
            Therefore, all other columns must be discarded before reformatting the dataframe!"""
            df = df[['start', 'end']]
            df.columns.name = 'start_or_end'

            df_formatted = (df.stack()
                            .to_frame('absolute_cutting_position')
                            .reset_index()
                            # change column order
                            .loc[:, ['bsseq_strand', 'start_or_end', 'flen', 'absolute_cutting_position']]
                            )

            # Replace string dimension labels with integers
            df_formatted.bsseq_strand.replace(b_inds._fields, b_inds, inplace=True)
            df_formatted.start_or_end.replace(['start', 'end'], [0, 1], inplace=True)

            arr = np.zeros([4, 2, self.max_flen_considered_for_trimming + 1])
            # +1 so that indices correspond to 1-based fragment lengths

            for row in df_formatted.itertuples():
                # row is a named tuple, first element is the index
                bsseqstrand_startend_flen_tuple = row[1:4]
                cutting_pos = row[4]
                arr[bsseqstrand_startend_flen_tuple] = cutting_pos
            self._adjusted_cutting_sites_array = arr
            return arr

    def get_df(self) -> pd.DataFrame:
        """ Compute M-bias cutting sites from smoothed M-bias stats (where available)

        Iterates over all (bsseq_strand, flen) strata in the M-bias stats dataframe (must have smoothed beta values).
        If smoothed beta values are available for a given stratum (if the coverage is too low, the stratum will have np.nan instead of beta values),
        cutting sites will be determined using the calling function specified in the config file.

        All calling functions take a dataframe for a given stratum and return a Series with start and end of the plateau, the average height of the plateau as well as arbitrary other fields

        The returned dataframe looks like this (start and end columsn always included, others optional):

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

        beta_values = group_df['beta_values_smoothed']

        # TODO-doublecheck: is the isnull call robust against nullness, i.e. values which are not np.nan but recognized as such?
        #                   perhaps add check that the Series really has numeric values
        # Either all values for a given flen stratum are NAN, or none are NAN
        if beta_values.isnull().all():
            return pd.Series([0, 0, np.nan],
                             index=['start', 'end', 'average_methylation'])

        min_len = config["trimming"]["min_plateau_length"]
        max_std = config["trimming"]["max_std_within_plateau"]
        max_read_length = config['data_properties']['max_read_length_bp']

        flen = group_df.index.get_level_values('flen')[0]
        # TODO-doublecheck
        effective_read_length = (max_read_length if flen >= max_read_length
                                 else flen)


        # TODO-doublecheck: what happens if there is a NA value in the middle of the curve and std returns np.nan?
        for max_start_pos, plateau_length in enumerate(range(effective_read_length, min_len - 1, -1), 1):
            std_to_beat = max_std
            best_start = None
            for start_pos in range(0, max_start_pos):
                end_pos = start_pos + plateau_length
                # TODO-doublecheck
                curr_beta_values = beta_values.iloc[start_pos:end_pos]
                curr_std = curr_beta_values.std()
                if curr_std < std_to_beat:
                    plateau_height = curr_beta_values.mean()
                    left_end_bad = ( abs(curr_beta_values[0:3] - plateau_height) > curr_std).any()
                    right_end_bad = ( abs(curr_beta_values[-3:] - plateau_height) > curr_std).any()
                    if not (left_end_bad or right_end_bad):
                        std_to_beat = curr_std
                        best_start = start_pos
                        best_end = end_pos
                        best_plateau_height = plateau_height
            if best_start is not None:
                break
        else:
            best_start = 0
            best_end = 0
            best_plateau_height = 0

        # TODO-doublecheck

        return pd.Series([best_start, best_end, best_plateau_height],
                         index=['start', 'end', 'average_methylation'])

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

        plateau_start, plateau_end = get_plateau_start_and_end(group_df['beta_values_smoothed'])

        plateau_height = group_df['beta_values_smoothed'][plateau_start:plateau_end].mean()
        return pd.Series([plateau_start, plateau_end, plateau_height], index=['start', 'end', 'average_methylation'])


def cutting_sites_area_plot(mbias_cutting_sites: 'AdjustedMbiasCuttingSites', output_path):
    """
                              start    end   average_methylation
        bsseq_strand  flen
    """
    plot_df = mbias_cutting_sites.get_df().reset_index()
    # TODO-snippet: FacetGrid
    g = sns.FacetGrid(data=plot_df, col='bsseq_strand', col_wrap=2, col_order='w_bc w_bc_rv c_bc c_bc_rv'.split())
    g = g.map(plt.fill_between, 'flen', 'start', 'end')
    g.fig.savefig(output_path)


def cm_to_inch(cm):
    return cm / 2.54

def counter_array_to_dataframe(arr, labels, column_names,
                               increment_integer_labels=True):
    """

    Parameters
    ----------
    arr: multidimensional array
    labels: List[List] one list per dimension of the array.
        List can be empty (-> use integer labels for this dimension)
        or list may contain N labels, where N is the size of the dimension

    Returns
    -------
    pd.DataFrame

    """
    import itertools
    import numpy as np
    import pandas as pd


    indices = [np.arange(x) for x in arr.shape]

    if increment_integer_labels:
        integer_labels = [np.arange(1, x+1) for x in arr.shape]
    else:
        integer_labels = [np.arange(x) for x in arr.shape]

    filled_labels = [labels_for_dim if labels_for_dim else indices_for_dim
                     for labels_for_dim, indices_for_dim in zip(labels, integer_labels)]

    index_label_tuples = [list(zip(indices_for_dim, labels_for_dim))
                          for indices_for_dim, labels_for_dim
                          in zip(indices, filled_labels)]

    rows = []
    for idx_label_tuples in itertools.product(*index_label_tuples):
        idx_tuple = tuple(i for i,l in idx_label_tuples)
        labels = [l for i,l in idx_label_tuples]
        curr_row = labels + [arr[idx_tuple]]
        rows.append(curr_row)

    df = pd.DataFrame(rows, columns=column_names)

    return df
