import pandas as pd
import numpy as np
import mqc

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

b_flags = mqc.flag_and_index_values.bsseq_strand_flags
b_inds = mqc.flag_and_index_values.bsseq_strand_indices
m_flags = mqc.flag_and_index_values.methylation_status_flags


class MbiasCounter:
    """Count M-bias stats from a stream of motif pileups and provide as np.ndarray"""
    def __init__(self, max_read_length, min_phred_score, max_flen_considered_for_trimming):
        self.counter = np.zeros([4, max_flen_considered_for_trimming + 1, max_read_length + 1, 2], dtype='i4')
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

                    if meth_status_flag == 8:
                        meth_status_index = 0
                    elif meth_status_flag == 4:
                        meth_status_index = 1
                    else:  # SNP, Ref base
                        continue

                    if pileup_read.bs_seq_strand_flag == b_flags.c_bc:
                        bsseq_strand_index = b_inds.c_bc
                    elif pileup_read.bs_seq_strand_flag == b_flags.c_bc_rv:
                        bsseq_strand_index = b_inds.c_bc_rv
                    elif pileup_read.bs_seq_strand_flag == b_flags.w_bc:
                        bsseq_strand_index = b_inds.w_bc
                    elif pileup_read.bs_seq_strand_flag == b_flags.w_bc_rv:
                        bsseq_strand_index = b_inds.w_bc_rv
                    else:
                        continue

                    self.counter[bsseq_strand_index][tlen][pos_in_read][meth_status_index] += 1


class MbiasData:
    """Provide basic M-bias stats in dataframe format"""

    def __init__(self, max_flen_considered_for_trimming, max_read_length_bp, mbias_stats_array,
                 required_n_events_for_cutting_site_determination, max_window_size_for_smoothing):
        self._adjusted_mask_mbias_stats_df = None
        self._minimal_mask_mbias_stats_df = None
        self.max_window_size_for_smoothing = max_window_size_for_smoothing
        self.required_n_events_for_cutting_site_determination = required_n_events_for_cutting_site_determination
        self.mbias_stats_array = mbias_stats_array
        self.max_read_length_bp = max_read_length_bp
        self.max_flen_considered_for_trimming = max_flen_considered_for_trimming
        self.mbias_stats_df = self.convert_mbias_arr_info_to_df_format()
        self.add_smoothed_mbias_stats()
        self.add_beta_values()

    def convert_mbias_arr_info_to_df_format(self) -> pd.DataFrame:
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
        return (df.astype(dtype=np.float64, copy=True)
                .sort_index(axis='index', level=0, sort_remaining=True))

    def add_smoothed_mbias_stats(self):
        df = self.mbias_stats_df
        win_cols = ['smoothed_meth_events_per_pos', 'smoothed_unmeth_events_per_pos']
        flen_data_cols = ['meth_events_per_pos', 'unmeth_events_per_pos']
        for bsseq_strand in b_inds._fields:
            minimal_flen_with_enough_coverage_reached = False
            # TODO: allow all fragment lengths (currently computation takes a while,
            #       so I am not doing this while testing)
            for flen in range(self.max_flen_considered_for_trimming, 30, -1):
                # for flen in range(200, 40, -1):
                if minimal_flen_with_enough_coverage_reached:
                    break
                curr_flen_rows = pd.IndexSlice[bsseq_strand, flen, :]
                curr_win_data = df.loc[curr_flen_rows, flen_data_cols].values
                total_events = curr_win_data.sum().sum()
                if total_events < self.required_n_events_for_cutting_site_determination:
                    # stop if not enough coverage after smoothing ten adjacent flens
                    for ind in range(1, self.max_window_size_for_smoothing):
                        next_flen = flen - ind
                        if next_flen == 0:
                            minimal_flen_with_enough_coverage_reached = True
                            break
                        curr_win_data += df.loc[(bsseq_strand, next_flen, slice(None)), flen_data_cols].values
                        total_events = curr_win_data.sum().sum()
                        if total_events >= self.required_n_events_for_cutting_site_determination:
                            df.loc[curr_flen_rows, win_cols] = curr_win_data
                            break

    def add_beta_values(self):
        self.mbias_stats_df['smoothed_beta_values'] = self.mbias_stats_df['smoothed_meth_events_per_pos'] / (
            self.mbias_stats_df['smoothed_meth_events_per_pos'] + self.mbias_stats_df['smoothed_unmeth_events_per_pos'])

    def get_masked_mbias_df(self, mode, cutting_sites: 'MbiasCuttingSites'):
        """Get M-bias counts in dataframe format, with positions in M-bias trimming zones set to NA
        """
        if mode == 'minimal':
            if self._minimal_mask_mbias_stats_df is not None:
                return self._minimal_mask_mbias_stats_df
            else:
                cutting_sites_df = cutting_sites.get_minimal_cutting_sites_as_df()
        elif mode == 'adjusted':
            if self._adjusted_mask_mbias_stats_df is not None:
                return self._adjusted_mask_mbias_stats_df
            else:
                cutting_sites_df = cutting_sites.get_adjusted_cutting_sites_as_df(mbias_data=self,
                                                                                  calling_mode='adjusted')
        else:
            raise NotImplementedError('Trimming mode {} unknown'.format(mode))

        def mask_read_positions_in_trimming_zone(df):
            # Hierarchical index with three levels: bsseq_strand, flen, pos; get strand and flen in current group df
            labels = df.index.unique()[0][0:2]
            start, end = cutting_sites_df.loc[labels, ['start', 'end']]
            idx = pd.IndexSlice
            df.loc[idx[:, :, 1:start], :] = np.nan
            df.loc[idx[:, :, end:], :] = np.nan
            return df

        masked_df = (self.mbias_stats_df
                     .groupby(level=['bsseq_strand', 'flen'], axis='index')
                     .apply(mask_read_positions_in_trimming_zone))

        if mode == 'minimal':
            self._minimal_mask_mbias_stats_df = masked_df
        else:
            self._adjusted_mask_mbias_stats_df = masked_df

        return masked_df



class MbiasDataPlotter:
    def __init__(self, max_flen_considered_for_trimming, mbias_data: MbiasData):
        self.mbias_data = mbias_data
        self.max_flen = max_flen_considered_for_trimming

    def flen_strat_plot(self, output_path, mask_mode=None, cutting_sites=None):
        if mask_mode:
            if not cutting_sites:
                print("Can't use mask dataframe, missing cutting sites information")
                mbias_stats_df = self.mbias_data.mbias_stats_df
            else:
                mbias_stats_df = self.mbias_data.get_masked_mbias_df(mode=mask_mode, cutting_sites=cutting_sites)
        else:
            mbias_stats_df = self.mbias_data.mbias_stats_df

        # TODO: parameter for number of fragments (or: all fragments?)
        plotting_data = (mbias_stats_df
                         .loc[pd.IndexSlice[:, range(0, self.max_flen + 1, 30), :], 'smoothed_beta_values']
                         .reset_index())

        plotting_data = plotting_data.dropna(axis='index', how='any')
        g = (sns.FacetGrid(data=plotting_data, col='bsseq_strand',
                           col_order='w_bc w_bc_rv c_bc c_bc_rv'.split(),
                           col_wrap=2,
                           hue='flen')
             .map(plt.plot, 'pos', 'smoothed_beta_values')
             .add_legend())
        g.fig.savefig(output_path)


class MbiasCuttingSites:
    """Calculate and provide M-bias cutting sites in dataframe and array format"""

    def __init__(self):
        self.cutting_site_calling_mode_dict = {
            'standard': self.max_plateau_length_threshold_std
        }
        self.minimal_cutting_sites = None
        self._adjusted_cutting_sites = None

    def get_adjusted_cutting_sites_as_array(self, mbias_data: MbiasData, calling_mode, **kwargs):
        """cutting site df structure
                              start    end    optional_fields
        bsseq_strand  flen
        'c_bc'        1       0        0

        """
        raise NotImplementedError

    def get_adjusted_cutting_sites_as_df(self, mbias_data: MbiasData, calling_mode, **kwargs):
        """ Calling function takes mbias count data for one stratum (bsseq_strand, flen)
        and returns Series with start and end of plateau, as well as arbitrary other fields
        """
        if self._adjusted_cutting_sites is None:
            calling_function = self.cutting_site_calling_mode_dict[calling_mode]
            cutting_site_df = (mbias_data.mbias_stats_df.groupby(axis='index', level=['bsseq_strand', 'flen'])
                               .apply(calling_function, **kwargs))
            self._adjusted_cutting_sites = cutting_site_df
            return cutting_site_df
        else:
            return self._adjusted_cutting_sites

    @staticmethod
    def max_plateau_length_threshold_std(group_df, min_len, max_std):
        def get_plateau_start_and_end(beta_values: pd.Series):
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


    def get_minimal_cutting_sites_as_array(self, relative_cutting_site_dict,
                                           max_flen_considered_for_trimming, max_read_length_bp):
        """
        Relative cutting site dict:
        {
            'w_bc': [0, 9]
            'c_bc': [0, 9]
            'w_bc_rv': [9, 0]
            'c_bc_rv': [9, 0]
        }
        """

        """ Create res array

        bsseq_strand      start_or_end       flen
        'w_bc'            0                  1
                                             2
                                             ...
                                             max_flen_considered_for_trimming
                          1                  1
                                             ...
        'w_bc_rv'         ...                ...

        """
        if not self.minimal_cutting_sites:
            res = np.zeros([4, 2, max_flen_considered_for_trimming + 1], dtype=np.int32)

            for bsseq_strand_name, bsseq_strand_index in b_inds._asdict().items():
                # Set start position
                res[bsseq_strand_index, 0, :] = relative_cutting_site_dict[bsseq_strand_name][0]

                # Set end position
                for curr_flen in range(max_flen_considered_for_trimming + 1):
                    max_allowed_pos_in_fragment = curr_flen - relative_cutting_site_dict[bsseq_strand_name][1]
                    max_allow_pos_in_read = (max_allowed_pos_in_fragment
                                             if max_allowed_pos_in_fragment <= max_read_length_bp
                                             else max_read_length_bp)
                    res[bsseq_strand_index, 1, curr_flen] = max_allow_pos_in_read

            self.minimal_cutting_sites = res
            return res

    def get_minimal_cutting_sites_as_df(self):
        raise NotImplementedError


def cutting_sites_area_plot(mbias_cutting_sites_df: pd.DataFrame, output_path):
    plot_df = mbias_cutting_sites_df.reset_index()
    g = sns.FacetGrid(data=plot_df, col='bsseq_strand', col_wrap=2, col_order='w_bc w_bc_rv c_bc c_bc_rv'.split())
    g = g.map(plt.fill_between, 'flen', 'start', 'end')
    g.fig.savefig(output_path)
