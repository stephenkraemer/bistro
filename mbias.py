import pickle
import pandas as pd
import numpy as np
import pytoml

import matplotlib
import mqc

matplotlib.use('Agg')  # import before pyplot import!
import matplotlib.pyplot as plt
import seaborn as sns

b_flags = mqc.flag_and_index_values.bsseq_strand_flags
b_inds = mqc.flag_and_index_values.bsseq_strand_indices
m_flags = mqc.flag_and_index_values.methylation_status_flags


def main():
    with open('./config.default.toml') as f_toml:
        config_dict = pytoml.load(f_toml)

    mbias_stats_array_p = '/home/kraemers/projects/mqc/mqc/test/results/mbias_stats_array.p'
    with open(mbias_stats_array_p, 'rb') as f:
        mbias_stats_arr = pickle.load(f)

    # relative_cutting_site_dict=config_dict['trimming']['relative_to_fragment_ends'],
    mbias_data = MbiasData(
        max_flen_considered_for_trimming=config_dict['trimming']['max_flen_considered_for_trimming'],
        max_read_length_bp=config_dict['data_properties']['max_read_length_bp'],
        required_n_events_for_cutting_site_determination=(
            config_dict['trimming']['required_n_events_for_cutting_site_determination']),
        max_window_size_for_smoothing=config_dict['trimming']['max_window_size_for_smoothing'],
        mbias_stats_array=mbias_stats_arr)

    plotter = MbiasPlotter(
        max_flen_considered_for_trimming=config_dict['trimming']['max_flen_considered_for_trimming'],
        mbias_df=mbias_data.mbias_stats_df
    )

    plotter.flen_strat_plot(output_path='/home/kraemers/temp/mbias_flen_strat.tmp.png')


class MbiasData:
    """Basic formatting and processing of M-bias counts into dataframe structure"""

    def __init__(self, max_flen_considered_for_trimming, max_read_length_bp, mbias_stats_array,
                 required_n_events_for_cutting_site_determination, max_window_size_for_smoothing):
        self.max_window_size_for_smoothing = max_window_size_for_smoothing
        self.required_n_events_for_cutting_site_determination = required_n_events_for_cutting_site_determination
        self.mbias_stats_array = mbias_stats_array
        self.max_read_length_bp = max_read_length_bp
        self.max_flen_considered_for_trimming = max_flen_considered_for_trimming
        self.mbias_stats_df = pd.DataFrame()
        self.mbias_stats_array_to_df()

    def mbias_stats_array_to_df(self):
        self.convert_mbias_arr_info_to_df_format()
        self.add_smoothed_mbias_stats()
        self.add_beta_values()

    def convert_mbias_arr_info_to_df_format(self):
        rows = []
        for bsseq_strand_name, bsseq_strand_ind in b_inds._asdict().items():
            for flen in range(1, self.max_flen_considered_for_trimming + 1):
                for pos in range(1, self.max_read_length_bp):
                    row = dict()
                    row['bsseq_strand'] = bsseq_strand_name
                    row['flen'] = flen
                    row['pos'] = pos
                    row['meth_events_per_pos'] = self.mbias_stats_array[bsseq_strand_ind, flen, pos]
                    row['unmeth_events_per_pos'] = self.mbias_stats_array[bsseq_strand_ind + 1, flen, pos]
                    rows.append(row)

        df = pd.DataFrame(rows).set_index(['bsseq_strand', 'flen', 'pos'])
        df.loc[:, 'smoothed_meth_events_per_pos'] = np.nan
        df.loc[:, 'smoothed_unmeth_events_per_pos'] = np.nan
        self.mbias_stats_df = df.astype(dtype=np.float64, copy=True)

    def add_smoothed_mbias_stats(self):
        df = self.mbias_stats_df
        win_cols = ['smoothed_meth_events_per_pos', 'smoothed_unmeth_events_per_pos']
        flen_data_cols = ['meth_events_per_pos', 'unmeth_events_per_pos']
        for bsseq_strand in 'C_BC C_BC_RV W_BC W_BC_RV'.split():
            minimal_flen_with_enough_coverage_reached = False
            # TODO: allow all fragment lengths (currently computation takes a while,
            #       so I am not doing this while testing)
            # for flen in range(max_flen_considered_for_trimming, 0, -1):
            for flen in range(200, 40, -1):
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


class MbiasPlotter:
    def __init__(self, max_flen_considered_for_trimming, mbias_df: pd.DataFrame):
        self.df = mbias_df
        self.max_flen = max_flen_considered_for_trimming

    def flen_strat_plot(self, output_path):
        plotting_data = (self.df
                         .loc[pd.IndexSlice[:, range(1, self.max_flen + 1, 30), :], 'smoothed_beta_values']
                         .reset_index())

        plotting_data = plotting_data.dropna(axis='index', how='any')
        g = (sns.FacetGrid(data=plotting_data, col='bsseq_strand',
                           col_order='C_BC C_BC_RV W_BC W_BC_RV'.split(),
                           col_wrap=2,
                           hue='flen')
             .map(plt.plot, 'pos', 'smoothed_beta_values')
             .add_legend())
        g.fig.savefig(output_path)


class MbiasCounter:
    def __init__(self, max_read_length, min_phred_score, max_flen_considered_for_trimming):
        self.counter = np.zeros([8, max_flen_considered_for_trimming + 1, max_read_length + 1], dtype='i4')
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
                        row_index_selector = 0
                    elif meth_status_flag == 4:
                        row_index_selector = 1
                    else:  # SNP, Ref base
                        continue

                    if pileup_read.bs_seq_strand_flag == b_flags.c_bc:
                        strand_and_meth_status_based_index = b_inds.c_bc
                    elif pileup_read.bs_seq_strand_flag == b_flags.c_bc_rv:
                        strand_and_meth_status_based_index = b_inds.c_bc_rv
                    elif pileup_read.bs_seq_strand_flag == b_flags.w_bc:
                        strand_and_meth_status_based_index = b_inds.w_bc
                    elif pileup_read.bs_seq_strand_flag == b_flags.w_bc_rv:
                        strand_and_meth_status_based_index = b_inds.w_bc_rv
                    else:
                        continue

                    strand_and_meth_status_based_index += row_index_selector
                    self.counter[strand_and_meth_status_based_index][tlen][pos_in_read] += 1


if __name__ == '__main__':
    main()
