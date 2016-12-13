def view_df(df):
    import os
    import subprocess
    import time
    timestamp = time.strftime('%H_%M_%S')
    target_path = os.path.expanduser(os.path.join('~/temp/', timestamp + '_test.csv'))
    df.to_csv(target_path)
    subprocess.Popen(['localc', target_path])


import pickle
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Indices
C_BC_IND = 0
C_BC_RV_IND = 2
W_BC_IND = 4
W_BC_RV_IND = 6

max_read_length_bp = 101
max_flen_considered_for_trimming = 500
required_N_events = 3000

mbias_stats_array_p = '/home/kraemers/projects/mqc/mqc/test/results/mbias_stats_array.p'


def main():
    with open(mbias_stats_array_p, 'rb') as f:
        mbias_arr = pickle.load(f)
    df = mbias_arr_to_df(mbias_arr)
    df = call_flen_mbias_stats_in_windows(df)

    df['windowed_beta_values'] = df['windowed_meth_events_per_pos'] / (
        df['windowed_meth_events_per_pos'] + df['windowed_unmeth_events_per_pos'])

    plotting_data = df.loc[
        pd.IndexSlice[:, range(1, max_flen_considered_for_trimming + 1, 10), :], 'windowed_beta_values'].reset_index()

    plotting_data = plotting_data.dropna(axis='index', how='any')
    g = sns.FacetGrid(data=plotting_data, col='bsseq_strand',
                      col_order='C_BC C_BC_RV W_BC W_BC_RV'.split(),
                      col_wrap=2,
                      hue='flen')
    g.map(plt.plot, 'pos', 'windowed_beta_values')
    g.fig.savefig('/home/kraemers/temp/test.png')


def mbias_arr_to_df(mbias_arr):
    rows = []
    for bsseq_strand, bsseq_strand_ind in zip('C_BC C_BC_RV W_BC W_BC_RV'.split(),
                                              [C_BC_IND, C_BC_RV_IND, W_BC_IND, W_BC_RV_IND]):
        for flen in range(1, max_flen_considered_for_trimming + 1):
            for pos in range(1, max_read_length_bp):
                row = dict()
                row['bsseq_strand'] = bsseq_strand
                row['flen'] = flen
                row['pos'] = pos
                row['meth_events_per_pos'] = mbias_arr[bsseq_strand_ind, flen, pos]
                row['unmeth_events_per_pos'] = mbias_arr[bsseq_strand_ind + 1, flen, pos]
                rows.append(row)

    df = pd.DataFrame(rows).set_index(['bsseq_strand', 'flen', 'pos'])
    df.loc[:, 'windowed_meth_events_per_pos'], df.loc[:, 'windowed_unmeth_events_per_pos'] = np.nan, np.nan
    return df


def call_flen_mbias_stats_in_windows(df):
    win_cols = ['windowed_meth_events_per_pos', 'windowed_unmeth_events_per_pos']
    flen_data_cols = ['meth_events_per_pos', 'unmeth_events_per_pos']
    for bsseq_strand in 'C_BC C_BC_RV W_BC W_BC_RV'.split():
        # for flen in range(max_flen_considered_for_trimming, 0, -1):
        for flen in range(180, 120, -1):
            curr_flen_rows = pd.IndexSlice[bsseq_strand, flen, :]
            df.loc[curr_flen_rows, win_cols] = df.loc[curr_flen_rows, flen_data_cols].values
            total_N_events = df.loc[curr_flen_rows, win_cols].sum().sum()
            if total_N_events < required_N_events:
                total_N_events_sufficient = False
                last_added_flen = flen
                while not total_N_events_sufficient:
                    last_added_flen -= 1
                    # TODO: at max go down to flen - K, if then not enough coverage, set to na?
                    if last_added_flen == 0:
                        df.loc[curr_flen_rows, win_cols] = np.nan
                        break
                    df.loc[curr_flen_rows, win_cols] += df.loc[
                        (bsseq_strand, last_added_flen, slice(None)), flen_data_cols].values
                    total_N_events = df.loc[curr_flen_rows, win_cols].sum().sum()
                    if total_N_events >= required_N_events:
                        total_N_events_sufficient = True
    return df


if __name__ == '__main__':
    # main()
    pass
