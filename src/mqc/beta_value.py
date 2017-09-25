from mqc.visitors import Counter
from mqc.pileup.pileup import MotifPileup
import numpy as np
from typing import Dict
import pandas as pd
import matplotlib
matplotlib.use('Agg') # import before pyplot import!
import matplotlib.pyplot as plt
import seaborn as sns
import os
import mqc.flag_and_index_values as mfl

b_inds = mfl.bsseq_strand_indices
bstrat = mfl.beta_value_stratum_indices


class StratifiedBetaCounter(Counter):
    def __init__(self, config: Dict):

        save_stem = config['paths']['stratified_beta_counts']
        dim_names = ['strat_idx', 'beta_value']

        dim_levels = [range(0, 7), range(0, 1001)]
        array_shape = [7, 1001]

        counter_array = np.zeros(array_shape, dtype='u8')

        super().__init__(dim_names=dim_names,
                         dim_levels=dim_levels,
                         counter_array=counter_array,
                         save_stem=save_stem)

    def process(self, motif_pileup: MotifPileup):

        for idx, beta in enumerate(motif_pileup.strat_beta_arr):
            if np.isfinite(beta):
                self.counter_array[idx, int(beta * 1000)] += 1


def stratified_beta_hist(strat_beta_df, config):

    strat_df = pd.DataFrame({'strat_idx': [b_inds.c_bc,
                                           b_inds.c_bc_rv,
                                           b_inds.w_bc,
                                           b_inds.w_bc_rv,
                                           bstrat.mate1,
                                           bstrat.mate2,
                                           bstrat.all],
                             'strat_label': ['C-BC',
                                             'C-BC-RV',
                                             'W-BC',
                                             'W-BC-RV',
                                             'Mate 1',
                                             'Mate 2',
                                             'All']})

    strat_beta_df = strat_beta_df.reset_index().merge(strat_df, on='strat_idx')
    strat_beta_df['beta_value'] = strat_beta_df['beta_value']/1000

    #TODO: Split this into three plots: stratified by strand, mate, all

    strata_categories = {'strand-wise': list(b_inds), 'mate-wise': [bstrat.mate1, bstrat.mate2], 'total': [bstrat.all]}

    for strata_cat, strata in strata_categories.items():
        col_wrap = 2 if len(strata) > 1 else 1
        plot_df = strat_beta_df.loc[strat_beta_df['strat_idx'].isin(strata)]
        g = (sns.FacetGrid(plot_df, col='strat_label', col_wrap=col_wrap)
             .map(plt.plot, 'beta_value', 'counts')
             .set_titles("{col_name}")
             .set_axis_labels('Beta value', 'Frequency')
            )
        g.fig.tight_layout()
        g.fig.savefig(f"{config['paths']['stratified_beta_hist_trunk']}_{strata_cat}.png")


def analyze_stratified_beta_values(config):
    os.makedirs(config['paths']['qc_stats_dir'], exist_ok=True, mode=0o770)
    strat_beta_df = pd.read_pickle(config['paths']['stratified_beta_counts']+'.p')
    stratified_beta_hist(strat_beta_df, config)
    plt.close('all')
