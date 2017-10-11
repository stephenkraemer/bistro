from mqc.visitors import Counter
from mqc.pileup.pileup import MotifPileup
from typing import Dict
import numpy as np

import pandas as pd
import matplotlib
matplotlib.use('Agg') # import before pyplot import!
import matplotlib.pyplot as plt
import seaborn as sns
import os


class CoverageCounter(Counter):
    """
    Count total coverage of motifs
    Coverage above the config-set maximum will be pooled in the last element of the counter_array
    """

    def __init__(self, config: Dict):

        motifs = config['run']['motifs']
        save_stem = config['paths']['cov_counts']
        self.max_coverage = config['coverage_analysis']['max_per_motif_cov']
        dim_names = ['motif', 'coverage']

        idx_motif_tuples = enumerate(config['run']['motifs'])
        self.motif_idx_dict = {motif: i for i, motif in idx_motif_tuples}

        dim_levels = [motifs, range(0, self.max_coverage+1)]
        array_shape = [len(motifs), self.max_coverage+1]

        counter_array = np.zeros(array_shape, dtype='u4')

        super().__init__(dim_names=dim_names,
                         dim_levels=dim_levels,
                         counter_array=counter_array,
                         save_stem=save_stem)

    def process(self, motif_pileup: MotifPileup):
        if (motif_pileup.n_total) > self.max_coverage:
            coverage = self.max_coverage
        else:
            coverage = motif_pileup.n_total
        event_class = (self.motif_idx_dict[motif_pileup.idx_pos.motif], coverage)
        self.counter_array[event_class] += 1


def coverage_hist(coverage_df, config):
    plot_df = coverage_df.reset_index()
    last_bin = plot_df.loc[plot_df['counts'] > 0].max().loc['coverage']
    g = (sns.FacetGrid(plot_df, col='motif')
         .map(plt.plot, 'coverage', 'counts')
         .set_titles("{col_name}")
         .set_axis_labels('Coverage', 'Frequency')
         .set(xlim=(last_bin/(-15), last_bin+2))
         )
    g.fig.tight_layout()
    g.fig.savefig(config['paths']['coverage_hist'])


def analyze_coverage(config):
    os.makedirs(config['paths']['qc_stats_dir'], exist_ok=True, mode=0o770)
    coverage_df = pd.read_pickle(config['paths']['cov_counts']+'.p')
    coverage_hist(coverage_df, config)
    plt.close('all')

