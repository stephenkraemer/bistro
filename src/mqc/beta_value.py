from mqc.visitors import Counter
from mqc.pileup.pileup import MotifPileup
import numpy as np
from typing import Dict
import pandas as pd
import matplotlib
matplotlib.use('Agg') # import before pyplot import!
import matplotlib.pyplot as plt
import seaborn as sns
import os.path as op
import os
import mqc.flag_and_index_values as mfl


bstrands = mfl.bsseq_strand_indices
bstrat = mfl.beta_value_stratum_indices
roi_strats = mfl.region_of_interest_indices


class BedFile:
    def __init__(self, bedPath: str, chrom: str, idx: int):
        with open(bedPath, mode='r') as bedfile:
            self.idx = idx
            self.name = op.basename(op.splitext(bedPath)[0])
            self.intervals = [tuple(map(int, line.rstrip('\n').split('\t')[1:3]))
                              for line in bedfile.readlines()
                              if line.split('\t')[0].lstrip('chr') == chrom].__iter__()
            self.curr_interval = (-1, -1)
    def __iter__(self):
        return self
    def __next__(self):
        self.curr_interval = next(self.intervals)
        return self.curr_interval


class StratifiedBetaCounter(Counter):
    def __init__(self, config: Dict, chrom: str=''):


        save_stem = config['paths']['stratified_beta_counts']
        motifs = config['run']['motifs']

        self.motif_idx_dict = {motif: i for i, motif in enumerate(motifs)}
        self.roi_list = []
        self.bins = 1000

        if config['run']['roi_index_files'] and chrom:
            for idx, bedfile in enumerate(config['run']['roi_index_files'].split(',')):
                #Index of ROI namedtuples are 1-based, idx 0 is whole genome
                self.roi_list.append(BedFile(bedfile, chrom, idx+1))
        elif config['run']['roi_index_files'] and not chrom:
            raise RuntimeError('ROI bedfile has been specified but no chromosome')


        strat_dict = {bstrands.c_bc: 'c_bc',
                      bstrands.c_bc_rv: 'c_bc_rv',
                      bstrands.w_bc: 'w_bc',
                      bstrands.w_bc_rv: 'w_bc_rv',
                      bstrat.mate1: 'mate_1',
                      bstrat.mate2: 'mate_2',
                      bstrat.all: 'all'}

        strat_list_sorted = [strat_dict[x] for x in sorted(strat_dict.keys())]

        self.sample_meta = config['sample'].copy()

        dim_names = ['motif', 'roi_type', 'bs_strand', 'beta_value']

        dim_levels = [motifs,
                      ['whole_genome']+[roi.name for roi in self.roi_list],
                      strat_list_sorted,
                      [np.round(x / self.bins, decimals=3) for x in range(0, self.bins + 1)]]

        array_shape = [len(motifs),
                       1+len(self.roi_list),
                       7,
                       self.bins+1]

        counter_array = np.zeros(array_shape, dtype='u8')

        super().__init__(dim_names=dim_names,
                         dim_levels=dim_levels,
                         counter_array=counter_array,
                         save_stem=save_stem)

    def process(self, motif_pileup: MotifPileup):

        #Create a list containing the indices of each bedfile covering the MPU
        roi_hit_list = [roi_strats.whole_genome]

        start = motif_pileup.idx_pos.start
        end = motif_pileup.idx_pos.end

        for bedfile in self.roi_list:
            if start >= bedfile.curr_interval[1]:
                for interval in bedfile:
                    # MPU beyond current ROI tuple
                    if start >= interval[1]:
                        # bedfile.next_interval()
                        continue
                    # MPU before current ROI tuple
                    elif end <= interval[0]:
                        break
                    # MPU inside current ROI tuple
                    else:
                        roi_hit_list.append(bedfile.idx)
                        break

            # MPU before current ROI tuple
            elif end <= bedfile.curr_interval[0]:
                continue
            # MPU inside current ROI tuple
            else:
                roi_hit_list.append(bedfile.idx)

        for strat_idx, beta in enumerate(motif_pileup.strat_beta_arr):
            if np.isfinite(beta):
                self.counter_array[self.motif_idx_dict[motif_pileup.idx_pos.motif],
                                   roi_hit_list,
                                   strat_idx,
                                   int(beta * self.bins)] += 1

    def get_dataframe(self):
        if self._counter_dataframe.empty:
            counter_dataframe = super().get_dataframe().reset_index()
            for key in self.sample_meta.keys():
                counter_dataframe.insert(0, key, self.sample_meta[key])
            self._counter_dataframe = counter_dataframe.set_index(list(self.sample_meta.keys()) + self.dim_names)

        return self._counter_dataframe




def stratified_beta_hist(strat_beta_df, config):

    strat_df = strat_beta_df.reset_index()

    strata_categories = {'strand-wise': ['c_bc', 'c_bc_rv', 'w_bc', 'w_bc_rv'],
                         'mate-wise': ['mate_1', 'mate_2'],
                         'total': ['all']}

    for roi_type in strat_df['roi_type'].unique():

        for strata_cat, strata in strata_categories.items():
            plot_df = strat_df.copy().loc[(strat_df.roi_type == roi_type) & (strat_df.bs_strand.isin(strata))]
            plot_df['bs_strand'] = plot_df['bs_strand'].str.upper().str.replace('_', '-')
            #plot_df['roi_type'] = plot_df['roi_type'].str.replace('_', ' ')

            g = (sns.FacetGrid(plot_df, col='bs_strand', row='motif', sharey='row')
                 .map(plt.plot, 'beta_value', 'counts')
                 .set_titles("{row_name} | {col_name}")
                 .set_axis_labels('Beta value', 'Frequency'))
            g.fig.tight_layout()
            g.fig.savefig(f"{config['paths']['stratified_beta_hist_trunk']}_{roi_type}_{strata_cat}.png")


def analyze_stratified_beta_values(config):
    os.makedirs(config['paths']['qc_stats_dir'], exist_ok=True, mode=0o770)
    strat_beta_df = pd.read_pickle(config['paths']['stratified_beta_counts']+'.p')
    stratified_beta_hist(strat_beta_df, config)
    plt.close('all')
