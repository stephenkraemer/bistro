"""Sequence motif (CG, CH, CHH, CHG) coverage calculation from stream and plots"""

import numpy as np
import pandas as pd
from typing import List
import matplotlib
matplotlib.use('Agg')  # import before pyplot import!
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
sns.set_style('darkgrid')


class CoverageCounter:
    def __init__(self, config, trimming_status):
        self.trimming_status = trimming_status
        self.max_cov = config["coverage_analysis"]["max_per_cpg_cov"]
        self.per_cpg_cov_counts = np.zeros(self.max_cov + 1, dtype=np.int32)

    def update(self, per_cpg_cov):
        if per_cpg_cov > self.max_cov:
            self.per_cpg_cov_counts[self.max_cov] += 1
        else:
            self.per_cpg_cov_counts[per_cpg_cov] += 1


class CoverageData:
    def __init__(self, config):
        """

        Parameters
        ----------
        config

        Attributes
        ----------
        self.df: pd.DataFrame
                                                                      'Counts'    'Relative_frequency'
        'Sequence_motif'  'Trimming_status'    'CpG_coverage'
        'CG'              'adjusted'           0                      100         0.1
        'CG'              'adjusted'           1                      200         0.2
        """

        self.coverage_counters: List[CoverageCounter] = []
        self.df = pd.DataFrame()
        self.max_cov = config["coverage_analysis"]["max_per_cpg_cov"]
        self.config = config

    def add_counter(self, counter: CoverageCounter):
        self.coverage_counters.append(counter)
        self._add_arr_to_df(counter)

    def _add_arr_to_df(self, counter: CoverageCounter):
        # TODO-learn: columns='Counts'
        curr_df = pd.DataFrame({
            'Sequence_motif': 'CpG',
            'Trimming_status': counter.trimming_status,
            'CpG_coverage': range(0, self.max_cov + 1),
            'Counts': counter.per_cpg_cov_counts,
            'Relative_frequency': np.nan,
        })
        curr_df.set_index(['Sequence_motif', 'Trimming_status', 'CpG_coverage'], inplace=True)
        # TODO-learn: append takes index into account! also with empty frame?
        self.df = self.df.append(curr_df)


    def calculate_relative_frequencies(self) -> None:

        def get_relative_frequency(ser):
            return ser / ser.sum()

        self.df['Relative_frequency'] = (self.df
                                         .groupby(level=['Sequence_motif', 'Trimming_status'])['Counts']
                                         .transform(get_relative_frequency))

    # def save(self, out_paths):
    #     for k, v in self.config['sample'].items():
    #         out_df_counts = self.df.copy()
    #         out_df_counts[k] = v
    #     with open




class CoveragePlotter:
    def __init__(self, coverage_data: 'CoverageData', config):
        self.coverage_data = coverage_data
        self.max_cov = config["coverage_analysis"]["max_per_cpg_cov"]
        self.sample_name = config['sample']['name']

    def plot_cov_histogram(self, output_path, show_frequency=True):
        if show_frequency:
            y_column = 'Relative_frequency'
            y_label = 'Relative frequency'
        else:
            y_column = 'Counts'
            y_label = 'Counts'

        g = (sns.factorplot(x='CpG_coverage', y=y_column, data=self.coverage_data.df.reset_index(),
                            row='Sequence_motif', hue='Trimming_status', kind='bar',
                            sharex=True, sharey=True,
                            size=cm_to_inch(8), aspect=2, legend=False, legend_out=False)
            .set_axis_labels('Coverage', y_label)
            .add_legend(title='Trimming mode')
            .set_titles(row_template='{row_name}')
             )

        for curr_ax in g.axes.flat:
            curr_ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=30, integer=True))

        # TODO-snippet
        g.fig.subplots_adjust(top=0.85)
        g.fig.suptitle(self.sample_name, x=0.53, ha='center')

        g.fig.savefig(output_path)



def cm_to_inch(cm):
    return cm / 2.54
