import numpy as np
import mqc
import matplotlib
matplotlib.use('Agg')  # import before pyplot import!
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('darkgrid')


class CoverageCounter:
    def __init__(self, config):
        self.max_cov = config["coverage_analysis"]["max_per_cpg_cov"]
        self.per_cpg_cov_counts = np.zeros(self.max_cov + 1)
        self.per_cpg_cov_frequency = np.array([])

    def update(self, per_cpg_cov):
        if per_cpg_cov > self.max_cov:
            self.per_cpg_cov_counts[self.max_cov] += 1
        else:
            self.per_cpg_cov_counts[per_cpg_cov] += 1

    def calculate_relative_frequencies(self):
        self.per_cpg_cov_frequency = self.per_cpg_cov_counts / self.per_cpg_cov_counts.sum()


class CoveragePlotter:
    def __init__(self, coverage_counter: 'mqc.coverage.CoverageCounter'):
        self.max_cov = coverage_counter.max_cov
        self.coverage_counts = coverage_counter.per_cpg_cov_counts
        self.coverage_freqs = coverage_counter.per_cpg_cov_frequency

    def plot_cov_histogram(self, output_path, show_frequency=True):
        fig, ax = plt.subplots(1, 1)

        if show_frequency:
            y_values = self.coverage_freqs
            y_label = 'Relative frequency'
        else:
            y_values = self.coverage_counts
            y_label = 'Count'

        ax: plt.Axes
        x_values = np.arange(0, self.max_cov+1)
        ax.set_xlabel('Per CpG coverage')
        ax.set_ylabel(y_label)

        sns.barplot(x = x_values, y = y_values, ax=ax, color='grey', edgecolor='grey')
        fig.savefig(output_path)
