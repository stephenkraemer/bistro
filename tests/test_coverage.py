"""Tests for coverage calculation from stream and coverage plots"""

import pytest
import os
import mqc
import pandas as pd
import numpy as np
from unittest.mock import Mock

@pytest.fixture(scope='module')
def config_cov_test(config):
    config['coverage_analysis']['max_per_cpg_cov'] = 3
    return config

@pytest.fixture()
def coverage_data_df(config_cov_test):
    max_cov = config_cov_test['coverage_analysis']['max_per_cpg_cov']
    df = pd.DataFrame({
        'Sequence_motif': 'CG',
        'Trimming_status': 'adjusted',
        'CpG_coverage': range(0, max_cov + 1),
        'Counts': [0,1,2,3],
        'Relative_frequency': np.nan,
    })
    df2 = df.copy()
    df2.loc[:, 'Trimming_status'] = 'minimal'
    df2.loc[:, 'Counts'] = [0,2,4,6]
    df_tot = (pd.concat([df, df2])
              .set_index(['Sequence_motif', 'Trimming_status', 'CpG_coverage']))
    return df_tot


class TestCoverageCounter:
    def test_update(self, config_cov_test):
        coverage_counter = mqc.coverage.CoverageCounter(config=config_cov_test,
                                                        trimming_status='minimal')
        coverages = [0, 0,
                     1, 1, 1,
                     2, 2, 2,
                     3, 3, 3, 3,
                     4, 4, 4]
        for curr_cov in coverages:
            coverage_counter.update(curr_cov)
        assert coverage_counter.per_cpg_cov_counts[0] == 2
        assert coverage_counter.per_cpg_cov_counts[2] == 3
        assert coverage_counter.per_cpg_cov_counts[3] == 7

class TestCoverageData:
    def test_add_counter(self, config_cov_test):
        coverage_data = mqc.coverage.CoverageData(config=config_cov_test)
        counter_min = Mock(per_cpg_cov_counts = [0, 1, 2, 3],
                       trimming_status='minimal')
        counter_adjusted = Mock(per_cpg_cov_counts = [0, 1, 2, 3],
                                trimming_status='adjusted')
        coverage_data.add_counter(counter_min)
        coverage_data.add_counter(counter_adjusted)
        print(coverage_data.df)
        cov_1_2_counts_min = (coverage_data.df
                          .query('Trimming_status == "minimal" and '
                                 'CpG_coverage in [1,2] and '
                                 'Sequence_motif == "CpG"')
                          .loc[:, 'Counts'])
        assert (cov_1_2_counts_min == np.array([1,2])).all()

        cov_0_3_counts_adj = (coverage_data.df
                              .query('Trimming_status == "adjusted" and '
                                     'CpG_coverage in [0,3] and '
                                     'Sequence_motif == "CpG"')
                              .loc[:, 'Counts'])
        assert (cov_0_3_counts_adj == np.array([0,3])).all()

    def test_calculate_relative_frequencies(self, config_cov_test, coverage_data_df):
        coverage_data = mqc.coverage.CoverageData(config=config_cov_test)
        coverage_data.df = coverage_data_df
        print(coverage_data.df)
        coverage_data.calculate_relative_frequencies()
        assert (coverage_data.df
                .query('Trimming_status == "adjusted" and CpG_coverage in [0, 3]')
                .loc[:, 'Relative_frequency']== [0, 0.5]).all()


class TestCoveragePlotter:
    def test_cov_histogram(self, config_cov_test, coverage_data_df, test_output_dir):
        coverage_data = mqc.coverage.CoverageData(config=config_cov_test)
        coverage_data.df = coverage_data_df
        coverage_data.calculate_relative_frequencies()
        coverage_plotter = mqc.coverage.CoveragePlotter(coverage_data, config=config_cov_test)
        output_path = os.path.join(test_output_dir, 'cov_frequency_hist.png')
        coverage_plotter.plot_cov_histogram(output_path)
        output_path = os.path.join(test_output_dir, 'cov_counts_hist.png')
        coverage_plotter.plot_cov_histogram(output_path, show_frequency=False)
