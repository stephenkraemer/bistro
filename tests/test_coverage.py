import pytest
import os
import mqc
from unittest.mock import Mock


class TestCoverageCounter:
    def test_update(self, config):
        coverage_counter = mqc.coverage.CoverageCounter(config)
        coverage_counter.max_cov = 3
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

class TestCoveragePlotter:
    @pytest.fixture()
    def coverage_counter(self):
        coverage_counter = Mock(per_cpg_cov_counts=[2,2,2,2],
                                per_cpg_cov_frequency=[0.2,0.2,0.2,0.2],
                                max_cov = 3)
        return coverage_counter
    def test_cov_histogram(self, coverage_counter, test_output_dir):
        coverage_plotter = mqc.coverage.CoveragePlotter(coverage_counter)
        output_path = os.path.join(test_output_dir, 'cov_frequency_hist.png')
        coverage_plotter.plot_cov_histogram(output_path)
        output_path = os.path.join(test_output_dir, 'cov_counts_hist.png')
        coverage_plotter.plot_cov_histogram(output_path, show_frequency=False)
