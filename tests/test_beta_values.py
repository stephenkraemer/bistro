from unittest import mock
import numpy as np
import math
import pytest
import pandas as pd
import os

import mqc

IDX = pd.IndexSlice

class TestBetaValueCounter:
    def test_event_counting(self, pileup_motifs_list, config):
        beta_value_counter_minimal = mqc.beta_values.StratifiedBetaValueCounter(config)
        for motif_pileups, curr_idx_pos in pileup_motifs_list:
            n_meth, n_unmeth, beta = beta_value_counter_minimal.update_and_return_total_beta_value(
                motif_pileups, curr_idx_pos)
        assert beta_value_counter_minimal.beta_counter.any


@pytest.fixture
def beta_value_counters():
    beta_counter_array_1 = np.zeros([7, 1001])
    beta_counter_array_1[:, -1] = 100
    beta_counter_array_1[:, 0] = 50
    beta_counter_array_2 = beta_counter_array_1 * 2

    beta_value_counter_1 = mock.Mock(beta_counter=beta_counter_array_1)
    beta_value_counter_2 = mock.Mock(beta_counter=beta_counter_array_2)
    return beta_value_counter_1, beta_value_counter_2

@pytest.fixture()
def beta_value_data_df():
    idx = pd.MultiIndex.from_product(
            iterables=[ ['global', 'cpg_islands'],
                        ['minimal', 'adjusted'],
                        'c_bc c_bc_rv w_bc w_bc_rv mate1 mate2 total'.split(),
                        np.arange(0, 1.1, 0.1) ],
            names=['region', 'trimming_status', 'bsseq_strand', 'beta_value'])
    counts = [10, 8, 6, 0, 0, 5, 0, 0, 8, 12, 20]
    # sum(counts) -> 69
    counts_broadcast = counts * 2 * 2 * 7
    df = pd.DataFrame(counts_broadcast, index=idx, columns=['count'])
    df = df.sort_index()
    return df

class TestBetaValueData:
    def test_basic_function(self, config):
        beta_value_data = mqc.beta_values.BetaValueData(config)
        beta_counter_array = np.zeros([7, 1001])
        beta_counter_array[:, -1] = 100
        beta_counter_array[:, 0] = 50
        beta_value_counter = mock.Mock(beta_counter=beta_counter_array)
        beta_value_data.add_data_from_counter(beta_value_counter,
                                              region_str='global',
                                              trimming_status_str='minimal')
        beta_value_count_at_no_meth = beta_value_data.df.loc[('global', 'minimal', 'w_bc', 0)]
        # Returns a series, must take scalar for boolean test, check size first
        assert beta_value_count_at_no_meth.size == 1
        assert beta_value_count_at_no_meth[0] == 50
        assert beta_value_data.df.loc[('global', 'minimal', 'w_bc', 1)][0] == 100

    def test_adding_up_counters(self, beta_value_counters, config):
        """Test whether two beta value counters are concatenated correctly

        These will come from the minimal and adjusted beta value calling
        """

        beta_value_data = mqc.beta_values.BetaValueData(config)
        beta_value_counter_1, beta_value_counter_2 = beta_value_counters

        beta_value_data.add_data_from_counter(beta_value_counter_1,
                                              region_str='global',
                                              trimming_status_str='minimal')
        beta_value_data.add_data_from_counter(beta_value_counter_2,
                                              region_str='global',
                                              trimming_status_str='adjusted')


        beta_value_count_1_at_no_meth = beta_value_data.df.loc[('global', 'minimal', 'w_bc', 0)]
        # Returns a series, must take scalar for boolean test, check size first
        assert beta_value_count_1_at_no_meth.size == 1
        assert beta_value_count_1_at_no_meth[0] == 50
        assert beta_value_data.df.loc[('global', 'minimal', 'w_bc', 1)][0] == 100
        assert beta_value_data.df.loc[('global', 'adjusted', 'w_bc', 1)][0] == 200
        assert beta_value_data.df.loc[('global', 'adjusted', 'w_bc', 0)][0] == 100

    def test_compute_frequencies(self, beta_value_data_df, config):
        beta_data = mqc.beta_values.BetaValueData(config)
        beta_data.df = beta_value_data_df
        beta_data._compute_frequencies()
        row_idx1 =  ('global', 'minimal', 'w_bc', 0)
        assert math.isclose(beta_data.df.loc[row_idx1, 'Frequency'], 10/69)
        row_idx2 = ('cpg_islands', 'adjusted', 'c_bc', 1.0)
        assert math.isclose(beta_data.df.loc[row_idx2, 'Frequency'], 20/69)

    def test_get_smoothed_frequencies(self, beta_value_data_df, config):
        # Incremental test: will only work if ._compute_frequencies works
        # Just to test that this runs through,
        # no direct assertion of smoothing result, check in test plots
        beta_data = mqc.beta_values.BetaValueData(config)
        beta_data.df = beta_value_data_df
        beta_data.add_smoothed_beta_value_frequencies()

        # check again that frequency computation worked
        row_idx1 =  ('global', 'minimal', 'w_bc', 0)
        assert math.isclose(beta_data.df.loc[row_idx1, 'Frequency'], 10/69)

        assert (beta_data.df.loc[row_idx1, 'Smoothed_frequency'] < 10/69)
        assert isinstance(
                beta_data.df.loc[row_idx1, 'Smoothed_frequency'],
                np.float)




class TestBetaValuePlotter:
    def test_base_plot(self, beta_value_data_df, test_output_dir, config):
        beta_value_data = mqc.beta_values.BetaValueData(config)
        beta_value_data.df = beta_value_data_df
        beta_value_data.add_smoothed_beta_value_frequencies()

        beta_value_plotter = mqc.beta_values.BetaValuePlotter(
                beta_value_data, config)
        out_basepath = os.path.join(test_output_dir, 'test_beta_value_dist')
        beta_value_plotter.beta_value_dist_plot(out_basepath_abs=out_basepath)
