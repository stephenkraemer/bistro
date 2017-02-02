import textwrap
from unittest.mock import Mock
import io
import itertools
import mqc
import numpy as np
import os
import pandas as pd
import pytest
import time

import matplotlib
matplotlib.use('Agg')  # import before pyplot import!
import matplotlib.pyplot as plt


@pytest.fixture
def empty_mbias_data(config):
    max_flen_considered_for_trimming = config["trimming"]["max_flen_considered_for_trimming"]
    max_read_length = config["data_properties"]["max_read_length_bp"]
    mbias_counts = np.ones([4, max_flen_considered_for_trimming + 1, max_read_length + 1, 2], dtype='i4')
    mbias_counter = Mock(counter=mbias_counts)
    minimal_cutting_sites = Mock()
    minimal_cutting_sites.get_df.return_value = pd.DataFrame()
    mbias_data = mqc.mbias.MbiasData(mbias_counter=mbias_counter,
                                     config=config)
    return mbias_data


@pytest.fixture()
def mbias_stats_df():
    test_df_path = os.path.join(os.path.dirname(__file__), "mbias_data_frame.csv")
    df = pd.read_csv(test_df_path, sep='\t', header=0, index_col=[0, 1, 2],
                     nrows=None, usecols=None,
                     na_values=['#DIV/0!'], keep_default_na=True
                     )
    return df

# TODO-refactor: use this large dataframe in all tests
@pytest.fixture()
def mbias_stats_df_large_for_cutting_site_determin_tests(test_output_dir):
    """
    Flens 100 and 101 for all bsseq_strands
    Flen 100 is completely bad, while Flen 100 has a plateau between positions 21 and 80
    """

    bsseq_strands = 'c_bc c_bc_rv w_bc w_bc_rv'.split()
    flens = [100, 101]
    pos = np.arange(1, 102)


    data = itertools.product(bsseq_strands, flens, pos)
    df = pd.DataFrame(list(data), columns=['bsseq_strand', 'flen', 'pos'])
    # idx = pd.MultiIndex.from_product([bsseq_strands, flens, pos])

    bad_mbias_curve = np.random.normal(loc=0.3, scale=0.6, size=101)
    good_mbias_curve = np.random.normal(loc=0.7, scale=0.02, size=101)
    good_mbias_curve[0:20] = bad_mbias_curve[0:20]
    good_mbias_curve[80:] = bad_mbias_curve[80:]

    fig, axes = plt.subplots(1,2)
    axes[0].plot(bad_mbias_curve)
    axes[1].plot(good_mbias_curve)
    fig.savefig(os.path.join(test_output_dir, 'mbias_test_data.png'))

    df['beta_values_smoothed'] = 0
    for bsseq_strand in bsseq_strands:
        df.loc[(df.bsseq_strand == bsseq_strand) & (df.flen == 100), 'beta_values_smoothed'] = bad_mbias_curve
        df.loc[(df.bsseq_strand == bsseq_strand) & (df.flen == 101), 'beta_values_smoothed'] = good_mbias_curve

    df['beta_values'] = df['beta_values_smoothed']

    df = df.set_index(['bsseq_strand', 'flen', 'pos']).sort_index()
    return df

@pytest.fixture()
def mbias_stats_df_large_cutting_sites_df():
    # TODO-learn: textwrap.dedent
    # TODO-snippet: textwrap.dedent
    # TODO-document: cutting site df format
    df_file = io.StringIO(textwrap.dedent("""\
        bsseq_strand  flen  start  end  average_methylation
        'c_bc'        100   40     60   0.7
        'c_bc'        101   0      100  0.7
        'c_bc_rv'     100   40     60   0.7
        'c_bc_rv'     101   0      100  0.7
        'w_bc'        100   40     60   0.7
        'w_bc'        101   0      100  0.7
        'w_bc_rv'     100   40     60   0.7
        'w_bc_rv'     101   0      100  0.7"""))

    cutting_site_df = pd.read_csv(df_file, sep='\s+', header=0, index_col=[0,1])
    return cutting_site_df

@pytest.fixture()
def ya_mbias_stats_df():

    multi_idx = pd.MultiIndex.from_product(
            ['c_bc c_bc_rv w_bc w_bc_rv'.split(),
             range(1, 201),  # flen
             range(1, 101)   # pos
             ],
            names=['bsseq_strand', 'flen', 'pos'])

    n_meth_good = np.random.normal(70, 1, 100)
    n_unmeth_good = np.random.normal(30, 1, 100)

    df = pd.DataFrame(index=multi_idx)

    df['n_meth'] = np.tile(n_meth_good, 4 * 200)
    df['n_unmeth'] = np.tile(n_unmeth_good, 4 * 200)
    df['beta_value'] = (df['n_meth'] /
                        (df['n_unmeth'] + df['n_meth']))

    df['n_meth_smoothed'], df['n_unmeth_smoothed'] = df['n_meth'], df['n_unmeth']
    df['beta_value_smoothed'] = df['beta_value']

    return df



@pytest.fixture
def mbias_data_with_test_df(empty_mbias_data, mbias_stats_df):
    mbias_data = empty_mbias_data
    mbias_data.mbias_stats_df = mbias_stats_df

    # Adapt config parameters to test dataframe
    mbias_data.max_flen_considered_for_trimming = 5
    mbias_data.min_flen_considered_for_methylation_calling = 1
    mbias_data.required_n_events_for_cutting_site_determination = 12

    return mbias_data


@pytest.fixture()
def cutting_sites_df():
    # TODO-algorithm: are start and end positions zero-based, right-open?
    #                 I think it would good to see them as 1-based labels
    test_df = io.StringIO(textwrap.dedent("""\
        bsseq_strand  flen  start  end  average_methylation
        c_bc          1     2      3    0.7
        c_bc          2     2      3    0.7
        c_bc          3     2      3    0.7
        c_bc          4     2      3    0.7
        c_bc          5     2      3    0.7
        c_bc_rv       1     2      3    0.7
        c_bc_rv       2     2      3    0.7
        c_bc_rv       3     2      3    0.7
        c_bc_rv       4     2      3    0.7
        c_bc_rv       5     2      3    0.7
        w_bc          1     1      2    0.7
        w_bc          2     1      2    0.7
        w_bc          3     1      2    0.7
        w_bc          4     1      2    0.7
        w_bc          5     1      2    0.7
        w_bc_rv       1     1      2    0.7
        w_bc_rv       2     1      2    0.7
        w_bc_rv       3     1      2    0.7
        w_bc_rv       4     1      2    0.7
        w_bc_rv       5     1      2    0.7"""))

    cutting_sites_df = pd.read_csv(test_df, sep='\s+', header=0, index_col=[0, 1])
    return cutting_sites_df


class TestMbiasCounter:
    def test_update(self, pileup_motifs_list, config):
        mbias_counter = mqc.mbias.MbiasCounter(config)
        for motif_pileups, curr_idx_pos in pileup_motifs_list:
            mbias_counter.update(motif_pileups, curr_idx_pos)
        assert np.sum(mbias_counter.counter) > 30


@pytest.mark.incremental
class TestMbiasData:
    def test_counter_to_df_conversion(self, empty_mbias_data: 'mqc.mbias.MbiasData'):
        empty_mbias_data.convert_mbias_arr_info_to_df_format()
        assert empty_mbias_data.mbias_stats_df.loc[('c_bc', 1, 1), 'n_meth'] == 1

    def test_mbias_stats_smoothing(self, mbias_data_with_test_df: 'mqc.mbias.MbiasData'):
        mbias_data = mbias_data_with_test_df
        mbias_data.add_smoothed_mbias_stats()

        # Test that flens with enough coverage on their own are recognized
        assert mbias_data.mbias_stats_df.loc[('c_bc_rv', 1, 3), 'n_meth_smoothed'] == 4
        assert mbias_data.mbias_stats_df.loc[('w_bc_rv', 5, 1), 'n_meth_smoothed'] == 4

        # Test that flens with too low coverage are rescued by smoothing where possible
        assert mbias_data.mbias_stats_df.loc[('w_bc', 4, 1), 'n_meth_smoothed'] == 2

        # Test that flens with insufficient coverage and possibility for smoothing are discarded
        assert np.isnan(mbias_data.mbias_stats_df.loc[('w_bc', 2, 1), 'n_meth_smoothed'])


    def test_mask_df_computation(self, mbias_data_with_test_df: 'mqc.mbias.MbiasData', cutting_sites_df):
        mbias_data = mbias_data_with_test_df
        cutting_sites = Mock()
        cutting_sites.get_df.return_value = cutting_sites_df
        masked_df = mbias_data.get_masked_mbias_df(trimming_mode='minimal', cutting_sites=cutting_sites)
        assert np.isnan(masked_df.loc[('c_bc', 1, 1), 'n_meth'])
        assert np.isnan(masked_df.loc[('w_bc', 1, 3), 'n_meth'])
        assert masked_df.loc[('c_bc', 1, 3), 'n_meth'] == 1
        assert masked_df.loc[('w_bc', 1, 2), 'n_meth'] == 1


class TestMbiasPlotter:
    @pytest.fixture(scope='class')
    def config_plotter(self, config):
        config['data_properties']['max_read_length_bp'] = 3
        config['trimming']['flens_to_show_in_plots'] = [1,2,3,4,5]
        return config

    def test_unmasked_flen_plot(self, mbias_stats_df, test_output_dir,
                                config_plotter):
        mbias_stats_df['beta_values_smoothed'] = mbias_stats_df['beta_values']
        mbias_data = Mock(mbias_stats_df=mbias_stats_df,
                          min_flen_considered_for_methylation_calling = 1,
                          max_flen_considered_for_trimming = 5,
                          spec=mqc.mbias.MbiasData
                          )

        mbias_plotter = mqc.mbias.MbiasDataPlotter(mbias_data, config_plotter)
        mbias_plotter.distance_between_displayed_flen = 1
        mbias_plotter.min_mean_raw_data = 4

        output_path = os.path.join(test_output_dir, 'unmasked_mbias_plot.png')
        mbias_plotter.flen_strat_plot(output_path=output_path,
                                      plot_smoothed_values=True)
        output_path = os.path.join(test_output_dir, 'unmasked_mbias_plot.raw.png')
        mbias_plotter.flen_strat_plot(output_path=output_path,
                                      plot_smoothed_values=False)

    def test_masked_flen_plot(self, mbias_stats_df, test_output_dir,
                              config_plotter):

        idx = pd.IndexSlice

        # mbias_stats_df index levels:
        # bsseq_strand, flen, pos
        mbias_stats_df['beta_values_smoothed'] = mbias_stats_df['beta_values']
        mbias_stats_df.loc[idx[:, :, 1], :] = np.nan

        mbias_data = Mock(mbias_stats_df=None,
                          min_flen_considered_for_methylation_calling = 1,
                          max_flen_considered_for_trimming = 5,
                          spec=mqc.mbias.MbiasData
                          )

        mbias_data.get_masked_mbias_df.return_value = mbias_stats_df

        mbias_plotter = mqc.mbias.MbiasDataPlotter(mbias_data, config_plotter)
        mbias_plotter.distance_between_displayed_flen = 1
        mbias_plotter.min_mean_raw_data = 4

        output_path = os.path.join(test_output_dir, 'masked_mbias_plot.png')
        mbias_plotter.flen_strat_plot(output_path=output_path,
                                      cutting_sites=Mock(),
                                      plot_smoothed_values=True)

    # Note: this test will not fail when run as part of the test suite
    # HOWEVER, it will only produce a plot if it is run alone
    def test_total_plot(self,
                        ya_mbias_stats_df,
                        config, test_output_dir):

        df = ya_mbias_stats_df
        df.loc[df.eval('pos < 20 or pos > 80'), :] = np.nan

        mbias_data = Mock()
        mbias_data.get_masked_mbias_df.return_value = df

        mbias_plotter = mqc.mbias.MbiasDataPlotter(mbias_data, config)
        output_path = os.path.join(test_output_dir,
                                   'total_mbias_plot.masked.png')
        mbias_plotter.total_plot(output_path, adjusted_cutting_sites=Mock())


class TestAdjustedMbiasCuttingSites:
    def test_get_array(self, config, cutting_sites_df):
        adjusted_cutting_sites = mqc.mbias.AdjustedMbiasCuttingSites(
            mbias_data=Mock(),
            calling_mode='standard',
            config=config
        )
        adjusted_cutting_sites._adjusted_cutting_sites_df = cutting_sites_df
        arr = adjusted_cutting_sites.get_array()

        # bsseq_strand	flen	start	end
        # c_bc	4	2	3
        assert arr[0, 0, 4] == 2
        assert arr[0, 1, 4] == 3
        # w_bc	1	1	2
        assert arr[2, 0, 1] == 1
        assert arr[2, 1, 1] == 2

    def test_fit_normalvariate_plateau(
            self, mbias_stats_df_large_for_cutting_site_determin_tests, config):
        """Test and benchmark normalvariate-based plateau calling

        Test data curves are bad for flen 100, good for flen 101 with plateau from 21 to 80

        This only tests the plateau calling function directly, which acts on individual
        M-bias curves for a given (bsseq_strand, flen) stratum

        The benchmarks are extrapolations on the expected duration across all bsseq_strand, flen strata
        """

        df = mbias_stats_df_large_for_cutting_site_determin_tests

        def check_curve(curve_df, correct_cutting_sites, curve_name):
            t0 = time.time()
            ser = mqc.mbias.AdjustedMbiasCuttingSites.fit_normalvariate_plateau(curve_df, config)
            t1 = time.time()
            total_time = t1-t0
            if curve_name == 'good_curve':
                assert (ser[['start', 'end']] - correct_cutting_sites).abs().sum() < 8
            elif curve_name == 'bad_curve':
                assert (ser[['start', 'end']] == [0,0]).all()
            else:
                raise NotImplementedError
            print(f'Timing for {curve_name}')
            print(f'Individual run took {total_time}s')
            n_strands = 4
            n_flens = 500
            time_across_all_required_iterations_min = n_strands * n_flens * total_time / 60
            print(f'This is equivalent to a total time of {time_across_all_required_iterations_min}m')


        # good M-bias curve
        idx = pd.IndexSlice
        curve_df = df.loc[idx['w_bc':'w_bc', 101:101], :]
        check_curve(curve_df, [21, 80], 'good_curve')

        # bad M-bias curve
        curve_df = df.loc[idx['w_bc':'w_bc', 100:100], :]
        check_curve(curve_df, [0, 0], 'bad_curve')


    def test_max_plateau_length_threshold_std_cutting_site_calling(self, config,
                                                                   mbias_stats_df_large_for_cutting_site_determin_tests):
        """ Test determination of cutting sites df

        Test uses representative plateau calling function (normalvariate)

        Test data curves are bad for flen 100, good for flen 101 with plateau from 21 to 80
        """
        mbias_data = Mock(mbias_stats_df=mbias_stats_df_large_for_cutting_site_determin_tests)
        adjusted_cutting_sites = mqc.mbias.AdjustedMbiasCuttingSites(
            mbias_data=mbias_data,
            calling_mode='standard',
            config=config
        )
        cutting_sites_df = adjusted_cutting_sites.get_df()
        assert (cutting_sites_df.loc[('c_bc', 100), ['start', 'end']] == [0,0]).all()
        assert (cutting_sites_df.loc[('c_bc', 101), ['start', 'end']] - [21,80]).abs().sum() < 8
