import gzip
import os.path
import pickle
import shutil

import mqc


# TODO: there was a numpy warning: numpy-division-with-runtimewarning-invalid-value-encountered-in-double-scalars
# I think it was during beta value calculation in the second pass
# The warning went away when I scaled down the test size by changing some config values
# -> Reproduce the error with the normal config values and find out what is causing it
# based on this SO link, this error really makes no sense during my beta value calculation
# http://stackoverflow.com/questions/27784528/numpy-division-with-runtimewarning-invalid-value-encountered-in-double-scalars


def qc_run(bam_path, index_file_path, config, meth_metrics_dir_abs, sample_name):
    """Evaluate QC stats with minimal and adjusted trimming

    Runner function for the complete quality control analysis, including
    - beta value distribution statistics
    - all M-bias plots
    - coverage statistics
    - conversion error statistics

    If QC results from a previous run exist and the new QC stats are written to the old methylation metrics
    directory, all old results will be deleted
    """

    out_paths = _prepare_output_dirs_and_paths(meth_metrics_dir_abs, sample_name=sample_name)

    print('Calculating minimal cutting sites...')
    minimal_cutting_sites = mqc.counters.mbias.MinimalCuttingSites(config)
    print('done\n')

    mbias_counter, beta_value_counter_minimal, coverage_counter_minimal = (
        _first_pass_qc_stats_collection(
            bam_path=bam_path, index_file_path=index_file_path,
            meth_calls_minimal_trimming_gz_path=out_paths['meth_calls_minimal_trimming'],
            minimal_cutting_sites=minimal_cutting_sites, config=config))

    mbias_data, mbias_adjusted_cutting_sites = _calculate_mbias_stats(mbias_counter, config)

    beta_count_adjusted_trimm, coverage_counter_adjusted = _second_pass_qc_stats_collection(
            bam_path=bam_path, index_file_path=index_file_path,
            meth_calls_adjusted_trimming_gz_path=out_paths['meth_calls_adjusted_trimming'],
            config=config, mbias_adjusted_cutting_sites=mbias_adjusted_cutting_sites)

    beta_value_data = _calculate_beta_value_distribution_stats(
            beta_value_counter_minimal, beta_count_adjusted_trimm, config)

    coverage_data = mqc.counters.coverage.CoverageData(config)
    coverage_data.add_counter(coverage_counter_minimal)
    coverage_data.add_counter(coverage_counter_adjusted)
    coverage_data.calculate_relative_frequencies()

    _plot_mbias_stats(
            mbias_data, mbias_adjusted_cutting_sites, out_paths=out_paths, config=config)

    _plot_beta_value_distribution_stats(
            beta_value_data, out_paths=out_paths, config=config)

    _calculate_plot_save_coverage_stats(coverage_data, config, out_paths)

    _save_data(out_paths, mbias_data, beta_value_data, coverage_data)


def _prepare_output_dirs_and_paths(meth_metrics_dir_abs, sample_name) -> dict:

    if os.path.exists(meth_metrics_dir_abs):
        shutil.rmtree(meth_metrics_dir_abs)
    os.makedirs(meth_metrics_dir_abs, mode=0o770)
    print(f'Saving QC results in {meth_metrics_dir_abs}')

    d = meth_metrics_dir_abs

    os.mkdir(os.path.join(d, 'beta_values'))
    os.mkdir(os.path.join(d, 'mbias_stats'))
    os.mkdir(os.path.join(d, 'coverage'))
    os.mkdir(os.path.join(d, 'conversion'))
    os.mkdir(os.path.join(d, 'sampled_mcalls'))

    # TODO-format: refactor out_paths[] creation
    def create_path(category_str, suffix):
        return os.path.join(d, category_str, sample_name + suffix)

    out_paths = {
        'cutting_sites_plot': create_path('mbias_stats', '.cutting_sites.png'),
        'mbias_stats_no_mask_path': create_path('mbias_stats', '.mbias_plot.unmasked.png'),
        'mbias_stats_masked_adjusted': create_path('mbias_stats', '.mbias_plot.masked.png'),
        'mbias_stats_no_mask_path_raw': create_path('mbias_stats', '.mbias_plot.unmasked.raw.png'),
        'mbias_stats_masked_adjusted_raw': create_path('mbias_stats', '.mbias_plot.masked.raw.png'),
        'mbias_total.adjusted_trimming': create_path('mbias_stats', 'mbias_total.adjusted_trimming.png'),
        'mbias_total.minimal_trimming': create_path('mbias_stats', 'mbias_total.minimal_trimming.png'),
        'mbias_data_p': create_path('mbias_stats', '.mbias_data.p'),

        'beta_value_dist_basepath': create_path('beta_values', '.beta_value_dist'),
        'beta_value_data_p': create_path('beta_values', '.beta_value_data.p'),

        'coverage_cpg_freq': create_path('coverage', '.coverage.freq.cpg.png'),
        'coverage_cpg_count': create_path('coverage', '.coverage.count.cpg.png'),
        'coverage_data_p': create_path('coverage', '.coverage.data.p'),

        'meth_calls_minimal_trimming': create_path('sampled_mcalls', '.mcalls.min_trimm.bed.gz'),
        'meth_calls_adjusted_trimming': create_path('sampled_mcalls', '.mcalls.adjust_trimm.bed.gz'),
    }

    return out_paths


def _first_pass_qc_stats_collection(bam_path, index_file_path, meth_calls_minimal_trimming_gz_path,
                                    minimal_cutting_sites, config):
    print('Starting first pass collection...')
    mbias_counter = mqc.MbiasCounter(config)
    coverage_counter_minimal = mqc.CoverageCounter(config, trimming_status='minimal')
    beta_value_counter_minimal = mqc.counters.beta_values.StratifiedBetaValueCounter(config)
    meth_calls_fobj_minimal = gzip.open(meth_calls_minimal_trimming_gz_path, 'wt')

    motif_pileup_iter = mqc.motif_pileup_generator(bam_path, index_file_path)
    for motif_pileups, curr_idx_pos in motif_pileup_iter:
        mqc.pileup.annotate_pileupreads(motif_pileups=motif_pileups, index_position=curr_idx_pos,
                                        cutting_sites=minimal_cutting_sites,
                                        config=config)

        mbias_counter.process(motif_pileups, index_position=curr_idx_pos)

        n_meth, n_unmeth, beta = beta_value_counter_minimal.update_and_return_total_beta_value(
                motif_pileups, curr_idx_pos)

        per_cpg_cov = n_meth + n_unmeth
        coverage_counter_minimal.update(per_cpg_cov)

        print(n_meth, n_unmeth, beta, sep='\t', file=meth_calls_fobj_minimal)
    meth_calls_fobj_minimal.close()
    print('done\n')

    return mbias_counter, beta_value_counter_minimal, coverage_counter_minimal


def _calculate_mbias_stats(mbias_counter, config):
    print('Calculating M-bias stats...')
    mbias_data = mqc.counters.mbias.MbiasData(mbias_counter, config)
    mbias_data.convert_mbias_arr_info_to_df_format()
    mbias_data.add_smoothed_mbias_stats()

    mbias_adjusted_cutting_sites = mqc.counters.mbias.AdjustedMbiasCuttingSites(
            mbias_data, calling_mode='standard', config=config)

    print('Done')
    return mbias_data, mbias_adjusted_cutting_sites


def _second_pass_qc_stats_collection(bam_path, index_file_path, mbias_adjusted_cutting_sites,
                                     meth_calls_adjusted_trimming_gz_path, config):
    print('Starting second pass collection...')
    beta_count_adjusted_trimm = mqc.counters.beta_values.StratifiedBetaValueCounter(config)
    coverage_counter_adjusted = mqc.CoverageCounter(config, trimming_status='adjusted')
    meth_calls_fobj_adjusted = gzip.open(meth_calls_adjusted_trimming_gz_path, 'wt')

    motif_pileup_iter = mqc.motif_pileup_generator(bam_path, index_file_path)
    for motif_pileups, curr_idx_pos in motif_pileup_iter:
        mqc.pileup.annotate_pileupreads(motif_pileups=motif_pileups, index_position=curr_idx_pos,
                                        cutting_sites=mbias_adjusted_cutting_sites, config=config)
        n_meth, n_unmeth, beta = beta_count_adjusted_trimm.update_and_return_total_beta_value(
                motif_pileups, curr_idx_pos)

        per_cpg_cov = n_meth + n_unmeth
        coverage_counter_adjusted.update(per_cpg_cov)

        print(n_meth, n_unmeth, beta, sep='\t', file=meth_calls_fobj_adjusted)
    meth_calls_fobj_adjusted.close()

    print('done\n')
    return beta_count_adjusted_trimm, coverage_counter_adjusted


def _calculate_beta_value_distribution_stats(
        beta_value_counter_minimal, beta_count_adjusted_trimm, config):
    print('Calcuting beta value distribution stats...')
    beta_value_data = mqc.counters.beta_values.BetaValueData(config)
    beta_value_data.add_data_from_counter(beta_value_counter_minimal,
                                          region_str='global',
                                          trimming_status_str='minimal')
    beta_value_data.add_data_from_counter(beta_count_adjusted_trimm,
                                          region_str='global',
                                          trimming_status_str='adjusted')
    beta_value_data.add_smoothed_beta_value_frequencies()
    print('done\n')
    return beta_value_data


def _plot_mbias_stats(mbias_data, mbias_adjusted_cutting_sites, out_paths, config):
    print('Plotting Mbias stats')
    mqc.counters.mbias.cutting_sites_area_plot(mbias_cutting_sites=mbias_adjusted_cutting_sites,
                                               output_path=out_paths['cutting_sites_plot'])

    plotter = mqc.counters.mbias.MbiasDataPlotter(mbias_data, config)
    # TODO: make this paths arguments to qc_run function
    # TODO: add temp_out_dir folder
    plotter.flen_strat_plot(output_path=out_paths['mbias_stats_no_mask_path'])
    plotter.flen_strat_plot(output_path=out_paths['mbias_stats_masked_adjusted'],
                            cutting_sites=mbias_adjusted_cutting_sites)
    plotter.flen_strat_plot(output_path=out_paths['mbias_stats_no_mask_path_raw'],
                            plot_smoothed_values=False)
    plotter.flen_strat_plot(output_path=out_paths['mbias_stats_masked_adjusted_raw'],
                            cutting_sites=mbias_adjusted_cutting_sites,
                            plot_smoothed_values=False)
    plotter.total_plot(output_path=out_paths['mbias_total.adjusted_trimming'],
                       adjusted_cutting_sites=mbias_adjusted_cutting_sites)
    plotter.total_plot(output_path=out_paths['mbias_total.minimal_trimming'],
                       adjusted_cutting_sites=None)
    print('done\n')


def _plot_beta_value_distribution_stats(beta_value_data, out_paths, config):
    print('Plotting beta values')
    beta_value_plotter = mqc.counters.beta_values.BetaValuePlotter(
            beta_value_data, config)
    beta_value_plotter.beta_value_dist_plot(
            out_basepath_abs=out_paths['beta_value_dist_basepath'])
    print('done\n')


def _calculate_plot_save_coverage_stats(coverage_data, config, out_paths):
    coverage_data.save_aggregate_stats()
    coverage_data.save_counts()
    plotter = mqc.counters.coverage.CoveragePlotter(coverage_data, config)
    plotter.plot_cov_histogram(out_paths['coverage_cpg_freq'], show_frequency=True)
    plotter.plot_cov_histogram(out_paths['coverage_cpg_count'], show_frequency=False)


def _save_data(out_paths, mbias_data, beta_value_data, coverage_data):

    def dump(obj, path_name):
        with open(out_paths[path_name], 'wb') as f:
            pickle.dump(obj=obj, file=f)
    dump(mbias_data, 'mbias_data_p')
    dump(beta_value_data, 'beta_value_data_p')
    dump(coverage_data, 'coverage_data_p')
