import gzip
import shutil

import pickle

import mqc
import os.path


# TODO: there was a numpy warning: numpy-division-with-runtimewarning-invalid-value-encountered-in-double-scalars
# I think it was during beta value calculation in the second pass
# The warning went away when I scaled down the test size by changing some config values
# -> Reproduce the error with the normal config values and find out what is causing it
# based on this SO link, this error really makes no sense during my beta value calculation
# http://stackoverflow.com/questions/27784528/numpy-division-with-runtimewarning-invalid-value-encountered-in-double-scalars


def qc_run(bam_path, index_file_path, config, meth_metrics_dir_abs, sample_name, test=False):
    """Evaluate QC stats with minimal and adjusted trimming

    Runner function for the complete quality control analysis, including
    - beta value distribution statistics
    - all M-bias plots
    - coverage statistics
    - conversion error statistics

    If QC results from a previous run exist and the new QC stats are written to the old methylation metrics
    directory, all old results will be deleted
    """

    # TODO-format: remove this!
    if test:
        config['trimming']['min_flen_considered_for_methylation_calling'] = 100
        config['trimming']['max_flen_considered_for_trimming'] = 110
        config['trimming']['max_window_size_for_smoothing'] = 20
        config['trimming']['required_n_events_for_cutting_site_determination'] = 1

    out_paths = _prepare_output_dirs_and_paths(meth_metrics_dir_abs, sample_name=sample_name)

    print('Calculating minimal cutting sites...')
    minimal_cutting_sites = mqc.mbias.MinimalCuttingSites(config)
    print('done\n')

    mbias_counter, beta_value_counter_minimal = _first_pass_qc_stats_collection(
            bam_path=bam_path, index_file_path=index_file_path,
            meth_calls_minimal_trimming_gz_path=out_paths['meth_calls_minimal_trimming'],
            minimal_cutting_sites=minimal_cutting_sites, config=config)

    mbias_data, mbias_adjusted_cutting_sites = _calculate_mbias_stats(mbias_counter, config)

    beta_count_adjusted_trimm = _second_pass_qc_stats_collection(
            bam_path=bam_path, index_file_path=index_file_path,
            meth_calls_adjusted_trimming_gz_path=out_paths['meth_calls_adjusted_trimming'],
            config=config, mbias_adjusted_cutting_sites=mbias_adjusted_cutting_sites)

    beta_value_data = _calculate_beta_value_distribution_stats(beta_value_counter_minimal, beta_count_adjusted_trimm)

    _plot_mbias_stats(
            mbias_data, mbias_adjusted_cutting_sites, out_paths=out_paths, config=config)

    _plot_beta_value_distribution_stats(beta_value_data, out_paths=out_paths)

    _save_data(out_paths, mbias_data, beta_value_data)


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
    out_paths = dict()
    out_paths['cutting_sites_plot'] = os.path.join(d, 'mbias_stats', sample_name + '.cutting_sites.png')

    out_paths['mbias_stats_no_mask_path'] = os.path.join(d, 'mbias_stats', sample_name + '.mbias_plot.unmasked.png')
    out_paths['mbias_stats_masked_adjusted'] = os.path.join(d, 'mbias_stats', sample_name + '.mbias_plot.masked.png')
    out_paths['mbias_stats_no_mask_path_raw'] = os.path.join(d, 'mbias_stats', sample_name + '.mbias_plot.unmasked.raw.png')
    out_paths['mbias_stats_masked_adjusted_raw'] = os.path.join(d, 'mbias_stats', sample_name + '.mbias_plot.masked.raw.png')

    out_paths['beta_value_dist'] = os.path.join(d, 'beta_values', sample_name + '.beta_value_dist.png')


    out_paths['beta_value_data_p'] = os.path.join(d, 'beta_values', sample_name + '.beta_value_data.p')
    out_paths['mbias_data_p'] = os.path.join(d, 'mbias_stats', sample_name + '.mbias_data.p')


    out_paths['meth_calls_minimal_trimming'] = os.path.join(
            d, 'sampled_mcalls', sample_name + '.mcalls.min_trimm.bed.gz')
    out_paths['meth_calls_adjusted_trimming'] = os.path.join(
            d, 'sampled_mcalls',  sample_name + '.mcalls.adjust_trimm.bed.gz')

    return out_paths


def _first_pass_qc_stats_collection(bam_path, index_file_path, meth_calls_minimal_trimming_gz_path,
                                    minimal_cutting_sites, config):
    print('Starting first pass collection...')
    mbias_counter = mqc.MbiasCounter(config)
    beta_value_counter_minimal = mqc.beta_values.StratifiedBetaValueCounter(config)
    motif_pileup_iter = mqc.motif_pileup_generator(bam_path, index_file_path)
    meth_calls_fobj_minimal = gzip.open(meth_calls_minimal_trimming_gz_path, 'wt')

    for motif_pileups, curr_idx_pos in motif_pileup_iter:
        mqc.pileup.annotate_pileupreads(motif_pileups=motif_pileups, index_position=curr_idx_pos,
                                        cutting_sites=minimal_cutting_sites,
                                        config=config)

        mbias_counter.update(motif_pileups, index_position=curr_idx_pos)

        n_meth, n_unmeth, beta = beta_value_counter_minimal.update_and_return_total_beta_value(
                motif_pileups, curr_idx_pos)

        print(n_meth, n_unmeth, beta, sep='\t', file=meth_calls_fobj_minimal)
    meth_calls_fobj_minimal.close()
    print('done\n')

    return mbias_counter, beta_value_counter_minimal


def _calculate_mbias_stats(mbias_counter, config):
    print('Calculating M-bias stats...')
    mbias_data = mqc.mbias.MbiasData(mbias_counter, config)
    mbias_data.convert_mbias_arr_info_to_df_format()
    mbias_data.add_smoothed_mbias_stats()
    mbias_data.add_beta_values()

    mbias_adjusted_cutting_sites = mqc.mbias.AdjustedMbiasCuttingSites(
            mbias_data, calling_mode='standard', config=config)

    return mbias_data, mbias_adjusted_cutting_sites
    print('Done')


def _second_pass_qc_stats_collection(bam_path, index_file_path, mbias_adjusted_cutting_sites,
                                     meth_calls_adjusted_trimming_gz_path, config):
    print('Starting second pass collection...')
    beta_count_adjusted_trimm = mqc.beta_values.StratifiedBetaValueCounter(config)
    meth_calls_fobj_adjusted = gzip.open(meth_calls_adjusted_trimming_gz_path, 'wt')

    motif_pileup_iter = mqc.motif_pileup_generator(bam_path, index_file_path)
    for motif_pileups, curr_idx_pos in motif_pileup_iter:
        mqc.pileup.annotate_pileupreads(motif_pileups=motif_pileups, index_position=curr_idx_pos,
                                        cutting_sites=mbias_adjusted_cutting_sites, config=config)
        n_meth, n_unmeth, beta = beta_count_adjusted_trimm.update_and_return_total_beta_value(
                motif_pileups, curr_idx_pos)
        print(n_meth, n_unmeth, beta, sep='\t', file=meth_calls_fobj_adjusted)
    meth_calls_fobj_adjusted.close()

    print('done\n')
    return beta_count_adjusted_trimm


def _calculate_beta_value_distribution_stats(beta_value_counter_minimal, beta_count_adjusted_trimm):
    print('Calcuting beta value distribution stats...')
    beta_value_data = mqc.beta_values.BetaValueData()
    beta_value_data.add_data_from_counter(beta_value_counter_minimal,
                                          region_str='global',
                                          trimming_status_str='minimal')
    beta_value_data.add_data_from_counter(beta_count_adjusted_trimm,
                                          region_str='global',
                                          trimming_status_str='adjusted')
    beta_value_data.compute_frequencies()
    print('done\n')
    return beta_value_data


def _plot_mbias_stats(mbias_data, mbias_adjusted_cutting_sites, out_paths, config):
    print('Plotting Mbias stats')
    mqc.mbias.cutting_sites_area_plot(mbias_cutting_sites=mbias_adjusted_cutting_sites,
                                      output_path=out_paths['cutting_sites_plot'])

    plotter = mqc.mbias.MbiasDataPlotter(mbias_data, config)
    # TODO: make this paths arguments to qc_run function
    # TODO: add temp_out_dir folder
    plotter.flen_strat_plot(output_path=out_paths['mbias_stats_no_mask_path'])
    plotter.flen_strat_plot(output_path=out_paths['mbias_stats_masked_adjusted'],
                            cutting_sites=mbias_adjusted_cutting_sites, trimming_mode='adjusted')
    plotter.flen_strat_plot(output_path=out_paths['mbias_stats_no_mask_path_raw'],
                            plot_smoothed_values=False)
    plotter.flen_strat_plot(output_path=out_paths['mbias_stats_masked_adjusted_raw'],
                            cutting_sites=mbias_adjusted_cutting_sites, trimming_mode='adjusted',
                            plot_smoothed_values=False)
    print('done\n')


def _plot_beta_value_distribution_stats(beta_value_data, out_paths):
    print('Plotting beta values')
    beta_value_plotter = mqc.beta_values.BetaValuePlotter(beta_value_data)
    beta_value_plotter.beta_value_dist_plot(region_str='global',
                                            output_path=out_paths['beta_value_dist'])
    print('done\n')


def _save_data(out_paths, mbias_data, beta_value_data):
    with open(out_paths['mbias_data_p'], 'wb') as mbias_f:
        pickle.dump(obj=mbias_data, file=mbias_f)

    with open(out_paths['beta_value_data_p'], 'wb') as beta_f:
        pickle.dump(obj=beta_value_data, file=beta_f)
