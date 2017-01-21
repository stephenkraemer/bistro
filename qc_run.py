import gzip
import mqc
import os.path
# import logging



# TODO: there was a numpy warning: numpy-division-with-runtimewarning-invalid-value-encountered-in-double-scalars
# I think it was during beta value calculation in the second pass
# The warning went away when I scaled down the test size by changing some config values
# -> Reproduce the error with the normal config values and find out what is causing it
# based on this SO link, this error really makes no sense during my beta value calculation
# http://stackoverflow.com/questions/27784528/numpy-division-with-runtimewarning-invalid-value-encountered-in-double-scalars
def qc_run(bam_path, index_file_path, config, output_dir_abspath, sample_name):

    # logger = logging.getLogger(__name__)
    # logger.setLevel(logging.DEBUG)

    os.makedirs(output_dir_abspath, mode=0o770, exist_ok=True)

    print('Calculating minimal cutting sites...')
    minimal_cutting_sites = mqc.mbias.MinimalCuttingSites(config)
    print('done\n')

    meth_calls_minimal_trimming_gz_path = os.path.join(output_dir_abspath, sample_name + '.mcalls.min_trimm.bed.gz')
    meth_calls_adjusted_trimming_gz_path = os.path.join(output_dir_abspath, sample_name + '.mcalls.adjust_trimm.bed.gz')
    beta_value_plot_path = os.path.join(output_dir_abspath, sample_name + 'beta_value_dist.png')


    print('Starting first pass collection...')
    mbias_counter, beta_value_counter_minimal =_first_pass_qc_stats_collection(
        bam_path=bam_path, index_file_path=index_file_path,
        meth_calls_minimal_trimming_gz_path=meth_calls_minimal_trimming_gz_path,
        minimal_cutting_sites=minimal_cutting_sites, config=config)
    print('done\n')

    print('Calculating M-bias stats...')
    mbias_data, mbias_adjusted_cutting_sites = _calculate_mbias_stats(mbias_counter, config)
    print('done\n')

    print('Starting second pass collection...')
    beta_count_adjusted_trimm = _second_pass_qc_stats_collection(
        bam_path=bam_path, index_file_path=index_file_path,
        meth_calls_adjusted_trimming_gz_path=meth_calls_adjusted_trimming_gz_path,
        config=config, mbias_adjusted_cutting_sites=mbias_adjusted_cutting_sites)
    print('done\n')

    print('Calcuting beta value distribution stats...')
    beta_value_data = _calculate_beta_value_distribution_stats(beta_value_counter_minimal, beta_count_adjusted_trimm)
    print('done\n')

    print('Plotting...')
    _plot_mbias_stats(
        mbias_data, mbias_adjusted_cutting_sites, output_dir=output_dir_abspath, config=config)

    _plot_beta_value_distribution_stats(beta_value_data, beta_value_plot_path)
    print('done\n')

    # _save_objects(mbias_data, beta_value_data)


def _first_pass_qc_stats_collection(bam_path, index_file_path, meth_calls_minimal_trimming_gz_path,
                                    minimal_cutting_sites, config):

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

    return mbias_counter, beta_value_counter_minimal

def _calculate_mbias_stats(mbias_counter, config):
    mbias_data = mqc.mbias.MbiasData(mbias_counter, config)
    mbias_data.convert_mbias_arr_info_to_df_format()
    mbias_data.add_smoothed_mbias_stats()
    mbias_data.add_beta_values_from_smoothed_data()

    mbias_adjusted_cutting_sites = mqc.mbias.AdjustedMbiasCuttingSites(
        mbias_data, calling_mode='standard', config=config)

    return mbias_data, mbias_adjusted_cutting_sites

def _second_pass_qc_stats_collection(bam_path, index_file_path, mbias_adjusted_cutting_sites,
                                     meth_calls_adjusted_trimming_gz_path, config):

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

    return beta_count_adjusted_trimm


def _calculate_beta_value_distribution_stats(beta_value_counter_minimal, beta_count_adjusted_trimm):
    beta_value_data = mqc.beta_values.BetaValueData()
    beta_value_data.add_data_from_counter(beta_value_counter_minimal,
                                          region_str='global',
                                          trimming_status_str='minimal')
    beta_value_data.add_data_from_counter(beta_count_adjusted_trimm,
                                          region_str='global',
                                          trimming_status_str='adjusted')
    beta_value_data.compute_frequencies()
    return beta_value_data


def _plot_mbias_stats(mbias_data, mbias_adjusted_cutting_sites, output_dir, config):
    mqc.mbias.cutting_sites_area_plot(mbias_cutting_sites=mbias_adjusted_cutting_sites,
                                      output_path=os.path.join(output_dir, 'cutting_sites.png'))

    plotter = mqc.mbias.MbiasDataPlotter(mbias_data, config)
    # TODO: make this paths arguments to qc_run function
    # TODO: add temp_out_dir folder
    plotter.flen_strat_plot(output_path=os.path.join(output_dir, 'mbias_stats.no_mask.png'))
    plotter.flen_strat_plot(output_path=os.path.join(output_dir, 'mbias_stats_masked_adjusted.png'),
                            trimming_mode='adjusted')
    # Not implemented
    # plotter.flen_strat_plot(output_path='/home/kraemers/temp/mbias_stats_masked_minimal.png',
    #                         trimming_mode='minimal')


def _plot_beta_value_distribution_stats(beta_value_data, beta_value_plot_path):

    beta_value_plotter = mqc.beta_values.BetaValuePlotter(beta_value_data)
    beta_value_plotter.beta_value_dist_plot(region_str='global',
                                            output_path=beta_value_plot_path)


# def _save_objects(mbias_data, beta_value_data):
    # with open('/home/kraemers/projects/mqc/test_data/results/mbias_stats_array.p', 'wb') as fobj_array_dump:
    #     pickle.dump(mbias_counter.counter, fobj_array_dump)
    #
    # with open('/home/kraemers/projects/mqc/test_data/results/beta_value_dist_array.p', 'wb') as fobj:
    #     pickle.dump(beta_value_counter, fobj)
    # pass
