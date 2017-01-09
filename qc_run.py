import gzip
import mqc
import os.path

MIN_LEN = 50
MAX_STD = 10


def qc_run(bam_path, index_file_path, config_dict, output_dir_abspath, sample_name):

    max_read_length = config_dict['data_properties']['max_read_length_bp']
    min_phred_score = config_dict['basic_quality_filtering']['min_phred_score']
    max_flen = config_dict['trimming']['max_flen_considered_for_trimming']

    mbias_counter = mqc.MbiasCounter(
        max_read_length=max_read_length,
        min_phred_score=min_phred_score,
        max_flen_considered_for_trimming=max_flen
    )

    beta_value_counter_minimal = mqc.beta_values.BetaValueCounter(
        min_cov=config_dict['beta_value_dist_stats']['min_cov'],
    )

    minimal_cutting_sites = mqc.mbias.MinimalCuttingSites(
        relative_cutting_site_dict=config_dict['trimming']['relative_to_fragment_ends'],
        max_flen_considered_for_trimming=max_flen,
        max_read_length_bp=config_dict['data_properties']['max_read_length_bp']
    )
    cutting_site_array = minimal_cutting_sites.get_array()

    meth_calls_minimal_trimming_gz_path = os.path.join(output_dir_abspath, sample_name + '.mcalls.min_trimm.bed.gz')
    meth_calls_fobj = gzip.open(meth_calls_minimal_trimming_gz_path, 'wt')
    motif_pileup_iter = mqc.motif_pileup_generator(bam_path, index_file_path)
    for motif_pileups, curr_idx_pos in motif_pileup_iter:
        annotate_pileupreads(motif_pileups=motif_pileups, index_position=curr_idx_pos,
                             cutting_site_array=cutting_site_array,
                             max_flen_considered_for_trimming=max_flen)
        mbias_counter.update(motif_pileups, index_position=curr_idx_pos)
        n_meth, n_unmeth, beta = beta_value_counter_minimal.update_and_return_total_beta_value(
            motif_pileups, curr_idx_pos)
        print(n_meth, n_unmeth, beta, sep='\t', file=meth_calls_fobj)
    meth_calls_fobj.close()

    beta_value_data = mqc.beta_values.BetaValueData()
    beta_value_data.add_data_from_counter(beta_value_counter_minimal,
                                          region_str='global',
                                          trimming_status_str='minimal')

    mbias_data = mqc.mbias.MbiasData(
        max_flen_considered_for_trimming=max_flen,
        max_read_length_bp=max_read_length,
        # TODO: Ref. to take MbiasCounter instead of array
        mbias_stats_array=mbias_counter.counter,
        required_n_events_for_cutting_site_determination=(
            config_dict['trimming']['required_n_events_for_cutting_site_determination']),
        max_window_size_for_smoothing=config_dict['trimming']['max_window_size_for_smoothing'])

    mbias_adjusted_cutting_sites = mqc.mbias.AdjustedMbiasCuttingSites(
        mbias_data, calling_mode='standard', min_len=MIN_LEN, max_std=MAX_STD)

    adjusted_cutting_sites_df = mbias_adjusted_cutting_sites.get_df()
    # TODO: should accept CuttingSites
    cutting_sites_plot_path = os.path.join(output_dir_abspath, sample_name + '_adjusted_cutting_sites.png')
    mqc.mbias.cutting_sites_area_plot(mbias_cutting_sites_df=adjusted_cutting_sites_df,
                                      output_path=cutting_sites_plot_path)

    mbias_data._adjusted_cutting_sites_df = adjusted_cutting_sites_df

    plotter = mqc.mbias.MbiasDataPlotter(mbias_data=mbias_data)
    plotter.flen_strat_plot(output_path='/home/kraemers/temp/mbias_stats.no_mask.png')
    plotter.flen_strat_plot(output_path='/home/kraemers/temp/mbias_stats_masked_adjusted.png',
                            trimming_mode='adjusted')
    # Not implemented
    # plotter.flen_strat_plot(output_path='/home/kraemers/temp/mbias_stats_masked_minimal.png',
    #                         trimming_mode='minimal')

    beta_value_counter_adjusted_trimming = mqc.beta_values.BetaValueCounter(
        min_cov=config_dict['beta_value_dist_stats']['min_cov'],
    )

    meth_calls_adjusted_trimming_gz_path = os.path.join(output_dir_abspath, sample_name + '.mcalls.adjust_trimm.bed.gz')
    meth_calls_fobj = gzip.open(meth_calls_adjusted_trimming_gz_path, 'wt')
    motif_pileup_iter = mqc.motif_pileup_generator(bam_path, index_file_path)
    for motif_pileups, curr_idx_pos in motif_pileup_iter:
        annotate_pileupreads(motif_pileups=motif_pileups, index_position=curr_idx_pos,
                             cutting_site_array=cutting_site_array,
                             max_flen_considered_for_trimming=max_flen)
        n_meth, n_unmeth, beta = beta_value_counter_adjusted_trimming.update_and_return_total_beta_value(motif_pileups,
                                                                                                         curr_idx_pos)
        print(n_meth, n_unmeth, beta, sep='\t', file=meth_calls_fobj)
    meth_calls_fobj.close()

    beta_value_data.add_data_from_counter(beta_value_counter_adjusted_trimming,
                                          region_str='global',
                                          trimming_status_str='adjusted')

    beta_value_plotter = mqc.beta_values.BetaValuePlotter(beta_value_data)
    beta_value_plot_path = os.path.join(output_dir_abspath, sample_name + 'beta_value_dist.png')
    beta_value_plotter.beta_value_dist_plot(region_str='global',
                                            output_path=beta_value_plot_path)

    # with open('/home/kraemers/projects/mqc/test_data/results/mbias_stats_array.p', 'wb') as fobj_array_dump:
    #     pickle.dump(mbias_counter.counter, fobj_array_dump)
    #
    # with open('/home/kraemers/projects/mqc/test_data/results/beta_value_dist_array.p', 'wb') as fobj:
    #     pickle.dump(beta_value_counter, fobj)

    # print('Mbias:\n')
    # print(sum(sum(mbias_counter.counter > 0)))
    # print('-----------------')
    # print('Beta value dist:\n')
    # print(beta_value_counter.beta_counter.sum(axis=1))
    # print('-----------------')
    # print('Beta value df:\n')
    # print(beta_value_data)


def annotate_pileupreads(motif_pileups, index_position, cutting_site_array,
                         max_flen_considered_for_trimming):
    watson_motif_seq = index_position.watson_motif
    for motif_base, pileup_reads in zip(watson_motif_seq, motif_pileups):
        if motif_base in ['C', 'G']:
            mqc.trimming.set_trimming_flag(pileup_reads, cutting_site_array,
                                           max_flen_considered_for_trimming=max_flen_considered_for_trimming)
            mqc.overlap.tag_overlaps(pileup_reads)
