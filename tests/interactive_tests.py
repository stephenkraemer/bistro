print('reading')
import mqc

index_file_path = '/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/mcall_qc_test_data/indices/chr11.100_000.cg.bed.gz'
bam_path = '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/results_per_pid/mpp1_rep3/alignment/blood_mpp1_rep3_merged.mdup.bam'



args_dict = {'sample_name': 'hsc_rep1',
             'sample_meta': 'population=hsc,replicate=1',
             'output_dir': '/home/kraemers/temp/'}
config = mqc.config.get_config_dict(args_dict,
                                    config_file_path=None)


mbias_counter = mqc.MbiasCounter(config)
motif_pileup_iter = mqc.motif_pileup_generator(bam_path, index_file_path)
for motif_pileups, curr_idx_pos in motif_pileup_iter:
    mbias_counter.process(motif_pileups, index_position=curr_idx_pos)


mbias_data = mqc.counters.mbias.MbiasData(mbias_counter, config)
mbias_data.convert_mbias_arr_info_to_df_format()
plotter = mqc.counters.mbias.MbiasDataPlotter(mbias_data, config)
plotter.flen_strat_plot(output_path='/home/kraemers/temp/mbias.svg',
                        cutting_sites=None,
                        plot_smoothed_values=False)
