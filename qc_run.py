import pickle

import mqc
import pytoml
import sys


def main():
    with open('./config.default.toml') as f_toml:
        config_dict = pytoml.load(f_toml)

    max_flen = config_dict['trimming']['max_flen_considered_for_trimming']

    minimal_cutting_sites = mqc.mbias.MinimalCuttingSites(
        relative_cutting_site_dict=config_dict['trimming']['relative_to_fragment_ends'],
        max_flen_considered_for_trimming=max_flen,
        max_read_length_bp=config_dict['data_properties']['max_read_length_bp']
    )

    mode = sys.argv[1]
    if mode == 'small':
        qc_run(bam_path='../test_data/alignments/b_cells_rep1_chr11_16815793-16824254.bam',
               index_file_path='../test_data/indices/chr11_16815793-16824254.cg.bed.gz',
               cutting_site_array=minimal_cutting_sites.get_array(),
               config_dict=config_dict,
               max_flen_considered_for_trimming=max_flen)

    elif mode == 'medium':
        import time
        t0 = time.time()
        qc_run(bam_path=('/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs'
                         '/results_per_pid/hsc_rep2/alignment/blood_hsc_rep2_merged.mdup.bam'),
               index_file_path='/home/kraemers/projects/mqc/mqc/test/data/chr11.10000.cg.bed.gz',
               cutting_site_array=cutting_site_array,
               config_dict=config_dict)

        t1 = time.time()
        total = t1 - t0
        print('Time for 10,000 positions(s): ', total)
        print('Projected time(h): ', total / 10000 * 30 * 10 ** 6 / 60 / 60)

    elif mode == 'big':
        import time
        t0 = time.time()
        qc_run(bam_path=('/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs'
                         '/results_per_pid/hsc_rep2/alignment/blood_hsc_rep2_merged.mdup.bam'),
               index_file_path='/home/kraemers/projects/mqc/mqc/test/data/chr11.100000.cg.bed.gz',
               cutting_site_array=cutting_site_array,
               config_dict=config_dict)

        t1 = time.time()
        total = t1 - t0
        print('Time for 100,000 positions(s): ', total)
        print('Projected time(h): ', total / 100000 * 30 * 10 ** 6 / 60 / 60)

    else:
        print('Mode unknown')


def qc_run(bam_path, index_file_path, cutting_site_array, config_dict,
           max_flen_considered_for_trimming):

    mbias_counter = mqc.MbiasCounter(
        max_read_length=config_dict['data_properties']['max_read_length_bp'],
        min_phred_score=config_dict['basic_quality_filtering']['min_phred_score'],
        max_flen_considered_for_trimming=config_dict['trimming']['max_flen_considered_for_trimming']
    )

    beta_value_counter = mqc.beta_values.BetaValueCounter(
        min_cov=config_dict['beta_value_dist_stats']['min_cov'],
    )

    motif_pileup_iter = mqc.motif_pileup_generator(bam_path, index_file_path)
    for motif_pileups, curr_idx_pos in motif_pileup_iter:
        annotate_pileupreads(motif_pileups=motif_pileups, index_position=curr_idx_pos,
                             cutting_site_array=cutting_site_array,
                             max_flen_considered_for_trimming=max_flen_considered_for_trimming)
        mbias_counter.update(motif_pileups, index_position=curr_idx_pos)
        beta_value_counter.update_and_return_total_beta_value(motif_pileups, curr_idx_pos)

    beta_value_data = mqc.beta_values.BetaValueData()
    beta_value_data.add_data_from_counter(beta_value_counter,
                                          region_str='global',
                                          trimming_status_str='minimal')
    beta_value_plotter = mqc.beta_values.BetaValuePlotter(beta_value_data)
    beta_value_plotter.beta_value_dist_plot(region_str='global',
                                            output_path='/home/kraemers/projects/mqc/test_data/'
                                                        'results/beta_value_dist_mate_str.png')

    with open('/home/kraemers/projects/mqc/test_data/results/mbias_stats_array.p', 'wb') as fobj_array_dump:
        pickle.dump(mbias_counter.counter, fobj_array_dump)

    with open('/home/kraemers/projects/mqc/test_data/results/beta_value_dist_array.p', 'wb') as fobj:
        pickle.dump(beta_value_counter, fobj)

    print('Mbias:\n')
    print(sum(sum(mbias_counter.counter > 0)))
    print('-----------------')
    print('Beta value dist:\n')
    print(beta_value_counter.beta_counter.sum(axis=1))
    print('-----------------')
    print('Beta value df:\n')
    print(beta_value_data)


def annotate_pileupreads(motif_pileups, index_position, cutting_site_array,
                         max_flen_considered_for_trimming):
    watson_motif_seq = index_position.watson_motif
    for motif_base, pileup_reads in zip(watson_motif_seq, motif_pileups):
        if motif_base in ['C', 'G']:
            mqc.trimming.set_trimming_flag(pileup_reads, cutting_site_array,
                                           max_flen_considered_for_trimming=max_flen_considered_for_trimming)
            mqc.overlap.tag_overlaps(pileup_reads)


if __name__ == '__main__':
    main()
