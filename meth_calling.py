import gzip
import mqc
import mqc.pileup_methylation_calling as mc
import mqc.trimming
import pytoml


def main():
    with open('./config.default.toml') as f_toml:
        config_dict = pytoml.load(f_toml)

    import time
    t0 = time.time()
    print('Working')

    cutting_site_array = mqc.trimming.cutting_sites_array_from_flen_relative_minimal_cutting_sites(
        relative_cutting_site_dict=config_dict['trimming']['relative_to_fragment_ends'],
        max_flen_considered_for_trimming=config_dict['trimming']['max_flen_considered_for_trimming'],
        max_read_length_bp=config_dict['data_properties']['max_read_length_bp']
    )

    # tabulate_meth_calls(bam_path='/home/kraemers/projects/mqc/mqc/test/data/b_cells_rep1_chr11_16815793-16824254.bam',
    #                     index_file_path='./test/data/chr11_16815793-16824254.cg.bed.gz',
    #                     output_file='./test/results/methylation_calls_new.bed.gz',
    #                     cutting_site_array=cutting_site_array,
    #                     config_dict=config_dict)

    tabulate_meth_calls(bam_path=('/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs'
                                  '/results_per_pid/hsc_rep2/alignment/blood_hsc_rep2_merged.mdup.bam'),
                        index_file_path='/home/kraemers/projects/mqc/mqc/test/data/chr11.100000.cg.bed.gz',
                        cutting_site_array=cutting_site_array,
                        output_file='./test/results/methylation_calls_new.bed.gz',
                        config_dict=config_dict)

    t1 = time.time()
    total = t1 - t0
    print('Time for 100,000 positions(s): ', total)
    print('Projected time(h): ', total / 100000 * 30 * 10 ** 6 / 60 / 60)


def tabulate_meth_calls(bam_path, index_file_path, output_file, cutting_site_array, config_dict):
    motif_pileup_iter = mqc.motif_pileup_generator(bam_path, index_file_path)

    with gzip.open(output_file, 'wt') as fout:
        for motif_pileups, curr_idx_pos in motif_pileup_iter:
            beta, n_meth, n_unmeth = mc.call_meth_at_pileup(
                motif_pileups,
                index_position=curr_idx_pos,
                cutting_site_array=cutting_site_array)

            print(curr_idx_pos.chrom,
                  curr_idx_pos.start,
                  curr_idx_pos.end,
                  beta,
                  n_meth,
                  n_unmeth,
                  sep='\t', end='\n', file=fout)


if __name__ == '__main__':
    main()
