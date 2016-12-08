import pickle

import mqc
import pytoml
import sys

# Methylation status flags
IS_NA = 16
IS_METHYLATED = 8
IS_UNMETHYLATED = 4
IS_SNP = 2
IS_REF = 1

# BS-Seq strand flags
W_BC = 96
C_BC = 80
W_BC_RV = 144
C_BC_RV = 160
MATE_AND_DIR_BITS = 240

# Indices
C_BC_IND = 0
C_BC_RV_IND = 2
W_BC_IND = 4
W_BC_RV_IND = 6


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

    mode = sys.argv[1]
    if mode == "small":
        qc_run(bam_path='/home/kraemers/projects/mqc/mqc/test/data/b_cells_rep1_chr11_16815793-16824254.bam',
               index_file_path='./test/data/chr11_16815793-16824254.cg.bed.gz',
               cutting_site_array=cutting_site_array,
               config_dict=config_dict)
    else:
        qc_run(bam_path=('/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs'
                         '/results_per_pid/hsc_rep2/alignment/blood_hsc_rep2_merged.mdup.bam'),
               index_file_path='/home/kraemers/projects/mqc/mqc/test/data/chr11.100000.cg.bed.gz',
               cutting_site_array=cutting_site_array,
               config_dict=config_dict)

    t1 = time.time()
    total = t1 - t0
    print('Time for 100,000 positions(s): ', total)
    print('Projected time(h): ', total / 100000 * 30 * 10 ** 6 / 60 / 60)


def qc_run(bam_path, index_file_path, cutting_site_array, config_dict):
    mbias_counter = mqc.MbiasCounter(
        max_read_length=config_dict['data_properties']['max_read_length_bp'],
        min_phred_score=config_dict['basic_quality_filtering']['min_phred_score'],
        max_flen_considered_for_trimming=config_dict['trimming']['max_flen_considered_for_trimming']
    )

    motif_pileup_iter = mqc.motif_pileup_generator(bam_path, index_file_path)
    for motif_pileups, curr_idx_pos in motif_pileup_iter:
        annotate_pileupreads(motif_pileups=motif_pileups, index_position=curr_idx_pos,
                             cutting_site_array=cutting_site_array)
        mbias_counter.update(motif_pileups, index_position=curr_idx_pos)

    with open('./test/results/mbias_stats_array.p', 'wb') as fobj_array_dump:
        pickle.dump(mbias_counter.counter, fobj_array_dump)

    print(sum(sum(mbias_counter.counter > 0)))


def annotate_pileupreads(motif_pileups, index_position, cutting_site_array):
    watson_motif_seq = index_position.watson_motif
    for motif_base, pileup_reads in zip(watson_motif_seq, motif_pileups):
        if motif_base in ['C', 'G']:
            mqc.trimming.set_trimming_flag(pileup_reads, cutting_site_array)
            mqc.overlap.tag_overlaps(pileup_reads)


if __name__ == '__main__':
    main()
