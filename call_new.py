""" Parallelized methylation calling"""
import gzip
import math
import itertools
import pysam
import mqc
import multiprocessing as mp
from mqc.methylation_calling import get_total_methylation_stats_at_pileup
import mqc.methylation_calling.segment_methylation

bs_strand_dict_directional_protocol = {1: {'forward': 'W-BC',
                                           'reverse': 'C-BC'},
                                       2: {'forward': 'C-BC-Rv',
                                           'reverse': 'W-BC-Rv'}}

max_number_of_unconv_cytosines = 2

# Trimming site arrays will have to be defined per BS-Seq strand
# is the concept of BS-Seq strands meaningful in the undirectional protocol?
trimming_site_array = [[0] * 1000, [110] * 1000]


def parallel_call(bam_file_path, index_file_path, n_cores, output_file):
    # start a process for every piece, map results
    # write results to file
    index_file = mqc.IndexFile(index_file_path)
    index_positions = [x for x in index_file]
    chunk_size = math.ceil(len(index_positions) / n_cores)
    slices = [[i * chunk_size, (i + 1) * chunk_size] for i in range(n_cores)]
    if slices[-1][1] >= len(index_positions):
        slices[-1][1] = len(index_positions) - 1
    regions = [('11', index_positions[start].start, index_positions[end].end)
               for start, end in slices]
    pool = mp.Pool(n_cores)
    slice_iters = [iter(index_positions[start:end]) for start, end in slices]
    arg_tuples = zip([bam_file_path] * len(slices), regions, slice_iters)
    # TODO: pool.apply may be better suited
    list_of_results = pool.map(run_on_chunk, arg_tuples)
    with gzip.open(output_file, 'wt') as fout:
        for line in itertools.chain.from_iterable(list_of_results):
            fout.write(line + '\n')


def run_on_chunk(arg_tuple):
    bam_file_path, region, index_slice_iter = arg_tuple
    chrom, start, end = region
    sam_file = pysam.AlignmentFile(bam_file_path)
    idx_pos = next(index_slice_iter)
    pileup_iter = sam_file.pileup(chrom, start, end, truncate=True)
    # start is given in 1-based coordinates. we need it zero-based and will
    # increment by one at the beginning of the loop
    lines = []
    region_pos_zero_based = start - 1
    for pileup_column in pileup_iter:
        region_pos_zero_based += 1
        # TODO: don't skip index_slice_iter header
        # if it has no header, otherwise first CpG will be lost
        if region_pos_zero_based == idx_pos.start:
            # print('I think I am at ', region_pos_zero_based)
            # print('I am at ', pileup_column.reference_pos)
            # initialize new motif pileup
            curr_base = idx_pos.watson_motif[0]
            if curr_base in ['C', 'G']:
                motif_pileups = [pileup_column.pileups]
            else:
                motif_pileups = [[]]
            for curr_base in idx_pos.watson_motif[1:]:
                pileup_column = next(pileup_iter)
                region_pos_zero_based += 1
                # print('I think I am at ', region_pos_zero_based)
                if curr_base in ['C', 'G']:
                    motif_pileups.append(pileup_column.pileups)
                else:
                    motif_pileups.append([])

            beta, n_meth, n_unmeth = get_total_methylation_stats_at_pileup(
                motif_pileups,
                index_position=idx_pos,
                max_number_of_unconv_cytosines=max_number_of_unconv_cytosines,
                trimming_site_array=trimming_site_array)

            # print(beta, n_meth, n_unmeth, sep='\t', end='\n')
            # print(beta, n_meth, n_unmeth, sep='\t', end='\n', file=fout)
            lines.append('\t'.join([str(beta), str(n_meth), str(n_unmeth)]))
            try:
                idx_pos = next(index_slice_iter)
            except StopIteration:
                break
    return lines


def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--index_file')
    parser.add_argument('--bam', required=True)
    parser.add_argument('--output_file', required=True)
    parser.add_argument('--region', default=None)
    parser.add_argument('--cores', type=int, default=1)
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()

    import time

    print('1 core')
    t0 = time.time()
    parallel_call(bam_file_path=args.bam,
                  index_file_path=args.index_file,
                  output_file=args.output_file,
                  n_cores=1)
    t1 = time.time()
    total = t1-t0
    print(total * 1000)

    print('2 cores')
    t0 = time.time()
    parallel_call(bam_file_path=args.bam,
                  index_file_path=args.index_file,
                  output_file=args.output_file,
                  n_cores=2)
    t1 = time.time()
    total = t1-t0
    print(total * 1000)
"""
def call(bam_file, region, index_file_path, output_file):
    sam_file = pysam.AlignmentFile(bam_file)
    index_file = mqc.IndexFile(index_file_path)
    idx_pos = index_file.__next__()
    if region:
        chrom_str, start_str, end_str = re.split('[:-]', args.region)
        start = int(start_str)
    with gzip.open(output_file, 'wt') as fout:
        pileup_iter = sam_file.pileup(region=region, truncate=True)
        # start is given in 1-based coordinates. we need it zero-based and will
        # increment by one at the beginning of the loop
        region_pos_zero_based = start - 2 if region else -1
        for pileup_column in pileup_iter:
            region_pos_zero_based += 1
            # TODO: don't skip index header if it has no header, otherwise first CpG will be lost
            if region_pos_zero_based == idx_pos.start:
                # print('I think I am at ', region_pos_zero_based)
                # print('I am at ', pileup_column.reference_pos)
                # initialize new motif pileup
                curr_base = idx_pos.watson_motif[0]
                if curr_base in ['C', 'G']:
                    motif_pileups = [pileup_column.pileups]
                else:
                    motif_pileups = [[]]
                for curr_base in idx_pos.watson_motif[1:]:
                    pileup_column = next(pileup_iter)
                    region_pos_zero_based += 1
                    # print('I think I am at ', region_pos_zero_based)
                    if curr_base in ['C', 'G']:
                        motif_pileups.append(pileup_column.pileups)
                    else:
                        motif_pileups.append([])

                beta, n_meth, n_unmeth = call_methylation(motif_pileups,
                                                          idx_pos.watson_motif)

                # print(beta, n_meth, n_unmeth, sep='\t', end='\n')
                print(beta, n_meth, n_unmeth, sep='\t', end='\n', file=fout)
                try:
                    idx_pos = index_file.__next__()
                except StopIteration:
                    break
"""
