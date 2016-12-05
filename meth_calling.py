import mqc
import pysam
import gzip
import mqc.methylation_calling.pileup_methylation_calling as mc


def tabulate_meth_calls(bam_path, index_file_path, output_file):
    """
    Take an index file
    iterate over index file
    write output to file
    """
    with gzip.open(output_file, 'wt') as fout:
        sam_file = pysam.AlignmentFile(bam_path)

        index_file = mqc.IndexFile(index_file_path)
        curr_idx_pos = next(index_file)
        global_start = curr_idx_pos.start
        pileup_columns = sam_file.pileup(chrom=curr_idx_pos.chrom,
                                         start=curr_idx_pos.start,
                                         end=None,
                                         truncate=True)


        curr_pileup_pos = curr_idx_pos.start - 1
        arrived_at_start = False
        for pileup_column in pileup_columns:
            if pileup_column.reference_pos < global_start:
                continue
            if not arrived_at_start:
                print('arrived at start')
                arrived_at_start = True
            curr_pileup_pos += 1
            if curr_pileup_pos == curr_idx_pos.start:
                curr_base = curr_idx_pos.watson_motif[0]
                if curr_base in ['C', 'G']:
                    motif_pileups = [mqc.bsseq_pileup_read.pileups(pileup_column)]
                else:
                    motif_pileups = [[]]
                for curr_base in curr_idx_pos.watson_motif[1:]:
                    pileup_column = next(pileup_columns)
                    curr_pileup_pos += 1
                    if curr_base in ['C', 'G']:
                        motif_pileups.append(mqc.bsseq_pileup_read.pileups(pileup_column))
                    else:
                        motif_pileups.append([])

                beta, n_meth, n_unmeth = mc.call_meth_at_pileup(
                    motif_pileups,
                    index_position=curr_idx_pos)

                print(curr_idx_pos.chrom,
                      curr_idx_pos.start,
                      curr_idx_pos.end,
                      beta,
                      n_meth,
                      n_unmeth,
                      sep='\t', end='\n', file=fout)

                try:
                    curr_idx_pos = next(index_file)
                except StopIteration:
                    return

if __name__ == '__main__':
    """
    export PYTHONPATH=$PYTHONPATH:/home/kraemers/projects/mqc
    """
    import time
    print('Working')
    t0 = time.time()
    # tabulate_meth_calls(bam_path='/home/kraemers/projects/mqc/mqc/test/data/b_cells_rep1_chr11_16815793-16824254.bam',
    #                     index_file_path='./test/data/chr11_16815793-16824254.cg.bed.gz',
    #                     output_file='./test/results/methylation_calls_new.gz')
    tabulate_meth_calls(bam_path='/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/results_per_pid/hsc_rep2/alignment/blood_hsc_rep2_merged.mdup.bam',
                        index_file_path='/home/kraemers/projects/mqc/mqc/test/data/chr11.100000.cg.bed.gz',
                        output_file='./test/results/methylation_calls_new.gz')
    t1 = time.time()
    total = t1-t0
    print('Time for 100,000 positions(s): ', total)
    print('Projected time(h): ', total/100000 * 30 * 10**6 / 60 / 60)
