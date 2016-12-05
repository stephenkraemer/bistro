""" Parallelized methylation calling"""
import gzip
import queue
import pysam
import mqc
import multiprocessing as mp
import mqc.bsseq_pileup_read
import mqc.methylation_calling.pileup_methylation_calling as mc

LAST_CHROM = 0
LAST_START = 0


def call_methylation_parallelized(bam_file_path, index_file_path,
                                  n_cores, output_file, pos_per_process, lock):
    # TODO: pool.apply may be better suited
    with gzip.open(output_file, 'wt') as fout:

        results_queue = mp.Queue()  # returns process results in sorted manner
        process_queue = queue.Queue(maxsize=n_cores)  # FIFO for getting processes in order
        """
        Split index into chunks of equal size, then work through chunks in /n_cores/
        parallel processes, writing them to file in a sorted manner
        """
        process_kwargs = generate_process_args(index_file_path=index_file_path,
                                               pos_per_process=pos_per_process,
                                               bam_path=bam_file_path,
                                               results_queue=results_queue,
                                               lock=lock)

        for i in range(n_cores):
            # Does the number of processes have to be the same as the number of cores?
            kwargs = next(process_kwargs)
            start_new_process(process_queue, kwargs)

        while 1:
            # TODO: blocking probably not necessary
            next_process = process_queue.get(block=True, timeout=None)
            wrap_up_process(next_process, results_queue=results_queue, fout=fout)
            try:
                kwargs = next(process_kwargs)
                start_new_process(process_queue, kwargs)
            except StopIteration:
                break

        # Wait for remaining processes in queue (n_cores - 1 processes left, because
        # we could not replace one item in the while loop above
        for i in range(n_cores - 1):
            next_process = process_queue.get(block=True, timeout=None)
            wrap_up_process(next_process, results_queue=results_queue, fout=fout)


def start_new_process(process_queue, kwargs):
    p = mp.Process(target=run_on_chunk, kwargs=kwargs)
    p.start()
    process_queue.put(p, block=True, timeout=None)


def wrap_up_process(proc, results_queue, fout):
    global LAST_START, LAST_CHROM
    proc.join()
    # TODO Could there be a delay in addition of the process result to the queue? If the
    # process result is added with a delay and a lower priority results is already
    # present, it may be returned instead
    start, chrom, lines = results_queue.get(block=True, timeout=15)
    """
    if chrom < LAST_CHROM or start < LAST_START:
        raise IndexError
    else:
        LAST_CHROM = chrom
        LAST_START = start
    """
    fout.write(lines)


def generate_process_args(index_file_path, pos_per_process, bam_path, results_queue, lock):
    # TODO: stop at chromosome transitions
    # TODO: deal with multiple index files (for CH context)
    # TODO: We may have to skipt the header of the index file
    with gzip.open(index_file_path, 'rt') as file_path_obj:
        still_lines_in_file = True
        while still_lines_in_file:
            index_positions = []
            for i in range(pos_per_process):
                try:
                    index_positions.append(mqc.IndexPosition(next(file_path_obj)))
                except StopIteration:
                    still_lines_in_file = False
                    if not index_positions:
                        return
                    else:
                        break
            kwargs = {
                'index_positions': index_positions,
                'chrom': index_positions[0].chrom,
                'start': index_positions[0].start,
                'end': index_positions[-1].end,
                'bam_path': bam_path,
                'results_queue': results_queue,
                'lock': lock,
            }
            yield kwargs


def run_on_chunk(index_positions, chrom, start, end, bam_path, results_queue, lock):
    sam_file = pysam.AlignmentFile(bam_path)
    pileup_columns = sam_file.pileup(chrom, start, end, truncate=True)
    lines = []
    curr_pileup_pos = start - 1
    index_position_iter = iter(index_positions)
    curr_idx_pos = next(index_position_iter)
    for pileup_column in pileup_columns:
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

            # print(beta, n_meth, n_unmeth, sep='\t', end='\n', file=fout)
            lines.append('\t'.join([str(beta), str(n_meth), str(n_unmeth)]))
            try:
                curr_idx_pos = next(index_position_iter)
            except StopIteration:
                break
    # TODO: convert chrom to integer, deal with X and Y (23, 24?)
    lines_str = '\n'.join(lines)
    results_queue.put((int(curr_idx_pos.chrom), curr_idx_pos.start, lines_str))


if __name__ == '__main__':
    lock = mp.Lock()
    call_methylation_parallelized(
        bam_file_path='./test/data/b_cells_rep1_chr11_16815793-16824254.bam',
        index_file_path='./test/data/chr11_16815793-16824254.cg.bed.gz',
        n_cores=2,
        output_file='./test/results/methylation_calls_new.gz',
        pos_per_process=3,
        lock=lock)
