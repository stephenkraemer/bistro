import pysam
import mqc


def motif_pileup_generator(bam_path, index_file_path):
    # TODO: close SAM file?
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
            yield motif_pileups, curr_idx_pos

            try:
                curr_idx_pos = next(index_file)
            except StopIteration:
                return
