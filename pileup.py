import pysam
import mqc


def motif_pileup_generator(bam_path, index_file_path):

    index_file = mqc.IndexFile(index_file_path)
    curr_idx_pos = next(index_file)
    pileup_columns = all_position_pileup_generator(bam_path=bam_path,
                                                   ref=curr_idx_pos.chrom, start=curr_idx_pos.start)

    assumed_pos = curr_idx_pos.start - 1
    for pileup_column, curr_pileup_pos in pileup_columns:
        assumed_pos += 1
        if pileup_column and assumed_pos != pileup_column.reference_pos:
            raise(ValueError('Positions dont mach'))

        if curr_pileup_pos == curr_idx_pos.start:
            curr_base = curr_idx_pos.watson_motif[0]
            if curr_base in ['C', 'G'] and pileup_column:
                motif_pileups = [mqc.bsseq_pileup_read.pileups(pileup_column)]
            else:
                motif_pileups = [[]]
            for curr_base in curr_idx_pos.watson_motif[1:]:
                pileup_column, curr_pileup_pos = next(pileup_columns)
                assumed_pos += 1
                if pileup_column and assumed_pos != pileup_column.reference_pos:
                    raise(ValueError('Positions dont mach'))
                if curr_base in ['C', 'G'] and pileup_column:
                    motif_pileups.append(mqc.bsseq_pileup_read.pileups(pileup_column))
                else:
                    motif_pileups.append([])
            yield motif_pileups, curr_idx_pos

            try:
                curr_idx_pos = next(index_file)
            except StopIteration:
                return

    # Fill remaining index positions with zeros
    motif_pileups = [[] for base in curr_idx_pos.watson_motif]
    yield motif_pileups

    for curr_idx_pos in index_file:
        motif_pileups = [[] for base in curr_idx_pos.watson_motif]
        yield motif_pileups


def all_position_pileup_generator(bam_path, ref, start):
    # TODO: close SAM file?
    sam_file = pysam.AlignmentFile(bam_path)
    pileup_columns = sam_file.pileup(reference=ref,
                                     start=start,
                                     end=None,
                                     truncate=True)

    last_pos = start - 1
    for pileup_column in pileup_columns:
        curr_pos = pileup_column.reference_pos
        while last_pos + 1 < curr_pos:
            yield (None, last_pos + 1)
            last_pos += 1
        last_pos = curr_pos
        yield (pileup_column, curr_pos)
