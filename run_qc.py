import pysam
import mqc
import mqc.methylation_calling.pileup_methylation_calling as mc
import numpy as np
"""
- run over index chunk
- calculate methylation status
- add coverage to predefined numpy array
- print array
"""


def get_cpg_coverage(index_positions, bam_path):
    sam_file = pysam.AlignmentFile(bam_path)

    chrom = index_positions[0].chrom
    start = index_positions[0].start
    # TODO: is the end coordinate included? Necessary to catch last motif
    end = index_positions[-1].end
    index_position_iter = iter(index_positions)

    pileup_columns = sam_file.pileup(chrom, start, end, truncate=True)

    curr_pileup_pos = start - 1
    curr_idx_pos = next(index_position_iter)

    cov_hist = np.zeros(1, 100)

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

            cov = n_meth + n_unmeth

            cov_hist[cov] += 1

            try:
                curr_idx_pos = next(index_position_iter)
            except StopIteration:
                return cov_hist


