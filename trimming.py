import mqc
from typing import List

b_inds = mqc.flag_and_index_values.bsseq_strand_indices
b_na_ind = mqc.flag_and_index_values.bsseq_strand_na_index


def set_trimming_flag(pileup_reads: 'List[mqc.bsseq_pileup_read.BSSeqPileupRead]',
                      trimm_site_array, max_flen_considered_for_trimming):
    """ Only one flag value. The idea is to start with minimal trimming at the first pass over the
    random index positions. This will allow for the calculation of beta value dists etc. for minimal trimming. Then,
    one or more cutting site determination functions may be called. The result of every cutting site determination
    variant can then be determined in a second run. It is also conceivable to improve cutting sites iteratively,
    by passing over the random index positions multiple times until all QC metrics (e.g. beta value distributions)
    are passed"""
    for read in pileup_reads:
        # TODO: make sure that ISNA flag is set for methylation status of indel, refskip and base N reads
        if not read.query_position:
            continue
        obs_tlen = abs(read.alignment.template_length)
        # TODO: tlen should not be smaller than length of query sequence
        tlen = obs_tlen if obs_tlen <= max_flen_considered_for_trimming else max_flen_considered_for_trimming

        bsseq_strand_ind = read.bsseq_strand_ind
        # TODO: make sure that ISNA flag is set for methylation status of indel, refskip and base N reads
        if bsseq_strand_ind == b_na_ind:
            continue
        start_of_plateau = trimm_site_array[bsseq_strand_ind, 0, tlen]
        end_of_plateau = trimm_site_array[bsseq_strand_ind, 1, tlen]

        if not (start_of_plateau <= read.query_position <= end_of_plateau):
            # one could use different trimming modes and set different flag values for the trimming flag
            # currently, only the first bit is used
            read.trimm_flag = 1
