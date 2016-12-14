import numpy as np
import mqc
from typing import List

b_flags = mqc.flag_and_index_values.bsseq_strand_flags
b_inds = mqc.flag_and_index_values.bsseq_strand_indices


def set_trimming_flag(pileup_reads: 'List[mqc.bsseq_pileup_read.BSSeqPileupRead]', trimm_site_array):
    """ Only one flag value. The idea is to start with minimal trimming at the first pass over the
    random index positions. This will allow for the calculation of beta value dists etc. for minimal trimming. Then,
    one or more cutting site determination functions may be called. The result of every cutting site determination
    variant can then be determined in a second run. It is also conceivable to improve cutting sites iteratively,
    by passing over the random index positions multiple times until all QC metrics (e.g. beta value distributions)
    are passed"""
    # TODO: avoid hardcoding
    max_tlen = 399
    for read in pileup_reads:
        # TODO: make sure that ISNA flag is set for methylation status of indel, refskip and base N reads
        if not read.query_position:
            continue
        obs_tlen = abs(read.alignment.template_length)
        # TODO: tlen should not be smaller than length of query sequence
        tlen = obs_tlen if obs_tlen <= max_tlen else max_tlen
        if read.bs_seq_strand_flag == b_flags.c_bc:
            start_of_plateau = trimm_site_array[b_inds.c_bc, 0, tlen]
            end_of_plateau = trimm_site_array[b_inds.c_bc, 1, tlen]
        elif read.bs_seq_strand_flag == b_flags.c_bc_rv:
            start_of_plateau = trimm_site_array[b_inds.c_bc_rv, 0, tlen]
            end_of_plateau = trimm_site_array[b_inds.c_bc_rv, 1, tlen]
        elif read.bs_seq_strand_flag == b_flags.w_bc:
            start_of_plateau = trimm_site_array[b_inds.w_bc, 0, tlen]
            end_of_plateau = trimm_site_array[b_inds.w_bc, 1, tlen]
        elif read.bs_seq_strand_flag == b_flags.w_bc_rv:
            start_of_plateau = trimm_site_array[b_inds.w_bc_rv, 0, tlen]
            end_of_plateau = trimm_site_array[b_inds.w_bc_rv, 1, tlen]
        else:  # bisulfite sequencing strand undetermined, no methylation calling possible
            # TODO: make sure that ISNA flag is set for methylation status of indel, refskip and base N reads
            continue

        if not (start_of_plateau <= read.query_position <= end_of_plateau):
            # one could use different trimming modes and set different flag values for the trimming flag
            # currently, only the first bit is used
            read.trimm_flag = 1


def cutting_sites_array_from_flen_relative_minimal_cutting_sites(relative_cutting_site_dict,
                                                                 max_flen_considered_for_trimming,
                                                                 max_read_length_bp):
    """
    Relative cutting site dict:
    {
        'w_bc': [0, 9]
        'c_bc': [0, 9]
        'w_bc_rv': [9, 0]
        'c_bc_rv': [9, 0]
    }
    """

    """ Create res array

    bsseq_strand      start_or_end       flen
    'w_bc'            0                  1
                                         2
                                         ...
                                         max_flen_considered_for_trimming
                      1                  1
                                         ...
    'w_bc_rv'         ...                ...

    """
    res = np.zeros([4, 2, max_flen_considered_for_trimming + 1], dtype=np.int32)

    for bsseq_strand_name, bsseq_strand_index in b_inds._asdict().items():
        # Set start position
        res[bsseq_strand_index, 0, :] = relative_cutting_site_dict[bsseq_strand_name][0]

        # Set end position
        for curr_flen in range(max_flen_considered_for_trimming + 1):
            max_allowed_pos_in_fragment = curr_flen - relative_cutting_site_dict[bsseq_strand_name][1]
            max_allow_pos_in_read = (max_allowed_pos_in_fragment
                                     if max_allowed_pos_in_fragment <= max_read_length_bp
                                     else max_read_length_bp)
            res[bsseq_strand_index, 1, curr_flen] = max_allow_pos_in_read

    return res
