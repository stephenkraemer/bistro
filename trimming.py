"""
Given a numpy array, iterate over a pileup of reads and set trimming flags
"""

import numpy as np

W_BC = 96
C_BC = 80
W_BC_RV = 144
C_BC_RV = 160

C_BC_IND = 0
C_BC_RV_IND = 2
W_BC_IND = 4
W_BC_RV_IND = 6


def cutting_sites_from_mbias_stats(mbias_stats_array):
    # Indexing will be done by integer TLEN - take index 0 as nonsense element into account
    res = np.ones([8, 400], dtype=np.int32)
    res[[0, 2, 4, 6], :] = 10
    res[[1, 3, 5, 7], :] = 90
    return res


def set_trimming_flag(pileup_reads: 'List[pysam.PileupRead]', trimm_site_array):
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
        tlen = obs_tlen if obs_tlen <= max_tlen else max_tlen
        if read.bs_seq_strand_flag == C_BC:
            start_of_plateau = trimm_site_array[C_BC_IND, tlen]
            end_of_plateau = trimm_site_array[C_BC_IND + 1, tlen]
        elif read.bs_seq_strand_flag == C_BC_RV:
            start_of_plateau = trimm_site_array[C_BC_RV_IND, tlen]
            end_of_plateau = trimm_site_array[C_BC_RV_IND + 1, tlen]
        elif read.bs_seq_strand_flag == W_BC:
            start_of_plateau = trimm_site_array[W_BC_IND, tlen]
            end_of_plateau = trimm_site_array[W_BC_IND + 1, tlen]
        elif read.bs_seq_strand_flag == W_BC_RV:
            start_of_plateau = trimm_site_array[W_BC_RV_IND, tlen]
            end_of_plateau = trimm_site_array[W_BC_RV_IND + 1, tlen]
        else:  # bisulfite sequencing strand undetermined, no methylation calling possible
            # TODO: make sure that ISNA flag is set for methylation status of indel, refskip and base N reads
            continue

        if not (start_of_plateau <= read.query_position <= end_of_plateau):
            # one could use different trimming modes and set different flag values for the trimming flag
            # currently, only the first bit is used
            read.trimm_flag = 1
