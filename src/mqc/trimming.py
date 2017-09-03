from mqc.mbias import CuttingSites
from mqc.visitors import Visitor
from mqc.pileup.pileup import MotifPileup
from mqc.flag_and_index_values import bsseq_strand_indices, methylation_status_flags
b_inds = bsseq_strand_indices
mflags = methylation_status_flags

class Trimmer(Visitor):
    """ Set trimming flag if position in read is outside of plateau positions

    Implementation notes:
    The cutting site array details zero-based positions (plateau_start, plateau_end)
    for every stratum (bsseq_strand, fragment length) with flen \in  [min_flen, max_flen]
    The (plateau_start, plateau_end) positions are the lowest/highest allowed position
    in the read.

    General thoughts:
    Only one flag value. The idea is to start with minimal trimming at the first pass over the
    random index positions. This will allow for the calculation of beta value dists etc. for minimal trimming. Then,
    one or more cutting site determination functions may be called. The result of every cutting site determination
    variant can then be determined in a second run. It is also conceivable to improve cutting sites iteratively,
    by passing over the random index positions multiple times until all QC metrics (e.g. beta value distributions)
    are passed
    """

    def __init__(self, config, cutting_sites: CuttingSites):
        self.cutting_sites = cutting_sites
        self.cutting_sites_array = cutting_sites.get_array()
        self.max_flen = config["trimming"]["max_flen_considered_for_trimming"]

    def process(self, motif_pileup: MotifPileup):
        for read in motif_pileup.reads:

            if (read.meth_status_flag == mflags.is_na
                or read.qc_fail_flag):
                continue

            obs_tlen = abs(read.alignment.template_length)
            tlen = obs_tlen if obs_tlen <= self.max_flen else self.max_flen

            bstrand_idx = read.bsseq_strand_ind
            start_of_plateau = self.cutting_sites_array[bstrand_idx, tlen, 0]
            end_of_plateau = self.cutting_sites_array[bstrand_idx, tlen, 1]

            if not (start_of_plateau <= read.pos_in_read <= end_of_plateau):
                # one could use different trimming modes and set different flag values for the trimming flag
                # currently, only the first bit is used
                read.trimm_flag = 1
