from mqc.mbias import CuttingSites
from mqc.visitors import Visitor
from mqc.pileup.pileup import MotifPileup
from mqc.flag_and_index_values import bsseq_strand_indices, methylation_status_flags
b_inds = bsseq_strand_indices
mflags = methylation_status_flags

class Trimmer(Visitor):
    """Set trimming flag if position in read is outside of plateau.

    Args:
        cutting_sites: start and end are slice boundaries according
            to the standard python slice notation. I.e. the plateau
            positions are all zero-based positions in [start, end[

    """

    def __init__(self, cutting_sites: CuttingSites) -> None:
        self.cutting_sites_array = cutting_sites.as_array()
        flen_idx = list(cutting_sites.df.index.names).index('flen')
        self.max_flen = self.cutting_sites_array.shape[flen_idx] - 1

    def process(self, motif_pileup: MotifPileup) -> None:
        """Add trimming flag to each read in MotifPileup.reads

        Trimming flags are given irrespective of other QC issues, to
        avoid biases in QC plots (e.g. because phred score failures
        would never be tagged as 'in trimming zone'.

        One could use different trimming modes and set different bits
        of the trimming flag. Currently, only the first bit
        is used.
        """
        for read in motif_pileup.reads:
            obs_tlen = abs(read.alignment.template_length)
            capped_tlen = obs_tlen if obs_tlen <= self.max_flen else self.max_flen

            bstrand_idx = read.bsseq_strand_ind
            start_of_plateau = self.cutting_sites_array[bstrand_idx, capped_tlen, 0]
            end_of_plateau = self.cutting_sites_array[bstrand_idx, capped_tlen, 1]

            if not start_of_plateau <= read.pos_in_read < end_of_plateau:
                read.trimm_flag = 1
