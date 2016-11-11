import mqc
from typing import List


class MethPileup:
    def __init__(self,
                 pileup_segments_per_pos: 'List[mqc.PileupSegment]',
                 index_position: mqc.IndexPosition):
        self.pileup_segments_per_pos = pileup_segments_per_pos
        self.index_position = index_position
        self._snp_score = None

    @property
    def snp_score(self):
        if self._snp_score:
            return self._snp_score
        else:
            self._snp_score = mqc.get_snp_score(self)
            return self._snp_score

    def tag_overlapping_reads(self):
        """See overlap handling module for more info, currently called as module function
        from analysis workflows"""
        return NotImplementedError

    def get_total_meth_stats(self, trimming_mode):
        """
        for observed_watson_base, pileup_segments_at_pos in zip(
            meth_pileup.index_position.watson_motif,
            meth_pileup.pileup_segments_per_pos):
            if
        """
        beta, n_meth, n_unmeth = mqc.methylation_calling.get_total_meth_stats_at_pileup(
            self, trimming_mode)
        return beta, n_meth, n_unmeth

    def __str__(self):
        lines = list()
        lines.append('MethPileup at {}'.format(self.index_position))
        for base, pileup_segments_at_base in zip(
                self.index_position.watson_motif, self.pileup_segments_per_pos):
            n_reads_at_base = len(pileup_segments_at_base)
            motif_base_on_strand = (base if self.index_position.strand == '+'
                                    else mqc.utilities.complement_seq(base))

            lines.append('{} reads at motif base {}'.format(
                n_reads_at_base,
                motif_base_on_strand
            ))

            lines.append('First read:')
            lines.append(str(pileup_segments_at_base[0]))
        return '\n'.join(lines)
