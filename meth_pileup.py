import itertools
import mqc
from typing import List


class MethPileup:
    def __init__(self,
                 pileup_segments_per_pos: List[mqc.PileupSegment],
                 index_position: mqc.IndexPosition,
                 run_mode,
                 overlap_handling_mode):
        self.pileup_segments_per_pos = pileup_segments_per_pos
        self.index_position = index_position
        self._snp_score = None
        self._overlap_handler = mqc.OverlapHandler(run_mode=run_mode,
                                                   overlap_handling_mode=overlap_handling_mode)

    @property
    def snp_score(self):
        if self._snp_score:
            return self._snp_score
        else:
            self._snp_score = mqc.get_snp_score(self)
            return self._snp_score

    def tag_overlapping_reads(self):
        for pileup_segments in self.pileup_segments_per_pos:
            self._overlap_handler.tag_overlapping_read_pairs(pileup_segments)

    def get_total_meth_stats(self):
        """
        for watson_base, pileup_segments_at_pos in zip(
            meth_pileup.index_position.watson_motif,
            meth_pileup.pileup_segments_per_pos):
            if
        """
        n_meth, n_unmeth = 0, 0
        for pileup_segment in itertools.chain.from_iterable(self.pileup_segments_per_pos):
            if not pileup_segment.is_ok:
                continue

            if pileup_segment.meth_status_str == 'methylated':
                n_meth += 1
            elif pileup_segment.meth_status_str == 'unmethylated':
                n_unmeth += 1

        n_total = n_meth + n_unmeth
        try:
            beta = n_meth / n_total
        except ZeroDivisionError:
            beta = None

        return beta, n_meth, n_total

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
