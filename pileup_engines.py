# TODO: is this a good structure for threading/multiprocessing/threading+multiprocessing combined?
import pysam
import mqc

from abc import ABCMeta, abstractmethod


class AbstractPileupEngine(metaclass=ABCMeta):
    """Abstract iterator over MethInfoAtPos instances"""

    @abstractmethod
    def __iter__(self) -> mqc.MethPileup:
        pass


class PysamPileupEngine(AbstractPileupEngine):
    """Iteratively return MethInfoAtPos via pysam.pileup function"""

    def __init__(self, index_file: mqc.IndexFile, bam_abspath: str,
                 mapq_threshold, phred_score_threshold,
                 run_mode, overlap_handling_mode):
        self.index_file = index_file
        self.pysam_alignment_file = pysam.AlignmentFile(bam_abspath, 'rb')
        self.mapq_threshold = mapq_threshold
        self.phred_score_threshold = phred_score_threshold
        self.run_mode = run_mode,
        self.overlap_handling_mode = overlap_handling_mode

    def __next__(self) -> mqc.MethPileup:
        index_position = next(self.index_file)
        pysam_pileup_segments_per_pos = self._get_list_of_pileup_segments_per_pos(index_position)
        meth_pileup = mqc.MethPileup(pysam_pileup_segments_per_pos,
                                     index_position,
                                     run_mode=self.run_mode,
                                     overlap_handling_mode=self.overlap_handling_mode)
        return meth_pileup

    def _get_list_of_pileup_segments_per_pos(self, index_position: mqc.index.IndexPosition):
        pileup_segments_per_pos = []
        pysam_pileup = self.pysam_alignment_file.pileup(reference=index_position.chrom,
                                                        start=index_position.start,
                                                        end=index_position.end,
                                                        truncate=True)
        for watson_ref_base, pysam_pileup_column in zip(index_position.watson_motif, pysam_pileup):
            if watson_ref_base in ['C', 'G']:
                pysam_pileup_reads = pysam_pileup_column.pileups
                pileup_segments = [mqc.PysamPileupSegment(
                                        x, watson_ref_base,
                                        mapq_threshold=self.mapq_threshold,
                                        phred_score_threshold=self.phred_score_threshold)
                                   for x in pysam_pileup_reads]
                pileup_segments_per_pos.append(pileup_segments)
            else:
                pileup_segments_per_pos.append([])

        return pileup_segments_per_pos

    def __iter__(self) -> mqc.MethPileup:
        return self
