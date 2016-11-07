import mqc

from abc import ABCMeta, abstractmethod

import pysam


class AlignedSegment(metaclass=ABCMeta):
    bs_strand_dict_directional_protocol = {1: {'forward': 'W-BC',
                                               'reverse': 'C-BC'},
                                           2: {'forward': 'C-BC-Rv',
                                               'reverse': 'W-BC-Rv'}}

    __slots__ = [
        'baseqs_watson_seq',
        'bs_seq_strand',
        'direction',
        'flag',
        'mapq',
        'mate',
        'name',
        'next_ref_name',
        'next_ref_start',
        'pos',
        'read_is_ok'
        'ref',
        'sam_tlen',
        'watson_seq',
    ]

    def get_tlen(self):
        return max(abs(self.sam_tlen), len(self.watson_seq))

    def set_read_is_ok(self, mapq_threshold):
        """ Basic quality test for segment

        Implementation notes:
        The flag value is compared with the flag value 3855, which has the binary flags for
        - read paired
        - read mapped in proper pair
        - read unmapped
        - mate unmapped
        - not primary alignment
        - read fails QC check
        - read is duplicate
        - read is supplementary

        If the flag value only has the first two of these bits set, the bitwise AND comparison will
        return 3

        Otherwise, it will return a number larger than 3, since another bit will survive the
        AND comparison
        """
        if self.mapq < mapq_threshold:
            return False
        if not self.flag & 3855 == 3:
            return False
        return True

    def set_bs_seq_strand(self):
        return self.bs_strand_dict_directional_protocol[self.mate][self.direction]

    @abstractmethod
    def set_segment_attributes(self):
        pass


class PysamAlignedSegment(AlignedSegment):
    def __init__(self, pysam_pileup_read, mapq_threshold):
        """
        Note:
            - this function assumes that reads aligning to crick were stored
            reverse-complemented by the aligner
        """
        self.pysam_aligned_segment = pysam_pileup_read.alignment
        self.mapq = self.pysam_aligned_segment.mapping_quality
        self.flag = self.pysam_aligned_segment.flag
        self.is_ok = self.set_read_is_ok(mapq_threshold)
        self.baseqs_watson_seq = self.pysam_aligned_segment.query_qualities

    def set_segment_attributes(self):
        self.mate = 1 if self.pysam_aligned_segment.is_read1 else 2
        self.direction = 'reverse' if self.pysam_aligned_segment.is_reverse else 'forward'
        self.bs_seq_strand = self.set_bs_seq_strand()
        self.name = self.pysam_aligned_segment.query_name
        self.next_ref_name = self.pysam_aligned_segment.next_reference_name
        self.next_ref_start = self.pysam_aligned_segment.next_reference_start
        self.watson_seq = self.pysam_aligned_segment.query_sequence


class PileupSegment(metaclass=ABCMeta):

    __slots__ = [
        'aligned_segment',
        'baseq',
        'better_of_overlapping_pair',
        'bs_seq_strand',
        'has_base_at_pos',
        'has_del_at_pos',
        'has_refskip_at_pos',
        'is_ok',
        'meth_status_str',
        'overlap_status_dict',
        'trimming_mode_status_dict',
        'same_mcall_from_mate',
        'snp_status_str',
        'watson_base',
        'watson_ref_base',
        'zero_based_pos_in_read'
    ]

    def set_default_attribute_values(self):
        self.trimming_mode_status_dict = {'minimal': False,
                                          'adjusted': False}
        self.overlap_status_dict = {'minimal': False,
                                    'adjusted': True}
        self.same_mcall_from_mate = False
        self.meth_status_str = ''
        self.snp_status_str = ''

    def add_methylation_status(self):
        mqc.methylation_calling.add_methylation_status(self)

    def __str__(self):
        if self.is_ok:
            return '''\
            Observed base: {}
            Reference base: {}
            BS-Seq strand: {}
            Methylation call: {}
            Name: {}
            Phred-Score'''.format(
                getattr(self, 'watson_base', None),
                getattr(self, 'watson_ref_base', None),
                getattr(self.aligned_segment, 'bs_seq_strand', None),
                getattr(self, 'meth_status_str', None),
                getattr(self, 'aligned_segment.name', None),
                getattr(self, 'baseq', None))
        else:
            return 'Bad quality read'


class PysamPileupSegment(PileupSegment):
    @profile
    def __init__(self, pysam_pileup_read: pysam.PileupRead, watson_ref_base,
                 mapq_threshold, phred_score_threshold):

        self.is_ok = False

        self.aligned_segment = PysamAlignedSegment(pysam_pileup_read, mapq_threshold)
        if not self.aligned_segment.is_ok:
            return

        # pysam returns None as read position if the read has a
        # refskip or indel, and an integer otherwise
        self.zero_based_pos_in_watson_seq = pysam_pileup_read.query_position

        self.baseq = self.aligned_segment.baseqs_watson_seq[self.zero_based_pos_in_watson_seq]
        if (not self.zero_based_pos_in_watson_seq) or (self.baseq < phred_score_threshold):
            return

        self.is_ok = True

        self.aligned_segment.set_segment_attributes()

        self.set_default_attribute_values()

        self.watson_base = self.aligned_segment.watson_seq[self.zero_based_pos_in_watson_seq]
        self.watson_ref_base = watson_ref_base

        self.methylation_status = mqc.methylation_calling.add_methylation_status(self)
