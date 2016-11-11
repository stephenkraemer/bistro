# cython: profile=True
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
        'read_is_ok',
        'read_length',
        'ref',
        'sam_tlen',
        'watson_seq',
    ]

    def get_tlen(self):
        return max(abs(self.sam_tlen), self.read_length)

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
        self.watson_seq = self.pysam_aligned_segment.query_sequence
        self.read_length = len(self.watson_seq)
        self.direction = 'reverse' if self.pysam_aligned_segment.is_reverse else 'forward'
        self.mate = 1 if self.pysam_aligned_segment.is_read1 else 2
        self.bs_seq_strand = self.set_bs_seq_strand()

    def set_segment_attributes(self):
        self.name = self.pysam_aligned_segment.query_name
        self.next_ref_name = self.pysam_aligned_segment.next_reference_name
        self.next_ref_start = self.pysam_aligned_segment.next_reference_start
        self.sam_tlen = self.pysam_aligned_segment.template_length


class PileupSegment(metaclass=ABCMeta):

    __slots__ = [
        'aligned_segment',
        'baseq',
        'better_of_overlapping_pair',
        'bs_seq_strand',
        'has_frag_conversion_error',
        'is_ok',
        'meth_status_str',
        'overlap_status_dict',
        'omit_due_to_overlap',
        'omit_due_to_trimming_dict',
        'read_length',
        'same_mcall_from_mate',
        'snp_status_str',
        'trimming_site_array',
        'observed_watson_base',
        'watson_ref_base',
        'zero_based_pos_in_read',
        'zero_based_pos_in_watson_seq'
    ]

    def set_default_attribute_values(self):
        self.omit_due_to_overlap = False
        self.meth_status_str = ''
        self.snp_status_str = ''

    def add_methylation_status(self):
        mqc.methylation_calling.add_methylation_status_to_segment(self)

    """
    - trimming may be necessary for both adjusted and minimal mode, or only for minimal mode
    - trimming may be done using constant trimming sites (especially in minimal mode) or flen
      dependent sites (either per fragment length or in bins)
    - either two functions for only adjust or adjusted+minimal trimming, or delegation function
      if more options are needed
    """
    def add_trimming_sites_for_adjusted_mode(self, trimming_site_array):
        try:
            lowest_allowed_pos = trimming_site_array[0][self.aligned_segment.get_tlen()]
            highest_allowed_pos = trimming_site_array[1][self.aligned_segment.get_tlen()]
        except IndexError:  # fragment exceptionally long
            lowest_allowed_pos = trimming_site_array[0][-1]
            highest_allowed_pos = trimming_site_array[1][-1]

        self.omit_due_to_trimming_dict = dict()
        if lowest_allowed_pos <= self.zero_based_pos_in_read <= highest_allowed_pos:
            self.omit_due_to_trimming_dict['adjusted'] = False
        else:
            self.omit_due_to_trimming_dict['adjusted'] = True

    def __str__(self):
        if self.is_ok:
            return '''\
            Observed base: {}
            Reference base: {}
            BS-Seq strand: {}
            Methylation call: {}
            Name: {}
            Phred-Score'''.format(
                getattr(self, 'observed_watson_base', None),
                getattr(self, 'watson_ref_base', None),
                getattr(self.aligned_segment, 'bs_seq_strand', None),
                getattr(self, 'meth_status_str', None),
                getattr(self, 'aligned_segment.name', None),
                getattr(self, 'baseq', None))
        else:
            return 'Bad quality read'


class PysamPileupSegment(PileupSegment):
    def __init__(self, pysam_pileup_read: pysam.PileupRead, watson_ref_base,
                 mapq_threshold, phred_score_threshold, trimming_site_array,
                 frag_conv_err_detect_fun, index_position, max_number_of_unconv_control_cyts):

        self.is_ok = False

        self.aligned_segment = PysamAlignedSegment(pysam_pileup_read, mapq_threshold)
        if not self.aligned_segment.is_ok:
            return

        # pysam returns None as read position if the read has a
        # refskip or indel, and an integer otherwise
        self.zero_based_pos_in_watson_seq = pysam_pileup_read.query_position
        if not self.zero_based_pos_in_watson_seq:
            return

        self.observed_watson_base = self.aligned_segment.watson_seq[self.zero_based_pos_in_watson_seq]
        if self.observed_watson_base == 'N':
            return

        self.baseq = self.aligned_segment.baseqs_watson_seq[self.zero_based_pos_in_watson_seq]
        if (not self.zero_based_pos_in_watson_seq) or (self.baseq < phred_score_threshold):
            return

        if self.aligned_segment.direction == 'forward':
            self.zero_based_pos_in_read = self.zero_based_pos_in_watson_seq
        else:
            self.zero_based_pos_in_read = (
                len(self.aligned_segment.watson_seq) - self.zero_based_pos_in_watson_seq)

        self.has_frag_conversion_error = frag_conv_err_detect_fun(
            pileup_segment=self,
            index_position=index_position,
            max_number_of_unconv_control_cyts=max_number_of_unconv_control_cyts
        )
        if self.has_frag_conversion_error:
            return

        self.is_ok = True

        self.aligned_segment.set_segment_attributes()
        self.set_default_attribute_values()


        self.add_trimming_sites_for_adjusted_mode(trimming_site_array)
        self.watson_ref_base = watson_ref_base
        # TODO: replace by call to call_methylation, then manually add attribute to segment
        self.methylation_status = mqc.methylation_calling.add_methylation_status_to_segment(self)
