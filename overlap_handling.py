import mqc
import random
from typing import List

"""
Note: overlap tags are a manipulation of the segments in a MethPileup
Therefore, ideally the MethPileup should provide the interface for tagging its segments
However, depending on the task, different ways of tagging fragments will be desirable
I still need to figure out a fast way of repeatedly choosing the correct tagging function
in every instance of the MethPileups created during an analysis. Until then, I will call the
functions directly as module functions from the analysis workflow

If i call this from a MethPileup, the MethPileup could iterate over the pileup positions, instead
of having to do this in the functions below
"""


def add_conservative_tags_for_calling(meth_pileup: 'mqc.MethPileup'):
    for pileup_segments_at_pos in meth_pileup.pileup_segments_per_pos:
        overlapping_read_pairs = overlapping_read_pairs_generator(
            pileup_segments_at_pos, trimming_mode='adjusted')
        for reads in overlapping_read_pairs:
            if reads[0].meth_status_str == reads[1].meth_status_str:
                reads[random.randint(0, 1)].omit_due_to_overlap = True
            else:
                reads[0].omit_due_to_overlap = True
                reads[1].omit_due_to_overlap = True


def overlapping_read_pairs_generator(pileup_segments: 'List[mqc.PileupSegment]',
                                     trimming_mode='minimal'):
    pileup_segments_by_read_name_dict = dict()
    for curr_pileup_segment in pileup_segments:
        if (curr_pileup_segment.is_ok
                and not curr_pileup_segment.omit_due_to_trimming_dict[trimming_mode]):
            curr_read_name = curr_pileup_segment.aligned_segment.name
            if not pileup_segments_by_read_name_dict.get(curr_read_name, None):
                pileup_segments_by_read_name_dict[curr_read_name] = curr_pileup_segment
            else:
                yield [pileup_segments_by_read_name_dict[curr_read_name], curr_pileup_segment]

# class OverlapHandler:
#     def __init__(self, run_mode, overlap_handling_mode):
#         if run_mode == 'call' and overlap_handling_mode == 'discard_disagreeing_overlaps':
#             self.tag_overlapping_read_pairs = _tag_in_call_mode_conservatively
#         if run_mode == 'call' and overlap_handling_mode == 'better_call_from_overlap':
#             self.tag_overlapping_read_pairs = _tag_in_call_mode_take_better_read
#         else:
#             self.tag_overlapping_read_pairs = _tag_in_qc_mode
#
#
#
#
#
# def _tag_in_call_mode_conservatively(pileup_segments_at_pos):
#     trimming_mode = 'adjusted'
#     for read1, read2 in generate_overlapping_read_pairs(
#             pileup_segments_at_pos, trimming_mode):
#         tag_as_overlapping(read1, read2, trimming_mode)
#         add_methylation_comparison_tag(read1, read2)
#         """
#         tagging better read in pair is not necessary,
#         because pairs with unequal methylation calls will always be discarded, and otherwise
#         read1 from the pair will be used for methylation calling by default
#         """
#
#
# def _tag_in_call_mode_take_better_read(pileup_segments):
#     trimming_mode = 'adjusted'
#     for read1, read2 in generate_overlapping_read_pairs(
#             pileup_segments, trimming_mode):
#         tag_as_overlapping(read1, read2, trimming_mode)
#         add_methylation_comparison_tag(read1, read2)
#         # TODO: only tag better read if methylation status is not the same
#         tag_better_read_in_pair(read1, read2)
#
#
# def _tag_in_qc_mode(pileup_segments):
#     for read1, read2 in generate_overlapping_read_pairs(
#             pileup_segments, trimming_mode='minimal'):
#         tag_as_overlapping(read1, read2, trimming_mode='minimal')
#         add_methylation_comparison_tag(read1, read2)
#         tag_better_read_in_pair(read1, read2)
#         if not (read1.in_trimming_region['adjusted']
#                 or read2.in_trimming_region['adjusted']):
#             tag_as_overlapping(read1, read2, trimming_mode='adjusted')
#
#
#
#
# def tag_as_overlapping(read1: 'mqc.PileupSegment', read2: 'mqc.PileupSegment', trimming_mode: str):
#     read1.overlap_status_dict[trimming_mode] = True
#     read2.overlap_status_dict[trimming_mode] = True
#
#
# def add_methylation_comparison_tag(read1: 'mqc.PileupSegment', read2: 'mqc.PileupSegment'):
#     if read1.meth_status_str == read2.meth_status_str:
#         read1.same_mcall_from_mate = True
#         read2.same_mcall_from_mate = True
#
#
# def tag_better_read_in_pair(read1: 'mqc.PileupSegment', read2: 'mqc.PileupSegment'):
#     if read1.baseq > read2.baseq:
#         read1.better_of_overlapping_pair = True
#         read2.better_of_overlapping_pair = False
#     elif read2.baseq > read1.baseq:
#         read1.better_of_overlapping_pair = False
#         read2.better_of_overlapping_pair = True
#     elif read1.aligned_segment.mapq > read2.aligned_segment.mapq:
#         read1.better_of_overlapping_pair = True
#         read2.better_of_overlapping_pair = False
#     elif read2.aligned_segment.mapq > read1.aligned_segment.mapq:
#         read1.better_of_overlapping_pair = False
#         read2.better_of_overlapping_pair = True
#     else:
#         # everything the same
#         pass
#
#
# """
# Notes
#
# def tag_overlapping_reads_in_call_mode(read1: mqc.PileupSegment, read2: mqc.PileupSegment):
#     read1.
#     read1.overlaps_after_minimal_trimming = True
#     read2.overlaps_after_minimal_trimming = True
#     if read1.methylation_status != read2.methylation_status:
#         read1._overlapping_status = 'discard'
#         read2._overlapping_status = 'discard'
#     else:
# # to the methylation statuses agree?
# # which read has better quality at that position?
# # has this read an overlapping mate?
# # for both trimmed and untrimmed situation
#
# check_for_overlap_after_called_trimming()
#
# # discard from overlap comparison if baseq is below threshold
# """
