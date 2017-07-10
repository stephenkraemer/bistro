from collections import defaultdict
from typing import List
import mqc
import random
random.seed(1234)

def tag_overlaps(pileups: 'List[mqc.pileup.bsseq_pileup_read.BSSeqPileupRead]'):
    read_hash = defaultdict(list)
    for pileupread in pileups:
        if not (pileupread.trimm_flag or pileupread.qc_fail_flag):
            read_hash[pileupread.alignment.query_name].append(pileupread)
    for query_name, reads in read_hash.items():
        if len(reads) == 2:
            process_overlap(reads[0], reads[1])


def process_overlap(read1: 'mqc.BSSeqPileupRead', read2: 'mqc.BSSeqPileupRead'):
    if read1.observed_watson_base == read2.observed_watson_base:
        # both mates agree, but have to pick one event for usage in mcalling
        choice = random.randint(0,1)
        if choice == 0:
            read1.overlap_flag = 1
        else:
            read2.overlap_flag = 1
        """ Note on choosing the mate to tag:
        for mate- or strand-stratified statistics, the overlap flag will
        generally be ignored, because we want to count the events on both
        strands if both events are agreeing and have good quality. Randomly
        setting the overlap flag is just intended to make the flag more
        robust against read biases when the user forgets to use both events
        in stratified statistics
        """
    else:
        # Better safe than sorry. Statistically valid usage of phred
        # scores to select high confidence calls from disagreeing overlaps
        # will be implemented in the future
        read1.qc_fail_flag = 1
        read2.qc_fail_flag = 1
        """
        if read1.baseq_at_pos > read2.baseq_at_pos:
            read2.overlap_flag = 1
        elif read2.baseq_at_pos > read1.baseq_at_pos:
            read1.overlap_flag = 1
        else:  # different calls, same base qualities -> discard
            read1.qc_fail_flag = 1
            read2.qc_fail_flag = 1
        """

# def tag_overlaps(pileups: 'List[mqc.bsseq_pileup_read.BSSeqPileupRead]'):
#     read_hash = defaultdict(list)
#     for pileupread in pileups:
#         if not pileupread.trimm_flag:
#             read_hash[pileupread.alignment.query_name].append(pileupread)
#     for query_name, reads in read_hash.items():
#         if len(reads) == 2:
#             process_overlap(reads[0], reads[1])
#
#
# def process_overlap(read1: 'mqc.BSSeqPileupRead', read2: 'mqc.BSSeqPileupRead'):
#     if read1.observed_watson_base == reverse_complement_seq(read2.observed_watson_base):
#         read2.overlap_flag = 1
#     else:
#         if read1.baseq_at_pos > read2.baseq_at_pos:
#             read2.overlap_flag = 1
#         elif read2.baseq_at_pos > read1.baseq_at_pos:
#             read1.overlap_flag = 1
#         else:  # different calls, same base qualities -> discard
#             read1.qc_fail_flag = 1
#             read2.qc_fail_flag = 1
#
#
# def reverse_complement_seq(seq):
#     return seq.translate(str.maketrans('CGHD', 'GCDH'))[::-1]
