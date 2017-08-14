import random
random.seed(1234)

from collections import defaultdict

from mqc.visitors import Visitor
from mqc.pileup.pileup import MotifPileup
from mqc.pileup.bsseq_pileup_read import BSSeqPileupRead
from mqc.flag_and_index_values import qc_fail_flags
qflags = qc_fail_flags


class OverlapHandler(Visitor):
    def __init__(self):
        pass

    def process(self, motif_pileup: MotifPileup):
        read_hash = defaultdict(list)
        for curr_read in motif_pileup.reads:
            if not (curr_read.trimm_flag or curr_read.qc_fail_flag):
                read_hash[curr_read.alignment.query_name].append(curr_read)
        for query_name, reads in read_hash.items():
            if len(reads) == 2:
                process_overlap(reads[0], reads[1])


def process_overlap(read1: BSSeqPileupRead, read2: BSSeqPileupRead):
    if read1.meth_status_flag == read2.meth_status_flag:
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
        read1.qc_fail_flag |= qflags.overlap_fail
        read2.qc_fail_flag |= qflags.overlap_fail
