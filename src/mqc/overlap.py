"""Overlap handling"""

import random
from typing import DefaultDict, List

random.seed(1234)

from collections import defaultdict

from mqc.visitors import Visitor
from mqc.pileup.pileup import MotifPileup
from mqc.pileup.bsseq_pileup_read import BSSeqPileupRead
from mqc.flag_and_index_values import qc_fail_flags, methylation_status_flags
qflags = qc_fail_flags
mflags = methylation_status_flags


class OverlapHandler(Visitor):
    """Find overlapping reads, evaluate, pick one event or discard both"""

    def process(self, motif_pileup: MotifPileup) -> None:
        """Iterate over reads in MotifPileup to find and handle overlaps

        An overlap event is not called when one of the reads
        - is trimmed at the position
        - has a qc_fail_flag
        - has methylation status NA

        Reads with a phred_fail_flag are considered. Also below-threshold
        phred scores contain information and may be used to adjust the
        phred score of the better read.
        """
        read_hash: DefaultDict[str, List[BSSeqPileupRead]] = defaultdict(list)
        for curr_read in motif_pileup.reads:
            if (curr_read.trimm_flag
                or curr_read.qc_fail_flag
                or curr_read.meth_status_flag == mflags.is_na):
                continue
            read_hash[curr_read.alignment.query_name].append(curr_read)

        for unused_query_name, reads in read_hash.items():
            if len(reads) == 2:
                self._process_overlap(reads[0], reads[1])


    @staticmethod
    def _process_overlap(read1: BSSeqPileupRead, read2: BSSeqPileupRead) -> None:
        """Decide what to do with overlapping reads

        If meth. calls from both reads match, choose the better read (mate1
        if equal phred). Otherwise, discard both reads. Discards are
        implemented via overlap_fail bit of qc_fail_flag.
        """

        if read1.meth_status_flag == read2.meth_status_flag:
            # typically, mapq filtering will be applied before overlap handling
            # so if we get to here, both reads pass mapq filtering
            # by choosing the read with the better phred score, we make sure
            # that the 'higher confidence' event is handed to the phred score
            # filter, which will typically be applied next
            # if this step is neglected and a read is selected at random, it may
            # be that the selected read is discarded as phred score failure
            # even though its mate would have passed the phred score filter
            if read1.baseq_at_pos >= read2.baseq_at_pos:
                read2.overlap_flag = 1
            else:
                read1.overlap_flag = 1
        else:
            # Better safe than sorry. Statistically valid usage of phred
            # scores to select high confidence calls from disagreeing overlaps
            # will be implemented in the future
            read1.qc_fail_flag |= qflags.overlap_fail
            read2.qc_fail_flag |= qflags.overlap_fail
