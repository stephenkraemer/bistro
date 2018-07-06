"""MotifPileupProcessors for QC filtering"""

from mqc.visitors import Visitor
from mqc.pileup.pileup import MotifPileup

import mqc.flag_and_index_values as mfl
qflag = mfl.qc_fail_flags
mflag = mfl.methylation_status_flags


from typing import Dict, Any
ConfigDict = Dict[str, Any]

class PhredFilter(Visitor):
    def __init__(self, config: ConfigDict) -> None:
        # will keep phred >= min_phred_score
        self.min_phred_score = config['basic_quality_filtering']['min_phred_score']
    def process(self, motif_pileup: MotifPileup) -> None:
        for curr_read in motif_pileup.reads:

            if (curr_read.qc_fail_flag
                # phred filtering is usually done after overlap handling,
                # because overlap handling can be used to adjust phred scores
                or curr_read.overlap_flag
                # trimming is done before overlap handling
                or curr_read.trimm_flag
                or curr_read.meth_status_flag == mflag.is_na):
                continue

            if curr_read.baseq_at_pos < self.min_phred_score:
                curr_read.qc_fail_flag |= qflag.phred_score_fail

class MapqFilter(Visitor):
    def __init__(self, config: ConfigDict) -> None:
        # will keep mapq >= min_mapq
        self.min_mapq = config['basic_quality_filtering']['min_mapq']
    def process(self, motif_pileup: MotifPileup) -> None:
        for curr_read in motif_pileup.reads:

            if (curr_read.qc_fail_flag
                # mapq filter should be applied before overlap calling
                or curr_read.meth_status_flag == mflag.is_na):
                continue

            if curr_read.alignment.mapping_quality < self.min_mapq:
                curr_read.qc_fail_flag |= qflag.mapq_fail

