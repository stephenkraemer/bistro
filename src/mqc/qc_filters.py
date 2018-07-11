"""MotifPileupProcessors for QC filtering"""

# TODO: combine PhredFilter and MapqFilter to avoid cost for double
# iteration

from mqc.visitors import Visitor
from mqc.pileup.pileup import MotifPileup

import mqc.flag_and_index_values as mfl
qflag = mfl.qc_fail_flags
mflag = mfl.methylation_status_flags


from typing import Dict, Any
ConfigDict = Dict[str, Any]

class PhredFilter(Visitor):
    """Tag BSSeqPileupRead with low phred score qcfailflag bit

    Reads are tagged as qcfail with the appropriate flag value if
    BSSeqPileupRead.baseq_at_pos < min_phred_score.

    Notes:
    - phred filtering is usually done after overlap handling,
      because overlap handling can be used to adjust phred scores
    """
    def __init__(self, config: ConfigDict) -> None:
        self.min_phred_score = config['basic_quality_filtering']['min_phred_score']
    def process(self, motif_pileup: MotifPileup) -> None:
        """Check all reads and flag if appropriate

        The qcfail flag is set independent of the presence of other
        failure flags: overlap_flag, trimming flag and other qcfail
        flags are not considered. This is by design. In different
        scenarios, overlap flags are not considered (mate-stratified
        methylation calling, M-bias stats), so reads must have
        phred score fail informations independent of the overlap
        flag. Similar situations arise for the other fail flags.
        """
        for curr_read in motif_pileup.reads:
            if curr_read.baseq_at_pos < self.min_phred_score:
                curr_read.phred_fail_flag = 1

class MapqFilter(Visitor):
    """Tag BSSeqPileupRead with mapq fail flag

    Adds mapq_fail flag to BSSeqPileupRead.qflag if
    BSSeqPileupRead mapping quality is < min_mapq parameter
    """
    def __init__(self, config: ConfigDict) -> None:
        self.min_mapq = config['basic_quality_filtering']['min_mapq']
    def process(self, motif_pileup: MotifPileup) -> None:
        for curr_read in motif_pileup.reads:
            if curr_read.alignment.mapping_quality < self.min_mapq:
                curr_read.qc_fail_flag |= qflag.mapq_fail

