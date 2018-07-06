"""Tests for QC fail filters

- phred filter
- mapq filter
"""
from typing import cast

from mqc.pileup.pileup import MotifPileup
from mqc.qc_filters import PhredFilter, MapqFilter

import mqc.flag_and_index_values as mfl
qflag = mfl.qc_fail_flags
mflags = mfl.methylation_status_flags


class MotifPileupStub:
    def __init__(self, reads):
        self.reads = reads


def test_phred_filter_tags_reads_with_below_threshold_phred_score_with_specific_qc_fail_flag() -> None:

    class PileupreadStub:
        def __init__(self, phred_score, qc_fail_flag=0):
            self.qc_fail_flag = qc_fail_flag
            self.baseq_at_pos = phred_score
            self.overlap_flag = 0
            self.meth_status_flag = mflags.is_methylated
            self.trimm_flag = 0
    pread = PileupreadStub

    motif_pileup = MotifPileupStub([pread(1), pread(10), pread(20), pread(30)])

    config = { 'basic_quality_filtering': {'min_phred_score': 20}}
    phred_filter = PhredFilter(config)

    phred_filter.process(cast(MotifPileup, motif_pileup))
    for curr_read in motif_pileup.reads:
        if curr_read.baseq_at_pos < 20:
            assert curr_read.qc_fail_flag & qflag.phred_score_fail > 0
        else:
            assert curr_read.qc_fail_flag & qflag.phred_score_fail == 0

def test_mapq_filter_tags_reads_with_below_threshold_mapq_with_specific_qc_fail_flag() -> None:

    class AlignedSegmentStub:
        def __init__(self, mapq):
            self.mapping_quality = mapq

    class PileupreadStub:
        def __init__(self, mapq, qc_fail_flag=0):
            self.alignment = AlignedSegmentStub(mapq)
            self.qc_fail_flag = qc_fail_flag
            self.meth_status_flag = mflags.is_unmethylated
    pread = PileupreadStub

    motif_pileup = MotifPileupStub([pread(1), pread(10), pread(20), pread(30)])

    config = {'basic_quality_filtering': {'min_mapq': 20}}
    mapq_filter = MapqFilter(config)

    mapq_filter.process(cast(MotifPileup, motif_pileup))

    for curr_read in motif_pileup.reads:
        if curr_read.alignment.mapping_quality >= 20:
            assert curr_read.qc_fail_flag & qflag.mapq_fail == 0
        else:
            assert curr_read.qc_fail_flag & qflag.mapq_fail > 0
